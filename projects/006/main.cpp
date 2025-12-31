#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "cpl_vsi.h"

#include <proj.h>

static void die(const std::string& msg, int code = 2) {
    std::cerr << msg << "\n";
    std::exit(code);
}

struct Args {
    std::string inPath;
    std::string outPath = "-";
    std::string dst = "EPSG:3857";
    std::string srcOverride;          // optional
    bool traditionalAxis = false;
    bool projCheck = false;
};

static Args parse_args(int argc, char** argv) {
    Args a;
    if (argc < 2) {
        std::cerr << "Usage: 006 <in> [out|-] --dst <SRS> [--src <SRS>] [--traditional-axis] [--proj-check]\n";
        std::exit(2);
    }
    a.inPath = argv[1];
    int i = 2;

    // 2番目が out（通常のパス or "-"）なら先に消費する
    if (i < argc) {
        const std::string cand = argv[i];
        if (cand == "-" || (!cand.empty() && cand[0] != '-')) {
            a.outPath = cand;
            i++;
        }
    }

    for (; i < argc; i++) {
        std::string k = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) die(std::string("Missing value for ") + name);
            return std::string(argv[++i]);
        };

        if (k == "--dst") a.dst = need("--dst");
        else if (k == "--src") a.srcOverride = need("--src");
        else if (k == "--traditional-axis") a.traditionalAxis = true;
        else if (k == "--proj-check") a.projCheck = true;
        else if (k == "--help") {
            std::cout << "Usage: 006 <in> [out|-] --dst <SRS> [--src <SRS>] [--traditional-axis] [--proj-check]\n";
            std::exit(0);
        } else {
            die("Unknown option: " + k);
        }
    }
    return a;
}

static void apply_axis_strategy(OGRSpatialReference& srs, bool traditionalAxis) {
    if (traditionalAxis) {
        // lat/lon 定義でもデータは lon/lat 扱いに寄せる :contentReference[oaicite:6]{index=6}
        srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    }
}

static void proj_sanity_check_one_point(
    const std::string& srcDef,
    const std::string& dstDef,
    bool normalizeForViz,
    double x_in, double y_in,
    double x_gdal, double y_gdal
) {
    PJ_CONTEXT* ctx = proj_context_create();
    if (!ctx) return;

    PJ* P = proj_create_crs_to_crs(ctx, srcDef.c_str(), dstDef.c_str(), nullptr);
    if (!P) {
        proj_context_destroy(ctx);
        return;
    }

    PJ* N = P;
    if (normalizeForViz) {
        // lon/lat順で投入しやすくする（可視化用途の正規化）:contentReference[oaicite:7]{index=7}
        N = proj_normalize_for_visualization(ctx, P);
        proj_destroy(P);
        if (!N) {
            proj_context_destroy(ctx);
            return;
        }
    }

    PJ_COORD c = proj_coord(x_in, y_in, 0, 0);
    PJ_COORD out = proj_trans(N, PJ_FWD, c);

    const double dx = out.xy.x - x_gdal;
    const double dy = out.xy.y - y_gdal;

    std::cerr << "[proj-check] input=(" << x_in << "," << y_in << ") "
              << "PROJ=(" << out.xy.x << "," << out.xy.y << ") "
              << "GDAL=(" << x_gdal << "," << y_gdal << ") "
              << "diff=(" << dx << "," << dy << ")\n";

    proj_destroy(N);
    proj_context_destroy(ctx);
}

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);

    GDALAllRegister();

    GDALDataset* inDs = static_cast<GDALDataset*>(
        GDALOpenEx(args.inPath.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr)
    );
    if (!inDs) die("Failed to open input: " + args.inPath);

    // 出力：ファイル or stdout(/vsimem)
    const bool toStdout = (args.outPath == "-" || args.outPath.empty());
    const std::string realOut = toStdout ? "/vsimem/out_006.geojson" : args.outPath;

    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GeoJSON");
    if (!drv) die("GeoJSON driver not found.");

    GDALDataset* outDs = drv->Create(realOut.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!outDs) die("Failed to create output: " + realOut);

    // dst SRS
    OGRSpatialReference dstSRS;
    if (dstSRS.SetFromUserInput(args.dst.c_str()) != OGRERR_NONE) {
        die("Failed to parse --dst SRS: " + args.dst);
    }
    apply_axis_strategy(dstSRS, args.traditionalAxis);

    const int layerCount = inDs->GetLayerCount();
    for (int li = 0; li < layerCount; li++) {
        OGRLayer* inLy = inDs->GetLayer(li);
        if (!inLy) continue;

        // src SRS（無ければ --src を要求）
        OGRSpatialReference srcSRS;
        const OGRSpatialReference* srcFromLayer = inLy->GetSpatialRef();
        if (srcFromLayer) {
            srcSRS = *srcFromLayer; // copy
        } else if (!args.srcOverride.empty()) {
            if (srcSRS.SetFromUserInput(args.srcOverride.c_str()) != OGRERR_NONE) {
                die("Failed to parse --src SRS: " + args.srcOverride);
            }
        } else {
            die("Input layer has no SRS. Provide --src <SRS> (e.g., EPSG:4326).");
        }
        apply_axis_strategy(srcSRS, args.traditionalAxis);

        // 変換器（OGRCreateCoordinateTransformation の世界）:contentReference[oaicite:8]{index=8}
        OGRCoordinateTransformation* ct = OGRCreateCoordinateTransformation(&srcSRS, &dstSRS);
        if (!ct) die("Failed to create coordinate transformation (src -> dst).");

        OGRFeatureDefn* inDefn = inLy->GetLayerDefn();
        const OGRwkbGeometryType outGeomType = inDefn->GetGeomType();

        OGRLayer* outLy = outDs->CreateLayer(inLy->GetName(), &dstSRS, outGeomType, nullptr);
        if (!outLy) die("Failed to CreateLayer().");

        // フィールドをコピー
        for (int fi = 0; fi < inDefn->GetFieldCount(); fi++) {
            OGRFieldDefn* fdef = inDefn->GetFieldDefn(fi);
            if (outLy->CreateField(fdef) != OGRERR_NONE) {
                die("Failed to CreateField() while copying schema.");
            }
        }

        // 1点だけPROJとGDALの結果を突き合わせ（任意）
        bool checked = false;

        inLy->ResetReading();
        while (OGRFeature* f = inLy->GetNextFeature()) {
            OGRFeature* outF = OGRFeature::CreateFeature(outLy->GetLayerDefn());
            outF->SetFrom(f, TRUE);

            const OGRGeometry* g = f->GetGeometryRef();
            if (g) {
                OGRGeometry* gg = g->clone();
                if (gg->transform(ct) != OGRERR_NONE) {
                    OGRGeometryFactory::destroyGeometry(gg);
                    OGRFeature::DestroyFeature(outF);
                    OGRFeature::DestroyFeature(f);
                    OCTDestroyCoordinateTransformation(ct);
                    die("Geometry transform failed (check SRS/axis order).");
                }
                outF->SetGeometryDirectly(gg); // ownership移譲 :contentReference[oaicite:9]{index=9}

                if (args.projCheck && !checked) {
                    // geometryの最初の座標（Point化できない場合は envelope 中心）で軽く比較
                    OGREnvelope env;
                    gg->getEnvelope(&env);
                    const double x_in = (env.MinX + env.MaxX) * 0.5;
                    const double y_in = (env.MinY + env.MaxY) * 0.5;

                    // GDALでその点を変換
                    double x_gdal = x_in, y_gdal = y_in;
                    ct->Transform(1, &x_gdal, &y_gdal);

                    // PROJで同じ点を変換（normalize_for_visualization を使うと lon/lat順で扱いやすい）
                    std::string srcDefWKT, dstDefWKT;
                    char* srcWkt = nullptr;
                    char* dstWkt = nullptr;
                    srcSRS.exportToWkt(&srcWkt);
                    dstSRS.exportToWkt(&dstWkt);
                    if (srcWkt && dstWkt) {
                        proj_sanity_check_one_point(srcWkt, dstWkt, true, x_in, y_in, x_gdal, y_gdal);
                    }
                    CPLFree(srcWkt);
                    CPLFree(dstWkt);

                    checked = true;
                }
            }

            if (outLy->CreateFeature(outF) != OGRERR_NONE) {
                OGRFeature::DestroyFeature(outF);
                OGRFeature::DestroyFeature(f);
                OCTDestroyCoordinateTransformation(ct);
                die("Failed to CreateFeature().");
            }

            OGRFeature::DestroyFeature(outF);
            OGRFeature::DestroyFeature(f);
        }

        OCTDestroyCoordinateTransformation(ct);
    }

    GDALClose(inDs);
    GDALClose(outDs);

    // stdoutモード：/vsimem から出力を取り出す（C++組み込み向き）:contentReference[oaicite:10]{index=10}
    if (toStdout) {
        vsi_l_offset nSize = 0;
        GByte* p = VSIGetMemFileBuffer(realOut.c_str(), &nSize, TRUE);
        if (!p) die("VSIGetMemFileBuffer() failed.");
        std::cout.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(nSize));
        CPLFree(p);
    }

    return 0;
}
