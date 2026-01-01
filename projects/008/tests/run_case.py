#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys

def approx_equal(a: float, b: float, eps: float) -> bool:
    return abs(a - b) <= eps

def require(cond: bool, msg: str):
    if not cond:
        raise AssertionError(msg)

def wkt_prefix(wkt: str) -> str:
    if not wkt:
        return ""
    return wkt.strip().split("(")[0].strip().upper()

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--eps", type=float, default=1e-9)
    ap.add_argument("--emit-wkt", action="store_true")
    args = ap.parse_args()

    with open(args.case, "r", encoding="utf-8") as f:
        case = json.load(f)

    op = case["op"]
    a = case["a"]
    b = case["b"]
    expect = case.get("expect", {})
    name = case.get("name", args.case)

    cmd = [args.exe, "--op", op, "--a", a, "--b", b]
    if args.emit_wkt:
        cmd.append("--emit-wkt")

    p = subprocess.run(cmd, text=True, capture_output=True)

    # 008は「片側失敗」で非0終了もありうる設計にしているので、
    # 期待終了コードを指定できるようにしておく（未指定は0）
    expect_rc = int(expect.get("exit_code", 0))
    if p.returncode != expect_rc:
        sys.stderr.write(f"[{name}] unexpected exit code: got={p.returncode} expect={expect_rc}\n")
        sys.stderr.write(p.stdout + "\n")
        sys.stderr.write(p.stderr + "\n")
        return 1

    try:
        out = json.loads(p.stdout)
    except Exception as e:
        sys.stderr.write(f"[{name}] JSON parse error: {e}\n")
        sys.stderr.write(p.stdout + "\n")
        return 1

    def check_side(side_key: str):
        side = out.get(side_key)
        require(side is not None, f"[{name}] missing key: {side_key}")

        # ok/valid は指定があるときだけ厳密にチェック
        if "ok" in expect:
            require(side.get("ok") == expect["ok"], f"[{name}] {side_key}.ok mismatch")
        else:
            require(side.get("ok") is True, f"[{name}] {side_key}.ok is not true")

        if "valid" in expect:
            require(side.get("valid") == expect["valid"], f"[{name}] {side_key}.valid mismatch")

        if "area" in expect:
            got = float(side.get("area"))
            exp = float(expect["area"])
            require(approx_equal(got, exp, args.eps),
                    f"[{name}] {side_key}.area mismatch: got={got} expect={exp} eps={args.eps}")

        if "perimeter_all_rings" in expect:
            got = float(side.get("perimeter_all_rings"))
            exp = float(expect["perimeter_all_rings"])
            require(approx_equal(got, exp, args.eps),
                    f"[{name}] {side_key}.perimeter_all_rings mismatch: got={got} expect={exp} eps={args.eps}")

        if "vertices" in expect:
            got = int(side.get("vertices"))
            exp = int(expect["vertices"])
            require(got == exp, f"[{name}] {side_key}.vertices mismatch: got={got} expect={exp}")

    check_side("geos")
    check_side("boost_geometry")

    # ---- WKT差分の検査（emit-wkt時のみ）----
    if args.emit_wkt and expect.get("wkt", {}) is not None:
        geos_wkt = out.get("geos", {}).get("wkt", "")
        boost_wkt = out.get("boost_geometry", {}).get("wkt", "")

        if expect["wkt"].get("equal") is not None:
            want_equal = bool(expect["wkt"]["equal"])
            require((geos_wkt == boost_wkt) == want_equal,
                    f"[{name}] wkt equality mismatch: want_equal={want_equal}")

        if expect["wkt"].get("len_diff_gt") is not None:
            lim = int(expect["wkt"]["len_diff_gt"])
            require(abs(len(geos_wkt) - len(boost_wkt)) > lim,
                    f"[{name}] wkt_len_diff not > {lim}")

        if expect["wkt"].get("geos_prefix") is not None:
            require(wkt_prefix(geos_wkt) == expect["wkt"]["geos_prefix"],
                    f"[{name}] geos wkt prefix mismatch")

        if expect["wkt"].get("boost_prefix") is not None:
            require(wkt_prefix(boost_wkt) == expect["wkt"]["boost_prefix"],
                    f"[{name}] boost wkt prefix mismatch")

    # diff側（指定があれば）
    diff = out.get("diff", {})
    if expect.get("diff_abs_area_le") is not None and "abs_area_diff" in diff:
        require(float(diff["abs_area_diff"]) <= float(expect["diff_abs_area_le"]),
                f"[{name}] diff.abs_area_diff too large")

    if expect.get("diff_abs_perimeter_le") is not None and "abs_perimeter_diff" in diff:
        require(float(diff["abs_perimeter_diff"]) <= float(expect["diff_abs_perimeter_le"]),
                f"[{name}] diff.abs_perimeter_diff too large")

    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as e:
        sys.stderr.write(str(e) + "\n")
        raise SystemExit(1)
