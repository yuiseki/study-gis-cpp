#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys

def require(cond: bool, msg: str):
    if not cond:
        raise AssertionError(msg)

def get(obj, path, default=None):
    cur = obj
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--case", required=True)
    args = ap.parse_args()

    with open(args.case, "r", encoding="utf-8") as f:
        case = json.load(f)

    name = case.get("name", args.case)
    op = case["op"]
    a = case["a"]
    b = case["b"]
    run = case.get("run", {})
    expect = case.get("expect", {})

    samples = int(run.get("samples", 600))
    seed = int(run.get("seed", 1))
    expand_deg = float(run.get("expand_deg", 1.0))

    cmd = [
        args.exe,
        "--op", op,
        "--a", a,
        "--b", b,
        "--samples", str(samples),
        "--seed", str(seed),
        "--expand-deg", str(expand_deg),
        "--no-emit-loops",
    ]

    p = subprocess.run(cmd, text=True, capture_output=True)

    exp_rc = int(expect.get("exit_code", 0))
    if p.returncode != exp_rc:
        sys.stderr.write(f"[{name}] unexpected exit code: got={p.returncode} expect={exp_rc}\n")
        sys.stderr.write("CMD: " + " ".join(cmd) + "\n")
        sys.stderr.write(p.stdout + "\n")
        sys.stderr.write(p.stderr + "\n")
        return 1

    try:
        out = json.loads(p.stdout)
    except Exception as e:
        sys.stderr.write(f"[{name}] JSON parse error: {e}\n")
        sys.stderr.write(p.stdout + "\n")
        return 1

    # 基本 sanity
    require(get(out, ["s2", "ok"]) is True, f"[{name}] s2.ok != true")
    require(get(out, ["cgal_nef_s2", "ok"]) is True, f"[{name}] cgal_nef_s2.ok != true")
    if expect.get("s2_valid_true", True):
        require(get(out, ["s2", "valid"]) is True, f"[{name}] s2.valid != true")

    mismatch_count = int(get(out, ["comparison", "containment_mismatch", "count"], -1))
    mismatch_rate = float(get(out, ["comparison", "containment_mismatch", "rate"], -1.0))
    require(mismatch_count >= 0 and mismatch_rate >= 0.0, f"[{name}] missing mismatch metrics")

    # 期待値チェック
    if "mismatch_count_eq" in expect:
        require(mismatch_count == int(expect["mismatch_count_eq"]),
                f"[{name}] mismatch_count mismatch: got={mismatch_count} expect={expect['mismatch_count_eq']}")

    if "mismatch_count_le" in expect:
        require(mismatch_count <= int(expect["mismatch_count_le"]),
                f"[{name}] mismatch_count too large: got={mismatch_count} limit={expect['mismatch_count_le']}")

    if "mismatch_rate_le" in expect:
        require(mismatch_rate <= float(expect["mismatch_rate_le"]),
                f"[{name}] mismatch_rate too large: got={mismatch_rate} limit={expect['mismatch_rate_le']}")

    if "mismatch_rate_ge" in expect:
        require(mismatch_rate >= float(expect["mismatch_rate_ge"]),
                f"[{name}] mismatch_rate too small: got={mismatch_rate} min={expect['mismatch_rate_ge']}")

    # 参考：S2の真の球面面積が出ていること（値自体はケース固定でなくても良い）
    s2_area_sr = float(get(out, ["s2", "area_sr"], float("nan")))
    require(s2_area_sr == s2_area_sr, f"[{name}] s2.area_sr is NaN")

    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as e:
        sys.stderr.write(str(e) + "\n")
        raise SystemExit(1)
