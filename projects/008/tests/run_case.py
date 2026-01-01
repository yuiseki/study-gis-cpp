#!/usr/bin/env python3
import argparse
import json
import math
import subprocess
import sys


def approx_equal(a: float, b: float, eps: float) -> bool:
    return abs(a - b) <= eps


def require(cond: bool, msg: str):
    if not cond:
        raise AssertionError(msg)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--eps", type=float, default=1e-9)
    args = ap.parse_args()

    with open(args.case, "r", encoding="utf-8") as f:
        case = json.load(f)

    op = case["op"]
    a = case["a"]
    b = case["b"]
    expect = case.get("expect", {})
    name = case.get("name", args.case)

    cmd = [args.exe, "--op", op, "--a", a, "--b", b]
    p = subprocess.run(cmd, text=True, capture_output=True)

    if p.returncode != 0:
        # 008は片側失敗でreturn 1になる設計なので、ここで即失敗にします
        sys.stderr.write(f"[{name}] command failed: {' '.join(cmd)}\n")
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

        require(side.get("ok") is True, f"[{name}] {side_key}.ok is not true")
        if "valid" in expect:
            require(side.get("valid") == expect["valid"],
                    f"[{name}] {side_key}.valid mismatch: got={side.get('valid')} expect={expect['valid']}")

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
            require(got == exp,
                    f"[{name}] {side_key}.vertices mismatch: got={got} expect={exp}")

    # 008の出力キーに合わせる
    check_side("geos")
    check_side("boost_geometry")

    # diff側の検査（指定があれば）
    diff = out.get("diff", {})
    if expect.get("diff_abs_area_le") is not None and "abs_area_diff" in diff:
        got = float(diff["abs_area_diff"])
        lim = float(expect["diff_abs_area_le"])
        require(got <= lim, f"[{name}] diff.abs_area_diff too large: got={got} limit={lim}")

    if expect.get("diff_abs_perimeter_le") is not None and "abs_perimeter_diff" in diff:
        got = float(diff["abs_perimeter_diff"])
        lim = float(expect["diff_abs_perimeter_le"])
        require(got <= lim, f"[{name}] diff.abs_perimeter_diff too large: got={got} limit={lim}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as e:
        sys.stderr.write(str(e) + "\n")
        raise SystemExit(1)
