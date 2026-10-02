import argparse
import re
import numpy as np

complex_re = re.compile(r"\(\s*([+-]?\d+\.\d+E[+-]\d+)\s*,\s*([+-]?\d+\.\d+E[+-]\d+)\s*\)")
meta_re = re.compile(r"([A-Z_]+)=([0-9]+)")


def read_loop_output(path):
    data = {}
    current_loop = None
    current_meta = {}
    matrix_rows = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if line.startswith("Loop "):
                current_loop = int(line.split()[1])
                current_meta = {}
                matrix_rows = []

            elif line.startswith("N4="):
                current_meta = {k: int(v) for k, v in meta_re.findall(line)}

            elif line.startswith("("):
                vals = complex_re.findall(line)
                row = [float(re) + 1j * float(im) for re, im in vals]
                matrix_rows.append(row)

                if len(matrix_rows) == 2:
                    key = (
                        current_meta.get("N4"),
                        current_meta.get("IBL"),
                        current_loop,
                    )
                    data[key] = {
                        "meta": current_meta,
                        "A11": np.array(matrix_rows, dtype=np.complex128),
                    }

    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("old_output")
    parser.add_argument("new_output")
    parser.add_argument("--tol", type=float, default=5e-6)
    args = parser.parse_args()

    old = read_loop_output(args.old_output)
    new = read_loop_output(args.new_output)

    ok = True

    old_keys = set(old)
    new_keys = set(new)

    for key in sorted(old_keys - new_keys):
        print(f"Missing in new output: {key}")
        ok = False

    for key in sorted(new_keys - old_keys):
        print(f"Missing in old output: {key}")
        ok = False

    for key in sorted(old_keys & new_keys):
        old_A11 = old[key]["A11"]
        new_A11 = new[key]["A11"]

        diff = np.linalg.norm(old_A11 - new_A11)

        if diff > args.tol:
            print("----------------------------------------")
            print(f"WARNING: mismatch at N4={key[0]} IBL={key[1]} Loop={key[2]}")
            print(f"norm(old-new) = {diff:.8e}")
            print("old A11 =")
            print(old_A11)
            print("new A11 =")
            print(new_A11)
            ok = False

        old_site = old[key]["meta"].get("OUT_SITE")
        new_site = new[key]["meta"].get("OUT_SITE")

        if old_site is not None and new_site is not None and old_site != new_site:
            print("----------------------------------------")
            print(f"WARNING: out_site mismatch at N4={key[0]} IBL={key[1]} Loop={key[2]}")
            print(f"old OUT_SITE = {old_site}")
            print(f"new OUT_SITE = {new_site}")
            ok = False

    if ok:
        print("All common loop outputs agree.")


if __name__ == "__main__":
    main()
