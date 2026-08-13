"""
validate_against_original.py
-----------------------------
Compares Output/<CaseName>/results.csv (produced by run_case.py) against
a reference results file, at matching reported time steps. Useful any
time you have an independent reference (e.g. an existing spreadsheet
model, a prior tool, or hand calculations) to check a case against.

Run AFTER run_case.py, from anywhere, e.g.:
    python Module/validate_against_original.py Template
(defaults to "Template" if no case name is given)
"""

import csv
import os
import sys

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(MODULE_DIR)
ORIGINAL_CSV = "/path/to/your/reference_results.csv"  # <-- EDIT this to your reference file

CASE_NAME = sys.argv[1] if len(sys.argv) > 1 else "Template"
OUR_CSV = os.path.join(ROOT, "Output", CASE_NAME, "Baseline", "results.csv")


def read_original():
    """NOTE: the row/column indices below (reader[14:63], row[1], row[2],
    etc.) were specific to one particular reference spreadsheet's export
    layout. Adapt these to match whatever reference file you're actually
    comparing against -- there's no way to make this generic without
    knowing your reference file's own column layout."""
    with open(ORIGINAL_CSV, encoding="utf-8-sig") as f:
        reader = list(csv.reader(f))
    rows = reader[14:63]
    data = []
    for row in rows:
        data.append({
            "t_h": float(row[1]), "H": float(row[2]), "Qin": float(row[3]),
            "Qwithdrawal": float(row[4]), "Qs_total": float(row[11]),
        })
    return data


def read_ours():
    with open(OUR_CSV, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = [row for row in reader]
    return rows


def main():
    orig = read_original()
    ours = read_ours()
    ours_by_th = {round(float(r["time_h"]), 2): r for r in ours}

    print(f"{'t[h]':>6} {'H_orig':>8} {'H_ours':>8} {'dH':>7}  {'Qtot_orig':>10} {'Qtot_ours':>10}")
    max_dH = 0.0
    for row in orig:
        t = round(row["t_h"], 2)
        match = ours_by_th.get(t)
        if match is None:
            continue
        H_ours = float(match["H_masl"])
        dH = H_ours - row["H"]
        max_dH = max(max_dH, abs(dH))
        print(f"{t:6.1f} {row['H']:8.2f} {H_ours:8.2f} {dH:7.3f}  "
              f"{row['Qs_total']:10.2f} {float(match['Qs_total']):10.2f}")

    print(f"\nMax |H_ours - H_original| over full run: {max_dH:.3f} m")

    peak_orig_H = max(r["H"] for r in orig)
    peak_orig_t = [r["t_h"] for r in orig if r["H"] == peak_orig_H][0]
    peak_ours = max(ours, key=lambda r: float(r["H_masl"]))
    print(f"Original peak: {peak_orig_H:.2f} m at t={peak_orig_t:.2f} h")
    print(f"Kernel peak:   {float(peak_ours['H_masl']):.2f} m at t={float(peak_ours['time_h']):.2f} h")


if __name__ == "__main__":
    main()