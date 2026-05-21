"""
regenerate_golden.py

Runs the EBDT pipeline with the paper's default parameters on the three
real datasets and overwrites the golden CSV files in tests/golden/.

Paper default parameters:
    inhibitionThreshold  = 0.5
    ratioThreshold       = 0.5
    probabilityThreshold = 0.75

Run from the project root:
    python tests/regenerate_golden.py
"""

import os
import sys
import csv
import shutil
import tempfile

# Make sure the project root is in the path
CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, CODE_DIR)

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "golden")

INHIBITION_THRESHOLD  = 0.5
RATIO_THRESHOLD       = 0.5
PROBABILITY_THRESHOLD = 0.75

SHEETS_TO_SAVE = [
    "PutativeKinaseSubstrates",
    "corrPPwithKinases",
    "nodes.edges",
    "ProbOfBeingKinaseSubs",
    "ratioSigPPoverSignInhSpeci",
    "tableEdgeSubs",
]

import openpyxl
from ebdtFunctions import GetExpectancyOfBeingDownstreamTarget


def save_sheet_as_csv(wb, sheet_name, csv_path):
    ws = wb[sheet_name]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        for row in ws.iter_rows(values_only=True):
            writer.writerow([v if v is not None else "" for v in row])


for cell_line in ["MCF7", "HL60", "NTERA2"]:
    print(f"\nProcessing {cell_line}...")

    run_dir = tempfile.mkdtemp(prefix=f"{cell_line}_golden_")
    try:
        xlsm_src = os.path.join(CODE_DIR, f"{cell_line}.xlsm")
        xlsm_dst = os.path.join(run_dir, f"{cell_line}.xlsm")
        shutil.copy(xlsm_src, xlsm_dst)
        shutil.copytree(
            os.path.join(CODE_DIR, "requiredData"),
            os.path.join(run_dir, "requiredData")
        )

        old_cwd = os.getcwd()
        try:
            os.chdir(run_dir)
            GetExpectancyOfBeingDownstreamTarget(
                inhibitionThreshold  = INHIBITION_THRESHOLD,
                ratioThreshold       = RATIO_THRESHOLD,
                probabilityThreshold = PROBABILITY_THRESHOLD,
                cellLinesFiles       = [f"{cell_line}.xlsm"],
            )
        finally:
            os.chdir(old_cwd)

        wb = openpyxl.load_workbook(xlsm_dst, read_only=True, keep_vba=False)
        for sheet_name in SHEETS_TO_SAVE:
            csv_path = os.path.join(GOLDEN_DIR, f"{cell_line}_{sheet_name}.csv")
            save_sheet_as_csv(wb, sheet_name, csv_path)
            print(f"  Saved {cell_line}_{sheet_name}.csv")
        wb.close()

    finally:
        shutil.rmtree(run_dir, ignore_errors=True)

print("\nAll golden files regenerated.")
