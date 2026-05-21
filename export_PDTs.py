"""
export_PDTs.py

Reads the PutativeKinaseSubstrates golden CSVs for each cell line,
filters out kinases with no substrates (n=0), and writes one Excel
file per cell line with columns:
    A = kinase
    B = n (number of PDTs)
    C = substrates (semicolon-separated)

Output files are saved in the same folder as this script.

Run from the project root:
    python export_PDTs.py
"""

import csv
import os
import openpyxl
from openpyxl.styles import Font

GOLDEN_DIR = os.path.join("tests", "golden")
CELL_LINES = ["MCF7", "HL60", "NTERA2"]

output_path = "PDTs.xlsx"
wb = openpyxl.Workbook()
wb.remove(wb.active)  # remove default empty sheet

for cell_line in CELL_LINES:
    input_path = os.path.join(GOLDEN_DIR, f"{cell_line}_PutativeKinaseSubstrates.csv")

    with open(input_path, "rt", encoding="utf-8") as f:
        reader = csv.reader(f)
        rows = list(reader)

    ws = wb.create_sheet(title=cell_line)

    # Header
    ws["A1"] = "kinase"
    ws["B1"] = "n"
    ws["C1"] = "substrates"
    for cell in [ws["A1"], ws["B1"], ws["C1"]]:
        cell.font = Font(bold=True)

    # Data: skip header row, skip kinases with n=0
    excel_row = 2
    for row in rows[1:]:
        if not row:
            continue
        kinase, n, substrates = row[0], row[1], row[2] if len(row) > 2 else ""
        if int(float(n)) == 0:
            continue
        ws.cell(row=excel_row, column=1).value = kinase
        ws.cell(row=excel_row, column=2).value = int(float(n))
        ws.cell(row=excel_row, column=3).value = substrates
        excel_row += 1

    print(f"{cell_line}: {excel_row - 2} kinases with PDTs")

wb.save(output_path)
print(f"Saved -> {output_path}")
