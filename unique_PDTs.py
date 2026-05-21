"""
unique_PDTs.py

Reads PDTs.xlsx and counts unique PDTs per cell line and across all three.

Run from the project root:
    python unique_PDTs.py
"""

import openpyxl

wb = openpyxl.load_workbook("PDTs.xlsx", read_only=True)

pdts_per_line = {}

for cell_line in ["MCF7", "HL60", "NTERA2"]:
    ws = wb[cell_line]
    unique_pdts = set()
    for row in ws.iter_rows(min_row=2, values_only=True):
        substrates = row[2]
        if substrates:
            for site in substrates.split(";"):
                site = site.strip()
                if site:
                    unique_pdts.add(site)
    pdts_per_line[cell_line] = unique_pdts
    print(f"{cell_line}: {len(unique_pdts)} unique PDTs")

all_pdts = pdts_per_line["MCF7"] | pdts_per_line["HL60"] | pdts_per_line["NTERA2"]
print(f"\nAll 3 cell lines combined: {len(all_pdts)} unique PDTs")
