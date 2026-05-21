"""
test_golden.py - regression tests (golden files)

run the full pipeline with the real data (MCF7.xlsm, HL60.xlsm and NTERA2.xlsm)
and verify results are identical to the stored reference results in tests/golden/

These tests are useful to:
  - They guarantee that the code will produce the same results on any machine or future version
  - If someone modifies the algorithm, these test will fail an alert you

*Note:These tests are slow bacuse they process bigger datsets
Run:
    pytest tests/python/test_golden.py -v

Or exclude in quick runs:
    pytest tests/python/ -v -m "not slow"
"""

import os
import csv
import math
import shutil
import pytest
import openpyxl

from ebdtFunctions import GetExpectancyOfBeingDownstreamTarget


# ==== PARAMETERS OF THE PIPELINE (must match with those used when generating golden) ====

INHIBITION_THRESHOLD   = 0.5
RATIO_THRESHOLD        = 0.5
PROBABILITY_THRESHOLD  = 0.75


# ==== FIXTURE: EXECUTED REAL PIPELINE ====

@pytest.fixture(scope="module")
def golden_results(tmp_path_factory):
    """
    Executes real pipeline with real MCF7, HL60 and NTERA2.
    returns dict: {'MCF7': workbook, 'HL60': workbook, 'NTERA2': workbook}
    """
    from conftest import CODE_DIR, GOLDEN_DIR
    run_dir = str(tmp_path_factory.mktemp("golden_run"))

    for cell_line in ["MCF7", "HL60", "NTERA2"]:
        shutil.copy(
            os.path.join(CODE_DIR, f"{cell_line}.xlsm"),
            os.path.join(run_dir, f"{cell_line}.xlsm")
        )
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
            cellLinesFiles       = ["MCF7.xlsm", "HL60.xlsm", "NTERA2.xlsm"],
        )
    finally:
        os.chdir(old_cwd)

    return {
        cell_line: openpyxl.load_workbook(
            os.path.join(run_dir, f"{cell_line}.xlsm"),
            read_only=True,
            keep_vba=False
        )
        for cell_line in ["MCF7", "HL60", "NTERA2"]
    }


# ==== HELPER: READ GOLDEN CSV ====

def _read_golden(filename):
    from conftest import GOLDEN_DIR
    golden_dir = GOLDEN_DIR
    #reads a golden file CSV and returns a list of lists
    path = os.path.join(golden_dir, filename)
    with open(path, encoding="utf-8") as f:
        return list(csv.reader(f))


def _sheet_to_rows(workbook, sheet_name):
    #Extracts all the rows of a sheet as a list of lists
    ws = workbook[sheet_name]
    return [list(row) for row in ws.iter_rows(values_only=True)]


# ==== TESTS: PutativeKinaseSubstrates ====

@pytest.mark.slow
@pytest.mark.parametrize("cell_line", ["MCF7", "HL60", "NTERA2"])
def test_golden_numero_quinasas_con_sustratos(cell_line, golden_results):
    #the number of kinases with at least one putative substrate must be identical to the one in the golden file
    
    golden = _read_golden(f"{cell_line}_PutativeKinaseSubstrates.csv")
    actual_wb = golden_results[cell_line]
    actual = _sheet_to_rows(actual_wb, "PutativeKinaseSubstrates")

    # Counts kinases with n>0
    golden_con_sust = sum(1 for r in golden[1:] if r[1] and int(float(r[1])) > 0)
    actual_con_sust = sum(1 for r in actual[1:] if r[1] is not None and int(float(r[1])) > 0)

    assert actual_con_sust == golden_con_sust, (
        f"{cell_line}: {actual_con_sust} kinases with substrates "
        f"(e {golden_con_sust})"
    )


@pytest.mark.slow
@pytest.mark.parametrize("cell_line", ["MCF7", "HL60", "NTERA2"])
def test_golden_sustratos_de_quinasa_top(cell_line, golden_results):
    """
    Substrates of the kinase with more substrates must be identical
    to the ones in the golden file (exact reproducibility test).
    """
    golden = _read_golden(f"{cell_line}_PutativeKinaseSubstrates.csv")
    actual_wb = golden_results[cell_line]
    actual = _sheet_to_rows(actual_wb, "PutativeKinaseSubstrates")

    # finds the kinase with more substrates in the golden
    golden_data = [(r[0], int(float(r[1])), r[2]) for r in golden[1:] if r[1]]
    golden_data.sort(key=lambda x: -x[1])
    top_kinase, top_n, top_subs_golden = golden_data[0]

    # search this kinase in the current result
    actual_row = next((r for r in actual[1:] if r[0] == top_kinase), None)
    assert actual_row is not None, f"Not found {top_kinase} in result"

    assert int(float(actual_row[1])) == top_n, (
        f"{cell_line} {top_kinase}: {actual_row[1]} substrates ({top_n})"
    )


# ==== TESTS: nodes.edges ====

@pytest.mark.slow
@pytest.mark.parametrize("cell_line", ["MCF7", "HL60", "NTERA2"])
def test_golden_numero_pares_kinasas(cell_line, golden_results):
    """
    the number of kinase pairs with more shared substrates must be
    identical to the one in the golden file.
    """
    golden = _read_golden(f"{cell_line}_nodes.edges.csv")
    actual_wb = golden_results[cell_line]
    actual = _sheet_to_rows(actual_wb, "nodes.edges")

    # -1 for the header 
    golden_pares = len([r for r in golden[1:] if r[0]])
    actual_pares = len([r for r in actual[1:] if r[0] is not None])

    assert actual_pares == golden_pares, (
        f"{cell_line}: {actual_pares} pairs ({golden_pares})"
    )


# ==== TESTS: corrPPwithKinases ====

@pytest.mark.slow
@pytest.mark.parametrize("cell_line", ["MCF7", "HL60", "NTERA2"])
def test_golden_correlaciones_primera_fila(cell_line, golden_results):
    """
    Correlations of the first row of psites must be
    numerically equal to the ones in the golden file (tolerance 1e-6).
    """
    golden = _read_golden(f"{cell_line}_corrPPwithKinases.csv")
    actual_wb = golden_results[cell_line]
    actual = _sheet_to_rows(actual_wb, "corrPPwithKinases")

    golden_row = golden[1]
    actual_row = actual[1]

    assert len(golden_row) == len(actual_row), (
        f"{cell_line}: number of different columns: {len(actual_row)} vs {len(golden_row)}"
    )

    for col_idx, (g_val, a_val) in enumerate(zip(golden_row[1:], actual_row[1:]), start=1):
        if g_val == "" or g_val is None:
            assert a_val is None or a_val == ""
        else:
            g_float = float(g_val) if g_val else 0.0
            a_float = float(a_val) if a_val is not None else 0.0
            assert math.isclose(g_float, a_float, rel_tol=1e-6, abs_tol=1e-9), (
                f"{cell_line} col {col_idx}: {a_float} != {g_float} (golden)"
            )


@pytest.mark.slow
def test_golden_MCF7_correlaciones_primera_fila(golden_results):
    """
    Correlations of the first row of psites (MCF7) must be
    numerically equal to the ones in the golden file (tolerance 1e-6).
    """
    golden = _read_golden("MCF7_corrPPwithKinases.csv")
    actual_wb = golden_results["MCF7"]
    actual = _sheet_to_rows(actual_wb, "corrPPwithKinases")

    # Comparare row 2 (first psite, index 1 in the list)
    golden_row = golden[1]
    actual_row = actual[1]

    assert len(golden_row) == len(actual_row), (
        f"Number of different columns: {len(actual_row)} vs {len(golden_row)}"
    )

    for col_idx, (g_val, a_val) in enumerate(zip(golden_row[1:], actual_row[1:]), start=1):
        if g_val == "" or g_val is None:
            assert a_val is None or a_val == ""
        else:
            g_float = float(g_val) if g_val else 0.0
            a_float = float(a_val) if a_val is not None else 0.0
            assert math.isclose(g_float, a_float, rel_tol=1e-6, abs_tol=1e-9), (
                f"Col {col_idx}: {a_float} != {g_float} (golden)"
            )
