import sys

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")

import pytest
import parse_alignment as pa

file = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/pytest/test_parse_samples.txt"

"""
P1:  A->T, index 29 n_replacement = 1
P2: coverage percent = 0.5, n_replacement = 0
P3:  G->R at 3, Y->V at 14 n_replacement = 2
P4: K->V, K->B, K->T at 3 n_replacement = 1
P5: K->V, K->B, K->T at 3 | G->T, G->H, G->A at 7 n_replacement = 2
P6: H->A, H->F, H->R at 22 | N->Q at 30, pcov = 27/36 = 3/4
"""

a1 = "MSKGPAIGIDLGTTYSCVGVFQHGKVEIIANDQGVD"
b1 = "MSKGPAIGIDLGTTYSCVGVFQHGKVEIITNDQGVD"

a2 = "MSKGPAIGIDLGTTYSCVGVFQHF"
b2 = "------------TTYSCVGVFQHF"

a3 = "MSKGPAIGIDLGTTYS"
b3 = "MSKRPAIGIDLGTTVS"

a4 = "MSKGPAIGIDLG"
b4 = "MS[VTR]GPAIGIDLG"

a5 = "MSKGPAIGIDLG"
b5 = "MS[VTR]GPAI[THA]IDLG"

a6 = "MSKGPAIGIDLGTTYSCVGVFQHGKVEIIANDQGVD"
b6 = "MSKGP---------YSCVGVFQ[AFR]GKVEIITQDQGNR"

a7 = "RHLQLAVRNDEELNKLLAGVTIAQGGVLPNIQVLLPKKTEKKQH"
b7 = "--LLLAVL[DN]DEEL[DN]KLLSGV[TDN]LAQG[GDN][VDN][LDN]P[DN][DN][DN][DN][DN][DN][DN][DN][DN][DN]EEKQA"

a8 = "LAVRNDEELNKLLAGVTIAQGGVLPNIQAVLLPKKTEKKQH"
b8 = "LAV[NLT][NLT][NLT]E[NLT]L[NLT]KLL[SNLT]G[VNLT][NLT][NLT][ANLT][QNLT]GGVLP[NLT][NLT][NLT][NLT][NLT][NLT]L[NLT][NLT][NLT][NLT][NLT]EKQA"


@pytest.mark.parametrize(
    "cases,n_replacements,changes,cov_percent",
    [
        ((a1, b1), 1, ("A->T",), 1),
        ((a2, b2), 0, [], 0.5),
        ((a3, b3), 2, ("G->R", "Y->V"), 1),
        ((a4, b4), 1, ("K->V", "K->R", "K->T"), 1),
        ((a5, b5), 2, ("K->V", "K->R", "K->T", "G->T", "G->H", "G->A"), 1),
        (
            (a6, b6),
            5,
            ("H->A", "H->F", "H->R", "A->T", "V->N", "D->R", "N->Q"),
            3 / 4,
        ),
    ],
)
def test_basic(cases, n_replacements, changes, cov_percent):
    _r, _m = pa.recordAlignment("p", cases[0], cases[1])
    assert _m["n_replacements"].item() == n_replacements
    assert _m["pcoverage_nmatch"].item() == cov_percent
    assert not changes == _r.empty
    if changes and not _r.empty:
        assert _r["change"].isin(changes).all()
        assert len(_r["change"]) == len(changes)


def test_file():
    m, r = pa.parseAlignmentFile(file)
    return m, r
