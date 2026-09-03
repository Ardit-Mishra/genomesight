"""API-level tests for the capabilities ported from the Streamlit app.

The ported modules keep their own unit tests (test_orf.py). These cover the
HTTP surface: that the endpoints exist, accept the same inputs the rest of the
API accepts, and reject bad input rather than returning a plausible wrong answer.
"""
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

# ATG + 25x GCT + TAA = 81 nt, a single unambiguous ORF on the plus strand.
ORF_SEQ = "ATG" + "GCT" * 25 + "TAA"


def test_orfs_finds_the_known_orf():
    r = client.post("/api/orfs", json={"sequence": ORF_SEQ, "min_length": 30})
    assert r.status_code == 200
    data = r.json()
    assert data["total"] == 1
    orf = data["orfs"][0]
    assert orf["length_nt"] == 81
    assert orf["start_codon"] == "ATG"
    assert orf["stop_codon"] == "TAA"
    assert orf["strand"] == "+"
    assert orf["protein"].startswith("M")
    assert orf["protein"].endswith("*")


def test_orfs_min_length_actually_filters():
    """A threshold above the only ORF must return nothing, not the ORF anyway."""
    r = client.post("/api/orfs", json={"sequence": ORF_SEQ, "min_length": 200})
    assert r.status_code == 200
    assert r.json()["total"] == 0


def test_orfs_reports_the_coordinate_frame_of_reference():
    """Minus-strand coordinates are not original-strand coordinates; say so."""
    r = client.post("/api/orfs", json={"sequence": ORF_SEQ, "min_length": 30})
    assert "reverse-complement" in r.json()["coordinate_note"]


def test_motif_search_finds_both_occurrences():
    r = client.post("/api/motifs", json={"sequence": "GAATTCATGCGAATTC", "pattern": "GAATTC"})
    assert r.status_code == 200
    data = r.json()
    assert data["total"] == 2
    assert [m["start_1based"] for m in data["matches"]] == [1, 11]


def test_motif_search_expands_iupac_ambiguity():
    """R must match A or G, and the expansion must be inspectable in the response."""
    r = client.post("/api/motifs", json={"sequence": "ATGAATGGATGC", "pattern": "ATGR"})
    data = r.json()
    assert data["regex"] == "ATG[AG]"
    assert sorted(m["matched_sequence"] for m in data["matches"]) == ["ATGA", "ATGG"]


def test_restriction_site_search():
    r = client.post("/api/restriction-sites", json={"sequence": "GAATTCAAAGAATTC", "enzyme": "EcoRI"})
    assert r.status_code == 200
    assert r.json()["total"] == 2


def test_unknown_enzyme_is_rejected_not_guessed():
    r = client.post("/api/restriction-sites", json={"sequence": "GAATTC", "enzyme": "NotAnEnzyme"})
    assert r.status_code == 422


def test_enzyme_list_is_served():
    r = client.get("/api/enzymes")
    assert r.status_code == 200
    enzymes = r.json()["enzymes"]
    assert enzymes["EcoRI"] == "GAATTC"


def test_endpoints_accept_bare_sequence_and_fasta_alike():
    """Both input styles must give the same answer, or one of them is lying."""
    bare = client.post("/api/orfs", json={"sequence": ORF_SEQ, "min_length": 30}).json()
    fasta = client.post("/api/orfs", json={"sequence": f">x\n{ORF_SEQ}", "min_length": 30}).json()
    assert bare["total"] == fasta["total"] == 1
    assert bare["orfs"][0]["length_nt"] == fasta["orfs"][0]["length_nt"]


def test_invalid_nucleotides_rejected_by_new_endpoints():
    r = client.post("/api/orfs", json={"sequence": "ATGCITGC", "min_length": 3})
    assert r.status_code == 422
