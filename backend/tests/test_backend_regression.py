import pytest
from fastapi.testclient import TestClient
from app.main import app
from app.core.validation import wrap_raw_as_fasta, normalize_sequence
from app.core.alignment import perform_alignment
from app.core.codons import calculate_rscu

client = TestClient(app)


def test_raw_sequence_auto_wrap_in_analyze():
    response = client.post("/api/analyze", json={"sequence": "ATGCATGC"})
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert data["statistics"]["total_length"] == 8


def test_alignment_identical_and_gapped():
    # Identical ACGT vs ACGT should give 100% identity
    res1 = perform_alignment("ACGT", "ACGT", mode="global")
    assert res1["identity_percentage"] == 100.0

    # One gap case
    res2 = perform_alignment("ACGT", "AGT", mode="global")
    assert res2["identity_percentage"] == 75.0


def test_rscu_calculation_alanine():
    # GCTGGC -> GCT (Ala), GGC (Gly)
    seq = "GCTGGCGCTGGC"
    rscu = calculate_rscu(seq)
    assert "GCT" in rscu
    assert rscu["GCT"] == 4.0


def test_invalid_nucleotide_rejection():
    response = client.post("/api/analyze", json={"sequence": "ATGCITGC"})
    assert response.status_code == 422
