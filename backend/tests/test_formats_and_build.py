"""GenBank parsing, and the health endpoint's deploy-identification fields.

The health assertions exist because a stale deployment once cost hours: the API
gave no way to tell which commit it was running, so a five-commit-behind deploy
was indistinguishable from a code defect.
"""
from fastapi.testclient import TestClient

from app.core.io_fasta import detect_format, parse_fasta_string
from app.main import app

client = TestClient(app)

SEQ = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

GENBANK = """LOCUS       TESTSEQ                 39 bp    DNA     linear   UNA 03-SEP-2026
DEFINITION  Synthetic test record.
ACCESSION   TESTSEQ
VERSION     TESTSEQ.1
FEATURES             Location/Qualifiers
     source          1..39
ORIGIN
        1 atggccattg taatgggccg ctgaaagggt gcccgatag
//
"""

FASTQ = "@r1\n" + SEQ + "\n+\n" + ("I" * len(SEQ)) + "\n"


def test_detect_format_by_content():
    assert detect_format(">x\nATGC") == "fasta"
    assert detect_format(FASTQ) == "fastq"
    assert detect_format(GENBANK) == "genbank"


def test_detect_format_falls_back_to_extension():
    assert detect_format("ATGC", "sample.gbk") == "genbank"
    assert detect_format("ATGC", "sample.fa") == "fasta"


def test_genbank_parses_into_records():
    records = parse_fasta_string(GENBANK)
    assert len(records) == 1
    assert str(records[0].seq).upper() == SEQ


def test_genbank_and_fasta_give_the_same_answer():
    """The same sequence in two formats must analyse identically."""
    gb = client.post("/api/analyze", json={"file_content": GENBANK}).json()
    fa = client.post("/api/analyze", json={"file_content": f">x\n{SEQ}"}).json()
    assert gb["statistics"]["total_length"] == fa["statistics"]["total_length"] == 39
    assert gb["statistics"]["gc_content_percentage"] == fa["statistics"]["gc_content_percentage"]
    assert gb["kmer_analysis"]["counts"] == fa["kmer_analysis"]["counts"]


def test_health_reports_the_deployed_commit():
    body = client.get("/api/health").json()
    for field in ("commit", "commit_short", "branch", "commit_source", "kmer_engine"):
        assert field in body, f"/api/health must expose {field}"
    # "unknown" is an acceptable honest answer; a missing field is not.
    assert body["commit"]
    assert body["kmer_engine"] in {"native", "python"}


def test_health_never_claims_a_kmer_engine_it_is_not_using():
    """The reported engine must match what actually ran, not a hoped-for value."""
    from app.core.kmer_native import backend_status
    assert client.get("/api/health").json()["kmer_engine"] == backend_status()["backend"]
    analyzed = client.post("/api/analyze", json={"sequence": SEQ}).json()
    assert analyzed["kmer_analysis"]["engine"] == backend_status()["backend"]
