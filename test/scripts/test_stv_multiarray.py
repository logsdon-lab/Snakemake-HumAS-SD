import pytest
from test.scripts.helpers import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected", "chrom"],
    [
        (
            "test/data/stv_multiarray/all/input.bed",
            "test/data/stv_multiarray/all/expected.bed",
            "chr19",
        ),
        (
            "test/data/stv_multiarray/S3CXH1L.4/input.bed",
            "test/data/stv_multiarray/S3CXH1L.4/expected.bed",
            "chrX",
        ),
        (
            "test/data/stv_multiarray/chrX_hybrid_ort/input.bed",
            "test/data/stv_multiarray/chrX_hybrid_ort/expected.bed",
            "chrX",
        ),
    ],
)
def test_stv_multiarray(infile: str, expected: str, chrom: str):
    run_integration_test(
        "python",
        "workflow/scripts/stv_multiarray.py",
        "-i",
        infile,
        "-c",
        chrom,
        expected_output=expected,
    )
