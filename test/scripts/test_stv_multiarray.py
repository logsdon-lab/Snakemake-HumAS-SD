import pytest
from test.scripts.helpers import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            "test/data/stv_multiarray/all/input.bed",
            "test/data/stv_multiarray/all/expected.bed",
        ),
        (
            "test/data/stv_multiarray/S3CXH1L.4/input.bed",
            "test/data/stv_multiarray/S3CXH1L.4/expected.bed",
        ),
    ],
)
def test_stv_multiarray(infile: str, expected: str):
    run_integration_test(
        "python",
        "workflow/scripts/stv_multiarray.py",
        "-i",
        infile,
        expected_output=expected,
    )
