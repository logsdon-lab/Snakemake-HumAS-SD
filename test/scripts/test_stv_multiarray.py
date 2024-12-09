import pytest
from test.scripts.helpers import run_integration_test


@pytest.mark.parametrize(
    ["infile", "expected"],
    [
        (
            "test/data/stv_multiarray/input.bed",
            "test/data/stv_multiarray/expected.bed",
        )
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
