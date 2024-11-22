import re
import sys
from os.path import join, dirname


def extract_fa_fnames_and_chr(
    input_dir: str, *, filter_chr: str | None = None
) -> tuple[list[str], list[str]]:
    fnames = glob_wildcards(join(input_dir, "{fname}.fa")).fname
    filtered_fnames, chrs = [], []
    for fname in fnames:
        if mtch_chr_name := re.search(RGX_CHR, fname):
            chr_name = mtch_chr_name.group().strip("_")

            if not filter_chr:
                filtered_fnames.append(fname)
                chrs.append(chr_name)
                continue

            # Filter by chr.
            if chr_name != filter_chr:
                continue

            filtered_fnames.append(fname)
            chrs.append(chr_name)
        else:
            print("Cannot find chr name in {fname}. Skipping.", file=sys.stderr)

    assert len(filtered_fnames) == len(
        chrs
    ), f"One or more fa files in {input_dir} does not contain a chromosome in its name."

    return filtered_fnames, chrs
