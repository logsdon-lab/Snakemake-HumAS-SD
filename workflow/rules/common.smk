import re
import sys
from os.path import join, dirname


def get_chrom(fname: str) -> str:
    return "-".join(re.findall(RGX_CHR, fname))


def extract_fa_fnames_and_chr(input_dir: str) -> tuple[list[str], list[str]]:
    fnames = glob_wildcards(join(input_dir, "{fname}.fa")).fname
    filtered_fnames, chrs = [], []
    for fname in fnames:
        chr_name = "-".join(re.findall(RGX_CHR, fname))
        if not chr_name:
            continue
        filtered_fnames.append(fname)
        chrs.append(chr_name)

    assert len(filtered_fnames) == len(
        chrs
    ), f"One or more fa files in {input_dir} does not contain a chromosome in its name."

    return filtered_fnames, chrs
