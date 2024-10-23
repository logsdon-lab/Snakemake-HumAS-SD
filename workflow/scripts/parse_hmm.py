import re
import os
import argparse
from collections import Counter
from typing import TextIO

RGX_NAME = re.compile(r"NAME\s\s(.*?)\n")
RGX_NT = re.compile(r"(\d)\s([atgc])")
RGX_CHRS = re.compile(r"C([\d\/XYM]+)H")

ap = argparse.ArgumentParser(description="Convert HMM model to fasta.")
ap.add_argument("-i", "--input_hmm", help="Input HMM model file.", required=True, type=argparse.FileType("r"))
ap.add_argument("-o", "--outdir", help="Output directory.", type=str, required=True)

args = ap.parse_args()

fh: TextIO = args.input_hmm
outdir = args.outdir
os.makedirs(outdir, exist_ok=True)

# http://eddylab.org/software/hmmer/Userguide.pdf
hmm = fh.read().split("//")

hors_done = Counter()
chr_fhs = {}
for hmm_rec in hmm:
    hor_name = re.search(RGX_NAME, hmm_rec)
    if not hor_name:
        continue
    hor_name = hor_name.group(1)

    # Add suffix number to create unique ids for HORs with multiple monomer sequences.
    if hor_name in hors_done:
        hor_cnt = hors_done[hor_name]
        new_hor_name = f"{hor_name}:{hor_cnt}"
    else:
        new_hor_name = hor_name

    mtch_chr_name = re.search(RGX_CHRS, hor_name)
    chr_names = mtch_chr_name.group(1).split("/") if mtch_chr_name else [hor_name]
    fa_fname_prefix = "chr" if mtch_chr_name else ""
    # Add multiple (ex. 1/15/19) to each chr file.
    for chr_name in chr_names:
        if not chr_fhs.get(chr_name):
            chr_fhs[chr_name] = open(
                os.path.join(outdir, f"{fa_fname_prefix}{chr_name}.fa"), "wt"
            )

        chr_fh = chr_fhs[chr_name]

        # TODO: Just take the first one.
        chr_fh.write(f">{new_hor_name}\n")
        hor_nts = ''.join(nt for _, nt in re.findall(RGX_NT, hmm_rec)).upper()
        chr_fh.write(f"{hor_nts}\n")

    hors_done[hor_name] += 1

for fh in chr_fhs.values():
    fh.close()
