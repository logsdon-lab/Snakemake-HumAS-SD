import re
import os
import argparse
import itertools
from collections import Counter
from typing import TextIO

RGX_NAME = re.compile(r"NAME\s\s(.*?)\n")
RGX_NT = re.compile(r"(\d)\s([atgc])")
RGX_CHRS = re.compile(r"C([\d\/XYM]+)H")
CHRS = [f"{i}" for i in itertools.chain(range(1, 23), ("X", "Y", "M"))]
ACRO_CHRS = {"13", "14", "15", "21", "22"}

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
all_chrom_fhs = {
    chrom: open(
        os.path.join(outdir, f"chr{chrom}.fa"), "wt"
    )
    for chrom in CHRS
}
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
    # Add multiple (ex. 1/15/19) to each chr file.
    for chr_name in chr_names:
        chrom_fh = all_chrom_fhs.get(chr_name)
        is_not_anc_mon = hor_name.startswith("S")
        hor_nts = ''.join(nt for _, nt in re.findall(RGX_NT, hmm_rec)).upper()

        chrom_fh = all_chrom_fhs.get(chr_name)
        # If not a chrom, is a ancestral mon.
        chrom_fhs = [chrom_fh] if chrom_fh else all_chrom_fhs.values()

        # Add all acro mons to each other.
        if chr_name in ACRO_CHRS:
            other_acro_chrs = ACRO_CHRS.difference(set([chr_name]))
            chrom_fhs.extend(all_chrom_fhs[acro_chr] for acro_chr in other_acro_chrs)

        for chrom_fh in chrom_fhs:
            chrom_fh.write(f">{new_hor_name}\n")
            chrom_fh.write(f"{hor_nts}\n")

    hors_done[hor_name] += 1

for chr_name, fh in all_chrom_fhs.items():
    fh.close()
