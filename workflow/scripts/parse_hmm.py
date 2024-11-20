import re
import os
import argparse
from collections import Counter, defaultdict
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
all_chrom_fhs = {}
all_anc_mons = defaultdict(list)
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
        chrom_fh = all_chrom_fhs.get(chr_name)
        is_not_anc_mon = hor_name.startswith("S")
        hor_nts = ''.join(nt for _, nt in re.findall(RGX_NT, hmm_rec)).upper()

        if not chrom_fh and is_not_anc_mon:
            all_chrom_fhs[chr_name] = open(
                os.path.join(outdir, f"{fa_fname_prefix}{chr_name}.fa"), "wt"
            )
            chrom_fh = all_chrom_fhs[chr_name]
            
        elif not is_not_anc_mon:
            all_anc_mons[chr_name].append(hor_nts)

        if is_not_anc_mon:
            chrom_fh.write(f">{new_hor_name}\n")
            chrom_fh.write(f"{hor_nts}\n")

    hors_done[hor_name] += 1

# Get ancestral monomers
for chr_name, chr_fh in all_chrom_fhs.items():
    # Write all ancestral monomers to all chr_fhs.
    for anc_mon_name, anc_mons in all_anc_mons.items():
        for anc_mon in anc_mons:
            chr_fh.write(f">{anc_mon_name}\n")
            chr_fh.write(f"{anc_mon}\n")

for chr_name, fh in all_chrom_fhs.items():
    fh.close()
