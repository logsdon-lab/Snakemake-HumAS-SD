import re
import gzip
import argparse
import itertools
from collections import Counter, defaultdict

RGX_NAME = re.compile(r"NAME\s\s(.*?)\n")
RGX_NT = re.compile(r"(\d)\s([atgc])")
RGX_CHRS = re.compile(r"C([\d\/XYM]+)H")
CHRS = [f"{i}" for i in itertools.chain(range(1, 23), ("X", "Y", "M"))]
# https://www.nature.com/articles/s41586-023-05976-y
ACRO_LD_CHRS = {"13", "14", "21", "22"}

ap = argparse.ArgumentParser(description="Convert HMM profile to fasta.")
ap.add_argument(
    "-i",
    "--input_hmm",
    help="Input HMM model file.",
    required=True,
    type=str,
)
ap.add_argument("-c", "--chromosomes", nargs="*", help="Chromosomes", default=RGX_CHRS)
ap.add_argument("-o", "--outfile", help="Output fasta.", type=str, required=True)

args = ap.parse_args()

if args.input_hmm.endswith(".gz"):
    with gzip.open(args.input_hmm) as fh:
        hmm = fh.read().decode().split("//")
else:
    with open(args.input_hmm) as fh:
        hmm = fh.read().split("//")

outfile = args.outfile
chromosomes = set(args.chromosomes)

# http://eddylab.org/software/hmmer/Userguide.pdf

hors_done = Counter()
all_chrom_mons = defaultdict(set)
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

    hor_seq = "".join(nt for _, nt in re.findall(RGX_NT, hmm_rec)).upper()
    # Is ancestral monomer. Add to all chrs.
    if not hor_name.startswith("S"):
        for chr_name in chromosomes:
            all_chrom_mons[chr_name].add((new_hor_name, hor_seq))
        continue

    mtch_chr_name = re.search(RGX_CHRS, hor_name)
    chr_names = mtch_chr_name.group(1).split("/") if mtch_chr_name else [hor_name]

    # Add multiple (ex. 1/15/19) to each chr file.
    for chr_name in chr_names:
        all_chrom_mons[chr_name].add((new_hor_name, hor_seq))
        # Add all acro mons with low LD in PHR regions that can recombine with one another.
        if chr_name in ACRO_LD_CHRS:
            for other_acro_chr in ACRO_LD_CHRS.difference(set([chr_name])):
                all_chrom_mons[other_acro_chr].add((new_hor_name, hor_seq))

    hors_done[hor_name] += 1

with open(outfile, "wt") as fh:
    for chr_name, hor_mons in all_chrom_mons.items():
        if chr_name not in chromosomes:
            continue
        for hor_mon, hor_mon_seq in hor_mons:
            fh.write(f">{hor_mon}\n")
            fh.write(f"{hor_mon_seq}\n")
