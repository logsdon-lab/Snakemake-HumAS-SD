import os
import intervaltree as it
import polars as pl


def main():
    wd = os.path.dirname(__file__)
    df_stringdecomposer = pl.read_csv(
        os.path.join(
            wd, "../data/HG01596_chr3_haplotype1-0000018:7760823-18212765_stringdecomposer.bed"
        ),
        separator="\t",
        new_columns=[
            "read-name",
            "best-monomer",
            "start-pos",
            "end-pos",
            "identity",
            "second-best-monomer",
            "second-best-monomer-identity",
            "homo-best-monomer",
            "homo-identity",
            "homo-second-best-monomer",
            "homo-second-best-monomer-identity",
            "reliability",
        ],
        has_header=False
    ).with_columns(
        # Remove other monomers count
        pl.col("best-monomer").str.replace(r":\d+", "")
    )
    df_humas_hmmer = pl.read_csv(
        os.path.join(
            wd, "../data/HG01596_chr3_haplotype1-0000018:7760823-18212765_humas_hmmer.bed"
        ),
        separator="\t",
        new_columns=[
            "chrom",
            "st",
            "end",
            "name",
            "score",
            "strand",
            "thick_st",
            "thick_end",
            "item_rgb"
        ],
        has_header=False
    ).with_columns(
        # Make similar format to stringdecomposer
        # Strand (-) -> Monomer'
        new_name=pl.col("name").cast(pl.String) + pl.when(pl.col("strand") == "-").then(pl.lit("'")).otherwise(pl.lit(""))
    )
    intervals_humas_hmmer = it.IntervalTree(
        it.Interval(*r) for r in df_humas_hmmer.select("st", "end", "new_name").rows()
    )
    intervals_stringdecomposer = it.IntervalTree(
        it.Interval(*r) for r in df_stringdecomposer.select("start-pos", "end-pos", "best-monomer").rows()
    ) 

    correct_l = 0
    total_l = 0
    correct_o = 0
    total_o = 0
    sum_len_diff = 0
    for interval in sorted(intervals_humas_hmmer.iter()):
        overlap_sd = intervals_stringdecomposer.overlap(interval)
        if not overlap_sd:
            continue
        overlap_interval = next(iter(overlap_sd))
        # Keep track of length difference.
        sum_len_diff += abs(interval.length() - overlap_interval.length())

        assignment_identical = overlap_interval.data == interval.data 
        if "L" in interval.data:
            correct_l += assignment_identical
            total_l += 1
        else:
            correct_o += assignment_identical
            total_o += 1

    num_hor = total_o + total_l
    avg_len_diff = sum_len_diff / num_hor

    print(f"Number of shared HORs annotated: {num_hor}")
    print(f"Average Length Difference: {avg_len_diff}")
    print("Live:")
    print(f"\tCorrect: {correct_l}\n\tTotal: {total_l}\n\tPerc: {(correct_l / total_l) * 100}")
    print("Other:")
    print(f"\tCorrect: {correct_o}\n\tTotal: {total_o}\n\tPerc: {(correct_o / total_o) * 100}")


if __name__ == "__main__": 
    raise SystemExit(main())