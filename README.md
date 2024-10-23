# Snakemake-HumAS-SD
Workflow to use [`stringdecomposer`](https://github.com/ablab/stringdecomposer) to annotate alpha-satellite HOR monomers.

HMM profile:
* https://github.com/fedorrik/HumAS-HMMER_for_AnVIL/blob/main/AS-HORs-hmmer3.4-071024.hmm

Based on annotations from:
* Uralsky, L. I., Shepelev, V. A., Alexandrov, A. A., Yurov, Y. B., Rogaev, E. I., & Alexandrov, I. A. (2019). Classification and monomer-by-monomer annotation dataset of suprachromosomal family 1 alpha satellite higher-order repeats in hg38 human genome assembly. Data in brief, 24, 103708. https://doi.org/10.1016/j.dib.2019.103708

### Usage
```bash
conda env create -f env.yaml --name humas_sd
conda activate humas_sd
snakemake -np --sdm conda --configfile config/config.yaml -c 1
```

### Config
```yaml
# Input directory with fa files.
input_dir: "cens"
# Various directories.
output_dir: "results"
benchmarks_dir: "benchmarks"
logs_dir: "logs"
# Input HMM profile with all possible monomers.
# May remove.
hmm_profile: "data/AS-HORs-hmmer.hmm"
# Threads passed to stringdecomposer
threads: 4
```

### Test
Two centromere HOR annotation bed files.
1. Output from HumAS-HMMER_for_AnVIL using https://github.com/logsdon-lab/CenMAP/blob/main/data/models/AS-HORs-hmmer3.0-170921.hmm.
    * See also https://github.com/logsdon-lab/Snakemake-HumAS-HMMER
2. Output from stringdecomposer.

> [!NOTE]
> From HGSVC3 sample , `HG01596_chr3_haplotype1-0000018:7760823-18212765`, which was extracted via CenMAP v0.1.5.

```bash
python test/scripts/compare.py
```
```
Number of shared HORs annotated: 22894
Average Length Difference: 1.1741504324277103
Live:
        Correct: 18245
        Total: 18327
        Perc: 99.55257270693512
Other:
        Correct: 3880
        Total: 4567
        Perc: 84.9573023866871
```
