

rule generate_monomers:
    input:
        script=workflow.source_path("../scripts/parse_hmm.py"),
        hmm=HMM_PROFILE,
    output:
        join(OUTPUT_DIR, "monomers", "{chrom}.fa"),
    benchmark:
        join(BMK_DIR, "generate_monomers_{chrom}.txt")
    params:
        # If multi chromosome, split and add to library. ex. chr1
        chrs=lambda wc: " ".join(wc.chrom.split("-")),
    log:
        join(LOG_DIR, "generate_monomers_{chrom}.log"),
    conda:
        "../envs/env.yaml"
    shell:
        """
        python {input.script} -i {input.hmm} -o {output} -c {params.chrs} 2> {log}
        """


def get_monomers(wc):
    if config.get("hmm_profile"):
        chr_name = get_chrom(wc.fname)
        monomer_fa = expand(rules.generate_monomers.output, chrom=chr_name)
    elif config.get("monomer_dir"):
        # Provide custom library matching sequence name.
        monomer_fa = os.path.join(config["monomer_dir"], f"{wc.fname}.fa")
    else:
        raise ValueError("No monomer source provided.")

    return monomer_fa


rule run_stringdecomposer:
    input:
        monomers=get_monomers,
        seq=os.path.join(INPUT_DIR, "{fname}.fa"),
    output:
        alt=join(OUTPUT_DIR, "{fname}", "final_decomposition_alt.tsv"),
        raw=join(OUTPUT_DIR, "{fname}", "final_decomposition_raw.tsv"),
        final=join(OUTPUT_DIR, "{fname}", "final_decomposition.tsv"),
        log=join(OUTPUT_DIR, "{fname}", "stringdecomposer.log"),
    resources:
        mem=config.get("mem", "8GB"),
    params:
        output_dir=lambda wc, output: dirname(str(output.alt)),
    benchmark:
        join(BMK_DIR, "run_stringdecomposer_{fname}.txt")
    log:
        join(LOG_DIR, "run_stringdecomposer_{fname}.log"),
    threads: config["threads"]
    conda:
        "../envs/env.yaml"
    shell:
        """
        ( stringdecomposer -t {threads} {input.seq} {input.monomers} -o {params.output_dir} || true ) &> {log}
        # Allow to continue on failures.
        touch {output}
        """


rule convert_to_bed9:
    input:
        rules.run_stringdecomposer.output.final,
    output:
        join(OUTPUT_DIR, "{fname}", "final_decomposition.bed"),
    log:
        join(LOG_DIR, "finalize_output_{fname}.log"),
    conda:
        "../envs/env.yaml"
    shell:
        """
        awk -v OFS="\\t" -v QT="'" '{{
            strand=($2 ~ QT) ? "-" : "+";
            # Remove HOR enumerator and single quote.
            gsub(":.+|"QT, "", $2);
            # Filter suboptimal calls.
            if ($12 != "+") {{next}};
            print $1, $3, $4, $2, $5, strand, $3, $4, "0,0,0"
        }}' {input} > {output} 2> {log}
        """


rule stringdecomposer_all:
    input:
        expand(rules.run_stringdecomposer.output, zip, fname=FNAMES),
        expand(rules.convert_to_bed9.output, zip, fname=FNAMES),
    default_target: True
