

checkpoint generate_monomers:
    input:
        script=workflow.source_path("../scripts/parse_hmm.py"),
        # TODO: Replace with fasta dir?
        hmm=HMM_PROFILE,
    output:
        directory(join(OUTPUT_DIR, "monomers")),
    benchmark:
        join(BMK_DIR, "generate_monomers.txt")
    log:
        join(LOG_DIR, "generate_monomers.log"),
    conda:
        "../envs/env.yaml"
    shell:
        """
        python {input.script} -i {input.hmm} -o {output} 2> {log}
        """


def get_monomer_by_chr(wc):
    output_dir = checkpoints.generate_monomers.get(**wc).output
    chr_name = re.search(RGX_CHR, str(wc.fname)).group()
    return os.path.join(str(output_dir), f"{chr_name}.fa")


rule run_stringdecomposer:
    input:
        monomers=get_monomer_by_chr,
        seq=os.path.join(INPUT_DIR, "{fname}.fa"),
    output:
        alt=join(OUTPUT_DIR, "{fname}", "final_decomposition_alt.tsv"),
        raw=join(OUTPUT_DIR, "{fname}", "final_decomposition_raw.tsv"),
        final=join(OUTPUT_DIR, "{fname}", "final_decomposition.tsv"),
        log=join(OUTPUT_DIR, "{fname}", "stringdecomposer.log"),
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
        stringdecomposer -t {threads} {input.seq} {input.monomers} -o {params.output_dir} &> {log}
        """


rule convert_to_bed9:
    input:
        rules.run_stringdecomposer.output.final,
    output:
        join(OUTPUT_DIR, "{fname}", "final_decomposition.bed"),
    log:
        join(LOG_DIR, "finalize_output_{fname}.log"),
    params:
        thr=IDENT_THR,
    conda:
        "../envs/env.yaml"
    shell:
        """
        awk -v OFS="\\t" -v QT="'" '{{
            strand=($2 ~ QT) ? "-" : "+";
            # Remove HOR enumerator and single quote.
            gsub(":.+|"QT, "", $2);
            # Filter suboptimal calls.
            if ($12 != "+" || $5 < {params.thr}) {{next}};
            print $1, $3, $4, $2, $5, strand, $3, $4, "0,0,0"
        }}' {input} > {output} 2> {log}
        """


rule stringdecomposer_all:
    input:
        rules.generate_monomers.output,
        expand(rules.run_stringdecomposer.output, zip, fname=FNAMES, chr=CHRS),
        expand(rules.convert_to_bed9.output, zip, fname=FNAMES, chr=CHRS),
    default_target: True
