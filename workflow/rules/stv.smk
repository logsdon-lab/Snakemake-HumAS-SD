

rule generate_stv:
    input:
        script=workflow.source_path("../scripts/stv_multiarray.py"),
        hor_bed=join(OUTPUT_DIR, "{fname}", "final_decomposition.bed"),
    output:
        stv_row_bed=join(OUTPUT_DIR, "{fname}", "stv_row.bed"),
    conda:
        "../envs/env.yaml"
    benchmark:
        join(BMK_DIR, "generate_stv_{fname}.txt")
    log:
        join(LOG_DIR, "generate_stv_{fname}.log"),
    shell:
        """
        python {input.script} -i {input.hor_bed} -o {output.stv_row_bed} 2> {log}
        """


rule stv_all:
    input:
        expand(rules.generate_stv.output, zip, fname=FNAMES, chr=CHRS),
    default_target: True
