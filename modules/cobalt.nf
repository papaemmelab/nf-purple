// See https://github.com/hartwigmedical/hmftools/tree/master/cobalt#mandatory-arguments
process runCobalt {
    tag "COBALT on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    tuple val(tumor), path(tumorBam), path(tumorBai)
    tuple val(normal), path(normalBam), path(normalBai)

    output:
    path "cobalt/${tumor}.cobalt.ratio.tsv.gz", emit: cobalt_tumor_ratio_tsv
    path "cobalt/${tumor}.cobalt.ratio.pcf", emit: cobalt_tumor_ratio_pcf
    path "cobalt/${normal}.cobalt.ratio.pcf", emit: cobalt_normal_ratio_pcf, optional: true

    script:
    def reference_args = normal ? """\\
            -reference ${normal} \\
            -reference_bam ${normalBam}""" : """\\
            -tumor_only_diploid_bed ${params.diploidRegions}"""

    """
    cobalt \\
        -tumor ${tumor} \\
        -tumor_bam ${tumorBam} ${reference_args} \\
        -output_dir cobalt \\
        -threads ${params.cores} \\
        -gc_profile ${params.gcProfile}
    """.stripIndent()
}

process binCobalt {
    // Only when Unmatched
    tag "COBALT BIN on ${params.tumor}" + (params.normal ? " and ${params.normal}" : "")
    publishDir "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR", mode: 'copy'
    cpus 1
    memory '4G'
    time '1h'

    input:
    path cobalt_tumor_ratio_tsv
    path cobalt_tumor_ratio_pcf

    output:
    path "${params.tumor}.cobalt.ratio.tsv.gz", emit: cobalt_tumor_ratio_tsv
    path "${params.tumor}.cobalt.ratio.pcf", emit: cobalt_tumor_ratio_pcf

    script:
    """
    bin_cobalt.py \\
        --in_pcf ${cobalt_tumor_ratio_pcf} \\
        --bin_probes ${params.binProbes} \\
        --bin_log_r ${params.binLogR}
    """.stripIndent()
}
