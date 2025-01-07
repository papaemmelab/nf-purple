// See https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#arguments
process runPurple {
    tag "PURPLE on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/purple/purple_${params.minPurity}_${params.maxPurity}", mode: 'copy'
    cpus Math.min(params.cores as int, 4)
    memory params.memory
    time '1h'

    input:
    path amber_baf_tsv
    path amber_baf_pcf
    path amber_qc
    path cobalt_tumor_ratio_tsv
    path cobalt_tumor_ratio_pcf
    path cobalt_path
    path somatic_vcf

    output:
    path "${params.tumor}.purple.*", emit: purple_outfiles
    path "plot/${params.tumor}.*.png", emit: purple_plots

    script:
    def reference_args = params.normal ? """\\
        -reference ${params.normal}""" : ""
    def somatic_vcf_args = params.normal && somatic_vcf ? """\\
        -somatic_vcf ${somatic_vcf}""" : ""

    """
    purple \\
        -tumor ${params.tumor} ${reference_args} ${somatic_vcf_args} \\
        -amber ${params.outdir}/amber \\
        -cobalt ${cobalt_path} \\
        -output_dir \$PWD \\
        -gc_profile ${params.gcProfile} \\
        -ref_genome ${params.refGenome} \\
        -ref_genome_version ${params.genomeVersion} \\
        -ensembl_data_dir ${params.ensemblDataDir} \\
        -circos ${params.circos} \\
        -min_purity ${params.minPurity} \\
        -max_purity ${params.maxPurity}
    """.stripIndent()
}
