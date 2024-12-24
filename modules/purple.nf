// See https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#arguments
process runPurple {
    tag "PURPLE on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/purple/purple_${params.minPurity}_${params.maxPurity}", mode: 'copy'
    cpus params.cores
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
    path "${params.tumor}.purple.cnv.gene.tsv"
    path "${params.tumor}.purple.cnv.somatic.tsv"
    path "${params.tumor}.purple.germline.deletion.tsv", optional: true
    path "${params.tumor}.purple.purity.range.tsv"
    path "${params.tumor}.purple.purity.tsv"
    path "${params.tumor}.purple.qc"
    path "${params.tumor}.purple.segment.tsv"
    path "${params.tumor}.purple.somatic.clonality.tsv"
    path "${params.tumor}.purple.somatic.hist.tsv", optional: true
    path "${params.tumor}.purple.somatic.vcf.gz", optional: true
    path "${params.tumor}.purple.somatic.vcf.gz.tbi", optional: true
    path "plot/${params.tumor}.circos.png"
    path "plot/${params.tumor}.copynumber.png"
    path "plot/${params.tumor}.input.png"
    path "plot/${params.tumor}.map.png"
    path "plot/${params.tumor}.purity.range.png"
    path "plot/${params.tumor}.segment.png"
    path "plot/${params.tumor}.somatic.clonality.png", optional: true
    path "plot/${params.tumor}.somatic.png", optional: true
    path "plot/${params.tumor}.somatic.rainfall.png", optional: true

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
