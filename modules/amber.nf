// See https://github.com/hartwigmedical/hmftools/tree/master/amber
process runAmber {
    tag "AMBER on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    tuple val(tumor), path(tumorBam), path(tumorBai)
    tuple val(normal), path(normalBam), path(normalBai)

    output:
    path "amber/${tumor}.amber.baf.tsv.gz", emit: amber_baf_tsv
    path "amber/${tumor}.amber.baf.pcf", emit: amber_baf_pcf
    path "amber/${tumor}.amber.qc", emit: amber_qc
    path "amber/${tumor}.amber.contamination.vcf.gz", emit: amber_contamination_vcf, optional: true

    script:
    def reference_args = normal ? """\\
        -reference ${normal} \\
        -reference_bam ${normalBam} """  : ""

    """
    amber \\
        -tumor ${tumor} \\
        -tumor_bam ${tumorBam} ${reference_args} \\
        -output_dir amber \\
        -threads ${params.cores} \\
        -loci ${params.loci} \\
        -ref_genome_version V${params.genomeVersion}
    """.stripIndent()
}
