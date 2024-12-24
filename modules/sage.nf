// See https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md#usage
process runSage {
    tag "SAGE on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '4h'

    input:
    tuple val(tumor), path(tumorBam), path(tumorBai)
    tuple val(normal), path(normalBam), path(normalBai)

    output:
    path "sage/${tumor}_vs_${normal}.vcf.gz", emit: sage_vcf

    script:
    """
    mkdir -p sage

    sage \\
        -tumor ${params.tumor} \\
        -tumor_bam ${tumorBam} \\
        -reference ${params.normal} \\
        -reference_bam ${normalBam} \\
        -ref_genome ${params.refGenome} \\
        -ref_genome_version ${params.genomeVersion} \\
        -output_vcf sage/${tumor}_vs_${normal}.vcf.gz \\
        -threads ${params.cores} \\
        -ensembl_data_dir ${params.ensemblDataDir}
    """.stripIndent()
}
