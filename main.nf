params.cores = 4

// Params Defaults in juno
params.refGenome = "/work/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta"
params.genomeVersion = "V37"
params.circos = "/opt/circos-0.69-2/bin/circos"
params.loci = "/data/copy_number/GermlineHetPon.37.vcf.gz"
params.gcProfile = "/data/copy_number/GC_profile.1000bp.37.cnp"
params.ensemblDataDir = "/data/common/ensembl_data"
params.diploidRegions = "/data/copy_number/DiploidRegions.37.bed.gz"
params.normal = null
params.normalBam = null


log.info """\
    HMFTOOLS - PURPLE
    ========================================
    Params:
    ----------------------------------------
    tumor      : ${params.tumor}
    tumorBam   : ${params.tumorBam}
    normal      : ${params.normal}
    normalBam   : ${params.normalBam}
    outdir     : ${params.outdir}
    cores      : ${params.cores}
    ========================================
    Workflow:
    ----------------------------------------
    Project    : ${workflow.projectDir}
    Cmd line   : ${workflow.commandLine}
    """
    .stripIndent()


process runAmber {
    tag "AMBER on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/amber", mode: 'copy'
    cpus params.cores
    memory '4 GB'
    time '1h'

    input:
    val tumor
    val normal
    path tumorBam
    path normalBam

    output:
    path "${tumor}.amber.baf.tsv.gz", emit: amber_baf_tsv
    path "${tumor}.amber.baf.pcf", emit: amber_baf_pcf
    path "${tumor}.amber.qc", emit: amber_qc
    path "${tumor}.amber.contamination.vcf.gz", emit: amber_contamination_vcf
    path "${normal}.amber.snp.vcf.gz", emit: amber_normal_snp_vcf, optional: true
    path "${normal}.amber.homozygousregion.tsv", emit: amber_normal_homozygousregion_tsv, optional: true

    script:
    def reference_args = normal ? "-reference ${normal} \\\n    -reference_bam ${normalBam}"  : ""
    """
    amber \
        -tumor ${tumor} \
        -tumor_bam ${tumorBam} \
        ${reference_args} \
        -output_dir \$PWD \
        -threads ${params.cores} \
        -loci ${params.loci} \
        -ref_genome_version ${params.genomeVersion}
    """.stripIndent()
}

process runCobalt {
    tag "COBALT on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/cobalt", mode: 'copy'
    cpus params.cores
    memory '4 GB'
    time '1h'

    input:
    val tumor
    val normal
    path tumorBam
    path normalBam

    output:
    path "${tumor}.cobalt.ratio.tsv.gz", emit: cobalt_ratio_tsv
    path "${tumor}.cobalt.ratio.pcf", emit: cobalt_ratio_pcf
    path "${normal}.cobalt.ratio.pcf", emit: cobalt_normal_ratio_pcf, optional: true

    script:
    def reference_args = normal ? "-reference ${normal} \\\n    -reference_bam ${normalBam}"  : "-tumor_only_diploid_bed ${params.diploidRegions}"
    """
    cobalt \
        -tumor ${tumor} \
        -tumor_bam ${tumorBam} \
        ${reference_args} \
        -output_dir \$PWD \
        -threads ${params.cores} \
        -gc_profile ${params.gcProfile}
    """.stripIndent()
}

process runPurple {
    tag "PURPLE on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/purple", mode: 'copy'
    cpus params.cores
    memory '4 GB'
    time '1h'

    input:
    val tumor
    val normal
    path amber_baf_tsv
    path amber_baf_pcf
    path amber_qc
    path amber_contamination_vcf
    path amber_normal_snp_vcf
    path amber_normal_homozygousregion_tsv
    path cobalt_ratio_tsv
    path cobalt_ratio_pcf
    path cobalt_normal_ratio_pcf

    output:
    path "${tumor}.purple.purity.tsv", emit: purple_purity_tsv
    path "${tumor}.purple.qc", emit: purple_qc
    path "${tumor}.purple.purity.range.tsv", emit: purple_purity_range_tsv
    path "${tumor}.purple.cnv.somatic.tsv", emit: purple_cnv_somatic_tsv
    path "${tumor}.purple.cnv.gene.tsv", emit: purple_cnv_gene_tsv
    path "${tumor}.purple.segment.tsv", emit: purple_segment_tsv
    path "${tumor}.purple.somatic.clonality.tsv", emit: purple_somatic_clonality_tsv
    path "plot/${tumor}.segment.png", emit: purple_segment_png
    path "plot/${tumor}.copynumber.png", emit: purple_copynumber_png
    path "plot/${tumor}.circos.png", emit: purple_circos_png
    path "plot/${tumor}.map.png", emit: purple_map_png
    path "plot/${tumor}.input.png", emit: purple_input_png
    path "plot/${tumor}.purity.range.png", emit: purple_purity_range_png

    script:
    def reference_args = normal ? "-reference ${normal}"  : ""
    """
    purple \
        -tumor ${tumor} \
        ${reference_args} \
        -amber ${params.outdir}/amber \
        -cobalt ${params.outdir}/cobalt \
        -output_dir \$PWD \
        -gc_profile ${params.gcProfile} \
        -ref_genome ${params.refGenome} \
        -ref_genome_version ${params.genomeVersion} \
        -ensembl_data_dir ${params.ensemblDataDir} \
        -circos ${params.circos}
    """.stripIndent()
}

workflow {
    tumor = Channel.value(params.tumor)
    normal = Channel.value(params.normal)
    tumorBam = Channel.fromPath(params.tumorBam)
    normalBam = Channel.fromPath(params.normalBam)

    runAmber(tumor, normal, tumorBam, normalBam)
    runCobalt(tumor, normal, tumorBam, normalBam)
    runPurple(tumor, normal, runAmber.out, runCobalt.out)
}


workflow.onComplete {
    log.info (
      workflow.success
        ? "\nDone! Purple ran successfully. See the results in: ${params.outdir}\n"
        : "\nOops .. something went wrong\n"
    )
}