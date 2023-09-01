params.cores = 16

// Params Defaults in juno
params.refGenome = "/work/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta"
params.circos = "/opt/circos-0.69-2/bin/circos"
params.loci = "reference_data/copy_number/GermlineHetPon.37.vcf.gz"
params.gcProfile = "reference_data/copy_number/GC_profile.1000bp.37.cnp"
params.ensemblDataDir = "reference_data/common/ensembl_data"
params.genomeVersion = "V37"
params.diploidRegions = "copy_number/DiploidRegions.37.bed.gz"


log.info """\
    HMFTOOLS - PURPLE
    ========================================
    Params:
    ----------------------------------------
    tumor      : ${params.tumor}
    tumorBam   : ${params.tumorBam}
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
    tag "AMBER on ${params.tumor}"
    publishDir "${params.outdir}/amber", mode: 'copy'
    cpus params.cores
    memory '4 GB'
    time '1h'

    input:
    path tumor
    path tumorBam

    output:
    path "${tumor}.amber.baf.tsv.gz"
    path "${tumor}.amber.baf.pcf"
    path "${tumor}.amber.qc"

    script:
    """
    hmftools amber \
        -tumor ${tumor} \
        -tumor_bam ${tumorBam} \
        -output_dir \$PWD \
        -threads ${params.cores} \
        -loci ${params.loci} \
        -ref_genome_version ${params.genomeVersion}
    """
}

process runCobalt {
    tag "COBALT on ${params.tumor}"
    publishDir "${params.outdir}/cobalt", mode: 'copy'
    cpus params.cores
    memory '4 GB'
    time '1h'

    input:
    path tumor
    path tumorBam

    output:
    path "${tumor}.cobalt.ratio.tsv.gz"
    path "${tumor}.cobalt.ratio.pcf"

    script:
    """
    hmftools cobalt \
        -tumor ${tumor} \
        -tumor_bam ${tumorBam} \
        -output_dir \$PWD \
        -threads ${params.cores} \
        -gc_profile ${params.gcProfile} \
        -tumor_only_diploid_bed ${params.diploidRegions}
    """
}

process runPurple {
    tag "PURPLE on ${params.tumor}"
    publishDir "${params.outdir}/purple", mode: 'copy'
    cpus params.cores
    memory '4 GB'
    time '1h'

    input:
    path tumor
    path amber_baf_tsv
    path amber_baf_pcf
    path amber_qc
    path cobalt_ratio_tsv
    path cobalt_ratio_pcf

    output:
    path "${tumor}.purple.purity.tsv"
    path "${tumor}.purple.qc"
    path "${tumor}.purple.purity.range.tsv"
    path "${tumor}.purple.cnv.somatic.tsv"
    path "${tumor}.purple.cnv.gene.tsv"
    path "${tumor}.purple.segment.tsv"
    path "${tumor}.purple.somatic.clonality.tsv"
    path "plot/${tumor}.segment.png"
    path "plot/${tumor}.copynumber.png"
    path "plot/${tumor}.circos.png"
    path "plot/${tumor}.map.png"
    path "plot/${tumor}.input.png"
    path "plot/${tumor}.purity.range.png"

    script:
    """
    hmftools purple \
    -tumor ${tumor} \
    -amber ${params.outdir}/amber \
    -cobalt ${params.outdir}/cobalt \
    -output_dir \$PWD \
    -gc_profile ${params.gcProfile} \
    -ref_genome ${params.refGenome} \
    -ref_genome_version ${params.genomeVersion} \
    -ensembl_data_dir ${params.ensemblDataDir} \
    -circos ${params.circos}
    """
}

workflow {
    tumor = Channel.fromPath(params.tumor)
    tumorBam = Channel.fromPath(params.tumorBam)

    runAmber(tumor, tumorBam)
    runCobalt(tumor, tumorBam)
    runPurple(tumor, runAmber.out, runCobalt.out)
}
