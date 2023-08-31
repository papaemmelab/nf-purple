params.cores = 16

// Params Defaults in juno
params.refGenome = "/work/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta"
params.circos = "/opt/circos-0.69-2/bin/circos"
params.loci = "reference_data/copy_number/GermlineHetPon.37.vcf.gz"
params.gcProfile = "reference_data/copy_number/GC_profile.1000bp.37.cnp"
params.ensemblDataDir = "reference_data/common/ensembl_data"
params.genomeVersion = "V37"


log.info """\
    A L L S O R T S    C L A S S I F I E R S
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
    dir 'amber'

    script:
    """
    hmftools amber \
        -tumor ${tumor} \
        -tumor_bam ${tumorBam} \
        -output_dir ${output} \
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
    path 'outputs/cobalt'

    script:
    """
    hmftools cobalt \
        -tumor ${tumor} \
        -tumor_bam ${tumorBam} \
        -output_dir ${output} \
        -threads ${params.cores} \
        -gc_profile ${params.gcProfile}
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
    path amberOutput
    path cobaltOutput

    output:
    path 'outputs/purple'

    script:
    """
    hmftools purple \
    -tumor ${tumor} \
    -amber ${amberOutput} \
    -cobalt ${cobaltOutput} \
    -output_dir ${outidr} \
    -gc_profile ${params.gcProfile} \
    -ref_genome ${params.refGenome} \
    -ref_genome_version ${params.genomeVersion} \
    -ensembl_data_dir ${params.ensemblDataDir} \
    -circos ${params.circos}
    """
}


workflow {
    runAmber(tumor, tumorBam, loci, genomeVersion)
    runCobalt(tumor, tumorBam, gcProfile)
    runPurple(tumor, runAmber.out, runCobalt.out, gcProfile, refGenome, ensemblDataDir, circos)
}