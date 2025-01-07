include { runAmber  }   from './modules/amber'
include { runCobalt }   from './modules/cobalt'
include { binCobalt }   from './modules/cobalt'
include { runSage   }   from './modules/sage'
include { runPurple }   from './modules/purple'

def logMessage = """\
    HMFTOOLS - PURPLE
    ========================================
    Running Mode   : ${params.normal ? 'Matched' : 'Unmatched'}
    ----------------------------------------
    Params:
    ----------------------------------------
    tumor          : ${params.tumor}
    tumorBam       : ${params.tumorBam}
"""
logMessage += (params.normal && params.normalBam) ? """\
    normal         : ${params.normal}
    normalBam      : ${params.normalBam}
""" : ""
logMessage += """\
    outdir         : ${params.outdir}
    cores          : ${params.cores}
    memory         : ${params.memory}
    binProbes      : ${params.binProbes}
    binLogR        : ${params.binLogR}
    minPurity      : ${params.minPurity}
    maxPurity      : ${params.maxPurity}
    ensemblDataDir : ${params.ensemblDataDir}
    ========================================
    Workflow:
    ----------------------------------------
    Project        : ${workflow.projectDir}
    Cmd line       : ${workflow.commandLine}
"""

log.info(logMessage.stripIndent())

def printHelpMessage() {
    helpMessage = """\
    Usage:  nextflow run papaemmelab/nf-purple -r main [options]

    PURPLE is a purity/ploidy estimator designed for whole genome sequenced (WGS) data.
    It runs: Amber, Cobalt, Sage (matched mode only), and Purple.
    For more: https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md


    * Options:
        --tumor             Name of tumor sample (required)
        --tumorBam          Path to indexed tumor bam/cram file (required).
        --outdir            Path to the output directory. This directory will be
                            created if it does not already exist.
        --refGenome         Path to the reference genome fasta file (required
                            only when using CRAM files).
        --loci              Path to BAF loci vcf file (required).
        --gcProfile         Location of HMFtools GC Profile (required).
        --ensemblDataDir    Ensembl data file directory (required).
        --circos            Location of circos binary (required)
        --genomeVersion     Ref genome version, 37 or 38 (default: 37).
        --cores             Cores/CPUs to use by executor (default: 1).
        --memory            Memory to allocate by executor (default: 4 GB).
        --minPurity         Minimum purity (default: 0.08).
        --maxPurity         Maximum purity (default: 1.0).

    * Required for Matched Mode:
        --normal            Name of normal sample (required on matched mode)
        --normalBam         Path to normal bam/cram file (required on matched mode)

    * Required for Unmatched Mode:
        --diploidRegions    Diploid regions (required for unmatched mode)
        --binProbes         Max probe bin size, for unmatched mode (default: 0).
        --binLogR           Max probe logR diff to bin, for unmatched mode (default: 0).
    """
    log.info(helpMessage.stripIndent())
}

def NO_FILE = "${projectDir}/assets/NO_FILE"

workflow {
    if ( 
        params.help == true 
        || params.tumor == false 
        || params.tumorBam == false 
        || params.refGenome == false 
    ){
        printHelpMessage()
        exit 1
    }
    // Input Bams
    tumorTuple = tuple(params.tumor, params.tumorBam, params.tumorBam + ".bai")
    normalTuple = tuple(params.normal, params.normalBam, params.normalBam + ".bai")

    // Run Amber and Cobalt
    amberOutput = runAmber(tumorTuple, normalTuple)
    cobaltOutput = runCobalt(tumorTuple, normalTuple)

    // Bin Cobalt, if unmatched and if any bin param is provided
    cobaltOutdir = "${params.outdir}/cobalt"
    if ((!params.normal) && (params.binProbes != 0 || params.binLogR != 0)) {
        cobaltOutput = binCobalt(cobaltOutput.cobalt_tumor_ratio_tsv, cobaltOutput.cobalt_tumor_ratio_pcf)
        cobaltOutdir = "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR"
    }

    // Run Sage, if matched
    if (params.normal) {
        sageOutput = runSage(tumorTuple, normalTuple)
        somatic_vcf = sageOutput.sage_vcf
    } else {
        somatic_vcf = file(NO_FILE)
    }

    // Run Purple
    runPurple(
        amberOutput.amber_baf_tsv,
        amberOutput.amber_baf_pcf,
        amberOutput.amber_qc,
        cobaltOutput.cobalt_tumor_ratio_tsv,
        cobaltOutput.cobalt_tumor_ratio_pcf,
        cobaltOutdir,
        somatic_vcf,
    )
}

workflow.onComplete {
    log.info (
      workflow.success
        ? "\nDone! Purple ran successfully. See the results in: ${params.outdir}\n"
        : "\nOops .. something went wrong\n"
    )
}
