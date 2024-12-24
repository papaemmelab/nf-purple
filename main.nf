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
    normal        : ${params.normal}
    normalBam     : ${params.normalBam}
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
        PURPLE is a purity ploidy estimator primarily designed for whole genome 
        sequenced (WGS) data. It combines B-allele frequency (BAF) from AMBER, 
        read depth ratios from COBALT, somatic variants and structural variants 
        to estimate the purity and copy number profile of a tumor sample. PURPLE 
        supports both grch 37 and 38 reference assemblies.

        For more info please see: 
        https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md
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
