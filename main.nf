#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
De Novo Assembly Pipeline:
    
Takes an input of paired-end fastq files from illumina
sequencers, processes the reads, removes host contamination,
and performs a de novo assembly.

USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR

OPTIONS:

--input INPUT_DIR - [Required] A directory containing paired-end fastq files

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

OPTIONAL:

    --host_fasta HOST_FASTA - A host fasta file that the reads will be aligned to to remove host contamination.
            
    --host_bt2_index INDEX_DIRECTORY - A directory containing an existing bowtie2 index to be used for host read removal. Must be in its own directory. If using a large genome, this option will greatly improve the pipeline runtime.

    --unicycler - Assembles reads using unicycler in place of spades

    --threads INT - The number of threads that can be use to run pipeline tools in parallel

    --ref REFERENCE_FASTA - The pipeline will align contigs produced by assembly to this reference

    --minLen INT - The minimum length of a read to keep post trimming [Default = 75bp]

    --minTrimQual INT - The average basecall quality threshold below which to trim a read. During trimming, trimmomatic performs a sliding window checking the average base quality, and removing the rest of the read if it drops below this treshold. [Default = 20]
    """
}

// Function that checks whether a directory ends in a trailing slash.
// This is useful for directory variables that are not parsed into
// file objects in the pipeline (such as the output directory).
def checkDirectoryEnding (fileName) {
    // Grabs the last character in the directory name.
    lastChar = fileName.substring(fileName.length() - 1, fileName.length())

    // Checks whether the last character is slash
    if (lastChar != "/") {
        // If it is not a slash, add that to the directory name.
        fileName = fileName + "/"
    }
    
    // Return the directory name.
    return fileName
}

// Function creates the header for the summary file based on the parameters
// supplied by the user. Because the pipeline dynamically changes based on
// what the user specifies, the summary file must also be alter to reflect
// the analysis. This will also make incorporating new modules easier.
def createSummaryHeader (hostRef, hostIdx) {
    
    // The header will always start with the sample.
    FinalHeader = 'Sample,'

    // The summary sheet will also always contain the Raw and trimmed read counts.
    FinalHeader = FinalHeader + "Raw Reads,Trimmed Reads,Deduped Reads,"

    // Next, the user may supply host removal, which would be the next useful
    // statistic to know. Thus, if a host reference or bowtie2 index is supplied,
    // add this field to the header.
    if (hostRef != false || hostIdx != false) {
        FinalHeader = FinalHeader + "Non-Host Reads,"
    }


    // Finally, the pipeline will always report the number of contigs and scaffolds
    FinalHeader = FinalHeader + "Contigs Generated,Scaffolds Generated"

    return FinalHeader
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// Defines input parameters. Setting to false by default
// allows us to check that these have been set by the user.
params.input = false
params.host_fasta = false
params.host_bt2_index = false
params.ref = false
params.output = false
params.threads = 1
params.unicycler = false
params.minLen = 75
params.minTrimQual = 20

// Inports modules
include { Setup } from "./modules.nf"
include { Index_Host_Reference } from './modules.nf'
include { QC_Report } from './modules.nf'
// The same module cannot be used more than once,
// thus it is aliased to be used multiple times.
include { QC_Report as QC_Report_Trimmed } from './modules.nf'
include { QC_Report as QC_Report_Deduped } from './modules.nf'
include { QC_Report as QC_Report_Host_Removed } from './modules.nf'
include { Trimming } from './modules.nf'
include { Remove_PCR_Duplicates } from './modules.nf'
include { Host_Read_Removal } from './modules.nf'
include { Spades_Assembly } from './modules.nf'
include { Unicycler_Assembly } from './modules.nf'
include { Contig_Alignment } from './modules.nf'
include { Write_Summary } from './modules.nf'

adapters = file("${baseDir}/adapters.fa")

// Checks the input parameter
inDir = ""
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}
else {
    inDir = file(params.input).toString()
}

println "Input Directory: ${inDir}"

// Create a channel for hte input files.
inputFiles_ch = Channel
    // Pull from pairs of files (illumina fastq files denoted by having R1 or R2 in
    // the file name).
    .fromFilePairs("${inDir}/*_R{1,2}*.fastq*")
    // The .fromFilePairs() function spits out a list where the first 
    // item is the base file name, and the second is a list of the files.
    // This command creates a tuple with the base file name and two files.
    .map { it -> [it[0], it[1][0], it[1][1]]}

// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = file(params.output).toString()
}
println outDir

// Parses the host options (--host_fasta and --host_bt2_index).
// For this, we cannot use a channel like we did for the input files
// because then Nextfow would only run other modules once.
// Thus, we need to manually create a tuple of input data to pass to the indexing
// step and the alignment step.
hostRefData = ''
hostRefIdxData = ''
hostRefName = 'NONE'
if (params.host_fasta != false && params.host_bt2_index != false) {
    // If both options are supplied, notify the user and exit.
    println "ERROR: you have specified both a host fasta file and bowtie2 index. Please only supply one."
    exit(1)
}
else if (params.host_fasta != false) {
    if (!(file(params.host_fasta).exists())) {
        // If the file supplied does not exist, notify the user and exit.
        println "ERROR: ${params.host_fasta} does not exist."
        exit(1)
    }
    else {
        // Parse the file into a file object
        hostRef = file(params.host_fasta)
        // Use the getBaseName() function to 
        // get the reference name. This will be
        // used to name the bowtie2 index.
        hostRefName = hostRef.getBaseName()
        // Place these both into a tuple.
        hostRefData = tuple(hostRefName, hostRef)
    }
}
// If the user supplied the --host_bt2_index
else if (params.host_bt2_index != false) {
    if (!(file(params.host_bt2_index).exists())) {
        // If the index provided does not exist, notify the user and exit.
        println "Error: ${params.host_bt2_idx} does not exist."
        exit(1)
    }
    else {
        // Parse the directory into a file object
        hostRefDir = file(params.host_bt2_index)
        println hostRefDir
        // Grab a list of file objects from the directory
        // ending in .bt2
        indexFiles = file("${hostRefDir}/*.bt2")
        if (indexFiles.size() == 0){
            // If there are no file in the directory ending in bt2, notify the user and exit
            println "Index Directory provided (${params.host_bt2_index}) does not contain any bt2 files"
            exit(1)
        }
        else {
            // Use the getSimpleName() function to grab the base name
            // of the index files (getSimpleName() removes anything following
            // the last . in a file name.)
            hostRefName = indexFiles[0].getSimpleName()
            println hostRefName
            // Place the index dir and name into a tuple.
            hostRefIdxData = tuple(hostRefDir, hostRefName)
        }
    }
}


// Parses the ref option.
refFile = ''
if (params.ref != false) {
    if (!(file(params.ref).exists())) {
        // If the reference file did not exist, notify the user and exit.
        println "ERROR: ${params.ref} does not exist."
        exit(1)
    }
    else {
        // Parse the provided file into a file object.
        refFile = file(params.ref)
    }
}

assembler = 'SPADES'
if (params.unicycler != false) {
    assembler = "Unicycler"
}

summaryHeader = createSummaryHeader(params.host_fasta, params.host_bt2_index)

workflow {

    Setup( summaryHeader, hostRefName, params.minLen, params.minTrimQual, assembler, outDir )

    // Use FASTQC to perform an initial QC check on the reads
    QC_Report( inputFiles_ch, outDir, "FASTQC-Pre-Processing", params.threads )

    // Perform adapter and quality trimming with Trimmomatic.
    Trimming( inputFiles_ch, outDir, adapters, params.minLen, params.minTrimQual)

    // Use FASTQC to perform a QC check on the trimmed reads.
    QC_Report_Trimmed( Trimming.out[0], outDir, "FASTQC-Trimmed", params.threads )

    // Perform PCR Duplicate removal using prinseq.
    Remove_PCR_Duplicates( Trimming.out[0], outDir, Trimming.out[2] )

    // Use FASTQC to perform a QC check on the deduped reads.
    QC_Report_Deduped( Remove_PCR_Duplicates.out[0], outDir, "FASTQC-Deduplicated", params.threads )

    // If the user supplied a reference fasta file, build a bowtie2 index and use
    // that for alignment.
    if (params.host_fasta) {
        Index_Host_Reference( hostRefData, outDir, params.threads )
        Host_Read_Removal( Remove_PCR_Duplicates.out[0], outDir, Index_Host_Reference.out, params.threads, Remove_PCR_Duplicates.out[2] )
        QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQC-Host-Removed", params.threads )

        // If the user supplied the --unicycler parameter, use unicycler for
        // assembly
        if (params.unicycler != false) {
            // Performs de novo assembly using unicycler
            Unicycler_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2])
            
            if (params.ref != false) {
                // Align the contigs to a reference genome using minimap2 and samtools
                Contig_Alignment( Unicycler_Assembly.out[0], outDir, refFile )
            }

            Write_Summary(Unicycler_Assembly.out[1], outDir )
        }
        else {
            // Perform de novo assembly using spades.
            Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2] )
            if (params.ref != false) {
                // Align the contigs to a reference genome using minimap2 and samtools
                Contig_Alignment( Spades_Assembly.out[0], outDir, refFile )
            }

            Write_Summary( Spades_Assembly.out[2], outDir )
        }
    }
    // If the user supplied an existing bowtie2 index, use that for alignment.
    else if (params.host_bt2_index) {
        Host_Read_Removal( Remove_PCR_Duplicates.out[0], outDir, hostRefIdxData, params.threads, Remove_PCR_Duplicates.out[2] )
        QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQC-Host-Removed", params.threads )

        // If the user supplied the --unicycler parameter, use unicycler for
        // assembly
        if (params.unicycler != false) {
            // Performs de novo assembly using unicycler
            Unicycler_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2])
            
            if (params.ref != false) {
                // Align the contigs to a reference genome using minimap2 and samtools
                Contig_Alignment( Unicycler_Assembly.out[0], outDir, refFile )
            }

            Write_Summary(Unicycler_Assembly.out[1], outDir )
        }
        else {
            // Perform de novo assembly using spades.
            Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, Host_Read_Removal.out[2] )
            if (params.ref != false) {
                // Align the contigs to a reference genome using minimap2 and samtools
                Contig_Alignment( Spades_Assembly.out[0], outDir, refFile )
            }

            Write_Summary( Spades_Assembly.out[2], outDir )
        }
    }
    else {

        // If the user supplied the --unicycler parameter, use unicycler for
        // assembly
        if (params.unicycler != false) {
            // Performs de novo assembly using unicycler
            Unicycler_Assembly( Remove_PCR_Duplicates.out[0], outDir, params.threads, Remove_PCR_Duplicates.out[2])
            
            if (params.ref != false) {
                // Align the contigs to a reference genome using minimap2 and samtools
                Contig_Alignment( Unicycler_Assembly.out[0], outDir, refFile )
            }

            Write_Summary(Unicycler_Assembly.out[1], outDir )
        }
        else {
            // Perform de novo assembly using spades.
            Spades_Assembly( Remove_PCR_Duplicates.out[0], outDir, params.threads, Remove_PCR_Duplicates.out[2] )
            if (params.ref != false) {
                // Align the contigs to a reference genome using minimap2 and samtools
                Contig_Alignment( Spades_Assembly.out[0], outDir, refFile )
            }

            Write_Summary( Spades_Assembly.out[2], outDir )
        }
    }

}