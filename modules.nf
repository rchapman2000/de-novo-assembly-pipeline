
// Builds a bowtie2 index for a provided reference file
process Index_Host_Reference {
    input:
        // Tuple contains the reference name and reference file.
        tuple val(refName), file(ref)
        // The name of the output directory
        val outDir
        // The number of threads provided
        val threads

    output:
        // Tuple containing the reference index directory and the index file name.
        tuple file("ref-idx/"), val(refName)

    publishDir "${outDir}", mode: 'copy'

    /*
    Creates an index Directory and then runs bowtie2-build to create the index
    */
    script:
    """
    #!/bin/bash
    
    mkdir ref-idx/

    bowtie2-build --threads ${threads} ${ref} ref-idx/${refName}given 
    """
}

// Creates a fastqc report for a set of paired-end reads provided.
process QC_Report {
    input:
        // Tuple contains the file basename as well as the paired-end read files
        tuple val(base), file(F1), file(F2)
        // The name of the output directory
        val outDir
        // The name of the directory to place the fastqc output into (allows
        // for this command to be used multiple times and to separate the output
        // i.e. pre-processed-reads vs trimmed-reads.)
        val dirName
        // The number of threads provided.
        val threads

    output:
        // The directory cotaining the fastqc files produced.
        file("${base}")
    
    publishDir "${outDir}${dirName}/", mode: 'copy'

    // Creates a Directory and then runs fastqc on a set of reads.
    script:
    """
    #!/bin/bash

    mkdir ${base}

    fastqc ${F1} ${F2} -o ${base} --threads ${threads}
    """
}


// Performs quality and adapter trimming on a set of paired-end reads.
process Trimming {
    input: 
        // Tuple cotains the file basename as well as the paired-end read files.
        tuple val(base), file(R1), file(R2)
        // The name of the output directory.
        val outDir
        // The adapter file in fasta format.
        file adapters

    output:
        // Tuple containing the file basename and the trimmed forward and reverse reads
        tuple val(base), file("${base}_1.trimmed.fastq"), file("${base}_2.trimmed.fastq")
        // Tuple containing the unpaired read files.
        tuple file("${base}_1.unpaired.fastq.gz"), file("${base}_2.unpaired.fastq.gz") 

    publishDir "${outDir}${base}-Processed-Reads", mode: 'copy'

    script:
    /*
    Uses trimmomatic to trim the read files.
    Trimmomatic performs trimming steps in the order provided in the 
    command line. Our steps:
    1. ILLUMINACLIP: Removes illumina adapters
    2. LEADING: Begins at the start of the read and trims bases with quality less than
        5 until it hits a base above that threshold.
    3. TRAILING: Begins at the end of the read and trims bases with quality less than 5
        until it hits a base above that threshold.
    4. SLIDINGWINDOW: Begins at the 5' end and scans in 4 base windows, trimming when it hits
        an average quality less than 20.
    5. MINLEN: Removes reads less than 75 bp long.

    Prinseq (used in the following process) cannot take gzipped reads as input. Thus,
    only the unpaired reads are gzipped to save space.
    */
    """
    #!/bin/bash

    trimmomatic PE ${R1} ${R2} ${base}_1.trimmed.fastq ${base}_1.unpaired.fastq \
    ${base}_2.trimmed.fastq ${base}_2.unpaired.fastq ILLUMINACLIP:${adapters}:2:30:10:1:true \
    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:75

    gzip ${base}_1.unpaired.fastq ${base}_2.unpaired.fastq
    """
        
}

// Removes PCR Duplicates from a set of reads using prinseq.
process Remove_PCR_Duplicates {
    input:
        // Tuple contains the file basename and the two read files.
        tuple val(base), file(R1), file(R2)
        // The name of the output directory.
        val outDir
    output:
        // Tuple containing the file basename and the deduped, paired-end reads.
        tuple val(base), file("${base}_deduped_1.fastq.gz"), file("${base}_deduped_2.fastq.gz")
        // The log file produced by prinseq for troubleshooting.
        file "${base}-prinseq-log.txt"

    publishDir "${outDir}${base}-Processed-Reads", mode: 'copy'

    /*
    Calls prinseq on the paired-end reads. The --derep command
    tells prinseq to deduplicate the reads. Prinseq has 5 
    derep modes:

    "1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 4
    (reverse complement exact duplicate), 5 (reverse complement
    5'/3' duplicate))"

    This pipelines only removes exact duplicates. (--derep 1)

    The fastq files produced are gzipped to save space.
    */
    script:
    """
    #!/bin/bash

    prinseq-lite.pl -fastq ${R1} -fastq2 ${R2} --out_good ${base}_deduped --derep 1 -log ${base}-prinseq-log.txt

    gzip ${base}_deduped_1.fastq ${base}_deduped_2.fastq
    """
}

// Removes host reads by aligning to a host genome using bowtie2.
process Host_Read_Removal {
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), file(R1), file(R2)
        // The output directory name
        val outDir
        // Tuple contains the bt2 index directory and basename of the index files.
        tuple file(refDir), val(refName)
        // The number of threads provided.
        val threads

    output:
        // Tuple contains the file basename and the paired-end read files with host reads removed.
        tuple val(base), file("${base}_host_removed_1.fq.gz"), file("${base}_host_removed_2.fq.gz")
        // A directory containing the alignment to the host file.
        file "host-reads/${base}-host.sam"
    
    publishDir "${outDir}${base}-Processed-Reads", mode: 'copy'

    script:
    /*
    Aligns the reads to a host reference genome using bowtie2. Local alignment is used (--local)
    to ensure the reads are aligne without gaps. The unaligned reads are sent back ot into a
    fastq files (--un-conc).

    The unaligned read files are then renamed to give them the typical paired end
    read name scheme.

    Finally, the reads are gzipped to preserve space.
    */
    """
    #!/bin/bash
    mkdir host-reads
    bowtie2 --threads ${threads} -x ${refDir}/${refName} -1 ${R1} -2 ${R2} --local -S host-reads/${base}-host.sam --un-conc ${base}_host_removed

    mv ${base}_host_removed.1 ${base}_host_removed_1.fq
    mv ${base}_host_removed.2 ${base}_host_removed_2.fq

    gzip ${base}_host_removed_1.fq ${base}_host_removed_2.fq
    """
}

// Uses spades to produce a de novo assembly.
process Assembly {
    input:
        // Tuple contains the file basename and paired-end reads.
        tuple val(base), file(R1), file(R2)
        // The name of the output directory.
        val outDir
        // The number of threads provided.
        val threads
    
    output:
        // Tuple contains the basename of the sample and the assembly directory 
        // produced by spades.
        tuple val(base), file("${base}-Assembly")
    
    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Runs spades using the provided paired-end reads.
    */
    """
    #!/bin/bash

    spades.py --threads ${threads} -1 ${R1} -2 ${R2} -o ${base}-Assembly
    """
}

// Generates an alignment of the asembly contigs to a reference genome
// using minimap.
process Contig_Alignment {
    input:
        // Tuple contains the sample basename
        // and the assembly directory
        tuple val(base), file(assemblyDir)
        // The name of the output directory
        val outDir
        // the reference fasta file to be aligned to.
        file ref

    output:
        // Tuple contains the file basename and the alignment bam file.
        tuple val(base), file("${base}-contig-align.bam")

    publishDir "${outDir}${assemblyDir}", mode: 'copy'
    
    script:
    /*
    Uses minimap2 to align the contigs to the reference fasta.

    Then uses samtools to conver the alignment sam into bam format
    (samtools view). The bam format is then sorted and stored in a bam
    file (samtools sort).
    */
    """
    #!/bin/bash

    minimap2 -ax asm5 ${ref} ${assemblyDir}/contigs.fasta > align.sam

    samtools view -b align.sam | samtools sort > ${base}-contig-align.bam
    """
}