// Creates a parameters file and a summary file to 
// be added to later
process Setup {
    input:
        // The name of the reference supplied (for use
        // in the parameters file)
        val refName
        // The minimum read length allowed post trimmming (for
        // use in the parameters file)
        val minLen
        // The name of the assembler used (for use in
        // the parameters file)
        val assembler
        // The output directory to be used.
        file outDir
        
    output:
        // The parameters file created.
        file "analysis-parameters.txt"
        // The blank summary file to be added to.
        file "stats-summary.csv"

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Creates a parameters file (in case the user needs to look back at how they ran the analysis)
    as well as a blank summary file, which will later be populated with data from each
    sample.

    The parameters file contains:
        1. The name of the reference supplied
        2. The minimum read length allowed after trimmming
        3. The assembler used

    The summary file contains:
        1. The sample
        2. Raw Read counts for R1 and R2
        3. Reads after trimming for R1 and R2
        4. Reads after deduplication for R1 and R2
        5. Reads after host removal for R1 and R2
        6. Number of contigs
        7. Number of scaffolds
    */
    """
    #!/bin/bash

    touch analysis-parameters.txt

    echo "Host Genome : ${refName}" >> analysis-parameters.txt
    echo "Minimum Read Length Allowed : ${minLen} bp" >> analysis-parameters.txt
    echo "Assember : ${assembler}" >> analysis-parameters.txt

    touch stats-summary.csv

    echo "Sample,Raw Reads R1,Trimmed Reads R1,Deduped Reads R1,Non-Host Reads R1,Raw Reads R2,Trimmed Reads R2,Deduped Reads R2,Non-Host Reads R2,Contigs Generated,Scaffolds Generated" > stats-summary.csv
    """
}


// Builds a bowtie2 index for a provided reference file
process Index_Host_Reference {
    input:
        // Tuple contains the reference name and reference file.
        tuple val(refName), file(ref)
        // The output directory
        file outDir
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

    bowtie2-build --threads ${threads} ${ref} ref-idx/${refName}
    """
}

// Creates a fastqc report for a set of paired-end reads provided.
process QC_Report {
    input:
        // Tuple contains the file basename as well as the paired-end read files
        tuple val(base), file(F1), file(F2)
        // The output directory
        file outDir
        // The name of the directory to place the fastqc output into (allows
        // for this command to be used multiple times and to separate the output
        // i.e. pre-processed-reads vs trimmed-reads.)
        val dirName
        // The number of threads provided.
        val threads

    output:
        // The directory cotaining the fastqc files produced.
        file("${base}")
    
    publishDir "${outDir}/${dirName}/", mode: 'copy'

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
        // The output directory
        file outDir
        // The adapter file in fasta format.
        file adapters
        //The minimum read lenths to be allowed post trimming
        val minLen

    output:
        // Tuple containing the file basename and the trimmed forward and reverse reads
        tuple val(base), file("${base}_1.trimmed.fastq"), file("${base}_2.trimmed.fastq")
        // Tuple containing the unpaired read files.
        tuple file("${base}_1.unpaired.fastq.gz"), file("${base}_2.unpaired.fastq.gz") 
        // Tuple containing the raw and trimmed read counts for the paired-end sample.
        tuple env(raw_reads_1), env(raw_reads_2), env(trimmed_reads_1), env(trimmed_reads_2)

    publishDir "${outDir}/${base}-Processed-Reads", mode: 'copy'

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

    Additionally, a count of the raw and trimmed reads are taking
    */
    """
    #!/bin/bash

    raw_reads_1=\$((\$(gunzip -c ${R1} | wc -l)/4))
    raw_reads_2=\$((\$(gunzip -c ${R1} | wc -l)/4))

    trimmomatic PE ${R1} ${R2} ${base}_1.trimmed.fastq ${base}_1.unpaired.fastq \
    ${base}_2.trimmed.fastq ${base}_2.unpaired.fastq ILLUMINACLIP:${adapters}:2:30:10:1:true \
    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:${minLen}

    gzip ${base}_1.unpaired.fastq ${base}_2.unpaired.fastq
    
    trimmed_reads_1=\$((\$(cat ${base}_1.trimmed.fastq | wc -l)/4))
    trimmed_reads_2=\$((\$(cat ${base}_2.trimmed.fastq | wc -l)/4))
    """
        
}

// Removes PCR Duplicates from a set of reads using prinseq.
process Remove_PCR_Duplicates {
    input:
        // Tuple contains the file basename and the two read files.
        tuple val(base), file(R1), file(R2)
        // The output directory
        file outDir
    output:
        // Tuple containing the file basename and the deduped, paired-end reads.
        tuple val(base), file("${base}_deduped_1.fastq.gz"), file("${base}_deduped_2.fastq.gz")
        // The log file produced by prinseq for troubleshooting.
        file "${base}-prinseq-log.txt"
        // Tuple containing the number of reads after deduplication for R1 and R2
        tuple env(deduped_reads_1), env(deduped_reads_2)

    publishDir "${outDir}/${base}-Processed-Reads", mode: 'copy'

    /*
    Calls prinseq on the paired-end reads. The --derep command
    tells prinseq to deduplicate the reads. Prinseq has 5 
    derep modes:

    "1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 4
    (reverse complement exact duplicate), 5 (reverse complement
    5'/3' duplicate))"

    This pipelines only removes exact duplicates. (--derep 1)

    The fastq files produced are gzipped to save space.

    Counts of the reads after deduplication are also taken.
    */
    script:
    """
    #!/bin/bash

    prinseq-lite.pl -fastq ${R1} -fastq2 ${R2} --out_good ${base}_deduped --derep 1 -log ${base}-prinseq-log.txt

    gzip ${R1} ${R2} ${base}_deduped_1.fastq ${base}_deduped_2.fastq

    deduped_reads_1=\$((\$(gunzip -c ${base}_deduped_1.fastq | wc -l)/4))
    deduped_reads_2=\$((\$(gunzip -c ${base}_deduped_2.fastq | wc -l)/4))
    """
}

// Removes host reads by aligning to a host genome using bowtie2.
process Host_Read_Removal {
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), file(R1), file(R2)
        // The output directory
        file outDir
        // Tuple contains the bt2 index directory and basename of the index files.
        tuple file(refDir), val(refName)
        // The number of threads provided.
        val threads

    output:
        // Tuple contains the file basename and the paired-end read files with host reads removed.
        tuple val(base), file("${base}_host_removed_1.fq.gz"), file("${base}_host_removed_2.fq.gz")
        // A directory containing the alignment to the host file.
        file "host-reads/${base}-host.sam"
        // Tuple contains the number of reads after host removal for R1 and R2
        tuple env(nonHost_reads_1), env(nonHost_reads_2)
    
    publishDir "${outDir}/${base}-Processed-Reads", mode: 'copy'

    script:
    /*
    Aligns the reads to a host reference genome using bowtie2. Local alignment is used (--local)
    to ensure the reads are aligne without gaps. The unaligned reads are sent back ot into a
    fastq files (--un-conc).

    The unaligned read files are then renamed to give them the typical paired end
    read name scheme.

    Additionally, the reads are gzipped to preserve space.

    Finally, the number of reads after host removal are recorded.
    */
    """
    #!/bin/bash
    mkdir host-reads
    bowtie2 --threads ${threads} -x ${refDir}/${refName} -1 ${R1} -2 ${R2} --local -S host-reads/${base}-host.sam --un-conc ${base}_host_removed

    mv ${base}_host_removed.1 ${base}_host_removed_1.fq
    mv ${base}_host_removed.2 ${base}_host_removed_2.fq

    gzip ${base}_host_removed_1.fq ${base}_host_removed_2.fq

    nonHost_reads_1=\$((\$(gunzip -c ${base}_deduped_1.fastq | wc -l)/4))
    nonHost_reads_2=\$((\$(gunzip -c ${base}_deduped_2.fastq | wc -l)/4))
    """
}

// Uses spades to produce a de novo assembly.
process Spades_Assembly {
    input:
        // Tuple contains the file basename and paired-end reads.
        tuple val(base), file(R1), file(R2)
        // The output directory
        file outDir
        // The number of threads provided.
        val threads
    
    output:
        // Tuple contains the basename of the sample and the assembled contigs 
        // produced by spades.
        tuple val(base), file("${base}-contigs.fasta")
        // The scaffolds produced by spades
        file "${base}-scaffolds.fasta"
        // Tuple contains the number of contigs and number of scaffolds
        // produced by spades.
        tuple env(numContigs), env(numScaffolds)
    
    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Runs spades using the provided paired-end reads.

    The contigs and scaffolds are renamed and moved, and
    the number of contigs/scaffolds are recorded.
    */
    """
    #!/bin/bash

    spades.py --threads ${threads} -1 ${R1} -2 ${R2} -o ${base}-Assembly

    mv ${base}-Assembly/contigs.fasta ./${base}-contigs.fasta

    numContigs=\$(grep ">" ${base}-contigs.fasta | wc -l)

    mv ${base}-Assembly/scaffolds.fasta ./${base}-scaffolds.fasta

    numScaffolds=\$(grep ">" ${base}-scaffolds.fasta | wc -l)
    """
}

process Unicycler_Assembly {
    input:
        // Tuple contains the file basename and paired-end reads.
        tuple val(base), file(R1), file(R2)
        // The output directory
        file outDir
        // The number of threads provided.
        val threads
    
    output:
        // Tuple contains the basename of the sample and the assembly directory 
        // produced by spades.
        tuple val(base), file("${base}-assembly.fasta")
        // Tuple contains the number of contigs and number of scaffolds
        // produced by spades.
        tuple env(numContigs), env(numScaffolds)
    
    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Runs unicycler using the provided paired-end reads.

    Unicycler only produces scaffolds, thus this file is moved/renamed,
    the number of scaffolds is recorded, and the number of contigs is set
    to a message letting the user know that these are not produced.
    */
    """
    #!/bin/bash

    unicycler --threads ${threads} -1 ${R1} -2 ${R2} -o ${base}-Assembly

    mv ${base}-Assembly/assembly.fasta ./${base}-assembly.fasta
    
    numScaffolds=\$(grep ">" ${base}-assembly.fasta | wc -l)

    numContigs="Unicycler only produces contigs"
    """
}

// Generates an alignment of the asembly contigs to a reference genome
// using minimap.
process Contig_Alignment {
    input:
        // Tuple contains the sample basename
        // and the assembly fasta to align
        tuple val(base), file(assembly)
        // The output directory
        file outDir
        // the reference fasta file to be aligned to.
        file ref

    output:
        // Tuple contains the file basename and the alignment bam file.
        tuple val(base), file("${base}-contig-align.bam")

    publishDir "${outDir}", mode: 'copy'
    
    script:
    /*
    Uses minimap2 to align the contigs to the reference fasta.

    Then uses samtools to conver the alignment sam into bam format
    (samtools view). The bam format is then sorted and stored in a bam
    file (samtools sort).
    */
    """
    #!/bin/bash

    minimap2 -ax asm5 ${ref} ${assembly} > align.sam

    samtools view -b align.sam | samtools sort > ${base}-contig-align.bam
    """
}

// Writes a line to the summary file for the sample.
process Write_Summary {
    input:
        // Tuple contains the sample basename and forward/reverse reads (the basename
        // is the only value important to this function).
        tuple val(base), file(F1), file(F2)
        // The summary file created earlier.
        file summary
        // Tuple contains raw read and trimmed read counts for R1 and R2
        tuple val(raw_reads_1), val(raw_reads_2), val(trimmed_reads_1), val(trimmed_reads_2)
        // Tuple contains read counts after deduplication for R1 and R2
        tuple val(deduped_reads_1), val(deduped_reads_2)
        // Tuple contains read counts after host removal for R1 and R2
        tuple val(nonHost_reads_1), val(nonHost_reads_2)
        // Tuple contains the number of contigs and number of scaffolds.
        tuple val(numContigs), val(numScaffolds)
        // The output directory.
        file outDir

    script:
    /*
    The values provided are written to the summary file.
    */
    """
    #!/bin/bash

    echo "${base},${raw_reads_1},${trimmed_reads_1},${deduped_reads_1},${nonHost_reads_1},${raw_reads_2},${trimmed_reads_2},${deduped_reads_2},${nonHost_reads_2},${numContigs},${numScaffolds}" >> ${outDir}/stats-summary.csv
    """  

}