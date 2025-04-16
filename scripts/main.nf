#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process REMOVE_HOST_READS {
    publishDir "${params.outdir}/host_filtered_reads", mode: params.publish_dir_mode
    module "samtools-1.19/python-3.12.0"
    cpus 8

    input: 
    tuple val(meta), path(bamfile)

    output: 
    tuple val(meta), path("*.fastq.gz")

    script:
    """
    samtools fastq -f 4 \
    --threads 4 \
    $bamfile \
    -1 "${meta.sample_id}_R1.fastq.gz" -2 "${meta.sample_id}_R2.fastq.gz"
    """

}

process RUN_KRAKEN2 {
    publishDir "${params.outdir}/kraken2", mode: params.publish_dir_mode
    module "kraken2/2.1.2"
    input: 
    tuple val(meta), path(reads)
    path(REF_DIR)
    val(C_SCORE)

    output:
    tuple val(meta), path("*.kraken"), emit: report

    script: 
    """
    kraken2 --paired \
    --gzip-compressed \
    --use-names \
    --confidence ${C_SCORE} \
    --db ${REF_DIR} \
    --report ${meta.sample_id}_${C_SCORE}.kraken \
    --report-minimizer-data \
    --output /dev/null ${reads} 
    """

}

process KRAKEN_TOOLS {
    publishDir "${params.outdir}/krakentools", mode: params.publish_dir_mode
    module "krakentools/1.2.4"
    input: 
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*.mpa"), emit: mpa_file

    script: 
    """
    kreport2mpa.py \
    --report ${report} \
    --output ${meta.sample_id}.kraken.mpa
    """
    

}


workflow {
    bam_ch = Channel.fromPath(params.bamfiles, checkIfExists: true)
    | map {file -> tuple([sample_id: file.simpleName], file)}
    REMOVE_HOST_READS(bam_ch)
    reference = file(params.reference_db, checkIfExists: true)
    RUN_KRAKEN2(REMOVE_HOST_READS.out, reference, params.c_score)
    KRAKEN_TOOLS(RUN_KRAKEN2.out.report)

}