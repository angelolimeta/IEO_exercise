nextflow.enable.dsl=2

// Define parameters
params.data_dir = '/home/angelol/projects/IEO_exercise/data/exercise'
params.filtered_dir = '/home/angelol/projects/IEO_exercise/data/filtered_reads'
params.results_dir = '/home/angelol/projects/IEO_exercise/results'
params.kallisto_index = '/home/angelol/projects/IEO_exercise/data/grch38/transcriptome.idx'
params.bowtie2_index = '/home/angelol/projects/IEO_exercise/data/grch38/Bowtie2_index/GRCh38_noalt_as'
params.include_fastp = false

// Create output directories if they don't exist
if (!file(params.results_dir).exists()) {
    file(params.results_dir).mkdirs()
}

// Define input channel based on the filtered reads directory or raw data
samples = Channel
    .fromPath("${params.data_dir}/samples.csv")
    .splitCsv(header: true)
    .map { row -> tuple(row['Sample_ID'], row['Sample_Folder'], row['R1_Fastq'], row['R2_Fastq']) }

// Workflow definition
workflow {

    if (params.include_fastp) {
        // Run fastp if include_fastp is true
        fastp_input_ch = samples.map { sample_id, sample_folder, r1_fastq, r2_fastq ->
            tuple(
                sample_id,
                file("${params.data_dir}/${sample_folder}/${r1_fastq}"),
                file("${params.data_dir}/${sample_folder}/${r2_fastq}")
            )
        }

        fastp_output_ch = QualityFilterFastp(fastp_input_ch)

        // Split RNA-seq and ChIP-seq samples
        rna_samples_ch = fastp_output_ch.filter { sample_id, r1, r2 -> sample_id.contains('IDR') }
        chip_samples_ch = fastp_output_ch.filter { sample_id, r1, r2 -> sample_id.contains('IDC') || sample_id == 'S48549_IDC_393' }

    } else {
        // Start from filtered reads if include_fastp is false
        filtered_reads = Channel
            .fromPath("${params.filtered_dir}/*_R1.fastq.gz")
            .map { r1 ->
                def sample_id = r1.baseName.replace('_R1.fastq', '')
                def r2 = file(r1.toString().replace('_R1.fastq.gz', '_R2.fastq.gz'))
                tuple(sample_id, r1, r2)
            }

        // Split RNA-seq and ChIP-seq samples
        rna_samples_ch = filtered_reads.filter { sample_id, r1, r2 -> sample_id.contains('IDR') }
        chip_samples_ch = filtered_reads.filter { sample_id, r1, r2 -> sample_id.contains('IDC') || sample_id == 'S48549_IDC_393' }
    }

    // Process the control sample first
    control_bam_ch = Bowtie2AlignControl(chip_samples_ch.filter { sample_id, r1, r2 -> sample_id == 'S48549_IDC_393' })

    // Process remaining ChIP-seq samples (excluding control)
    bam_files_ch = Bowtie2AlignTreatment(chip_samples_ch.filter { sample_id, r1, r2 -> sample_id != 'S48549_IDC_393' })

    // Combine BAM files from treatment samples and control sample for MACS2
    macs2_input_ch = bam_files_ch.combine(control_bam_ch)
        .view { "Combined BAM files for MACS2: $it" }  // Add a view to check combined BAM files
    MACS2Callpeak(macs2_input_ch)

    // Run RNA-seq alignment with Kallisto
    KallistoQuant(rna_samples_ch)
}

// Process for fastp
process QualityFilterFastp {
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz")

    publishDir "${params.filtered_dir}", mode: 'copy'

    script:
    """
    fastp -i ${r1} \\
          -I ${r2} \\
          -o ${sample_id}_R1.fastq.gz \\
          -O ${sample_id}_R2.fastq.gz \\
          -j ${sample_id}_fastp.json \\
          -h ${sample_id}_fastp.html
    """
}

// Process for Kallisto quantification (RNA-seq)
process KallistoQuant {
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "kallisto/${sample_id}"

    publishDir "${params.results_dir}/rna_data", mode: 'copy'

    script:
    """
    mkdir -p kallisto/${sample_id}
    kallisto quant -i ${params.kallisto_index} -o kallisto/${sample_id} -b 100 -t ${task.cpus} ${r1} ${r2}
    """
}

// Process for Bowtie2 alignment (ChIP-seq) for control
process Bowtie2AlignControl {
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path("${sample_id}_sorted.bam")

    publishDir "${params.results_dir}/chip_data", mode: 'copy'

    script:
    """
    bowtie2 -x ${params.bowtie2_index} -1 ${r1} -2 ${r2} -S ${sample_id}.sam --threads ${task.cpus}
    samtools view -bS ${sample_id}.sam > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}

// Process for Bowtie2 alignment (ChIP-seq) for treatment samples
process Bowtie2AlignTreatment {
    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    publishDir "${params.results_dir}/chip_data", mode: 'copy'

    script:
    """
    bowtie2 -x ${params.bowtie2_index} -1 ${r1} -2 ${r2} -S ${sample_id}.sam --threads ${task.cpus}
    samtools view -bS ${sample_id}.sam > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}

// Process for MACS2 peak calling (ChIP-seq)
process MACS2Callpeak {
    input:
    tuple val(sample_id), path(treatment_bam), path(control_bam)

    output:
    path "macs2/*_peaks.narrowPeak"
    path "macs2/*_peaks.xls"
    path "macs2/*_summits.bed"

    publishDir "${params.results_dir}/chip_data", mode: 'copy'

    script:
    """
    echo "Running MACS2 with treatment: ${treatment_bam} and control: ${control_bam}"
    macs2 callpeak -t ${treatment_bam} -c ${control_bam} -f BAM -g hs -n macs2_output_${sample_id} --outdir macs2/
    """
}
