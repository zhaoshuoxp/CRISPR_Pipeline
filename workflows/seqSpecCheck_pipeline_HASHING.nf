nextflow.enable.dsl=2

include { seqSpecCheck } from '../processes/seqSpecCheck.nf'

workflow guideWorkflow {

    Channel
        .of(params.fastq_files_guide[0])
        .set { raw_input }

    raw_input
        .map { line -> 
            def parts = line.trim().split(/\s+/)  
            if (parts.size() >= 2) {
                return [parts[0], parts[1]]  
            } else {
                return []
            }
        }
        .filter { parts -> parts.size() == 2 }  
        .map { r1_path, r2_path -> [file(r1_path), file(r2_path)] }
        .set { guide_fastq_ch }

    guide_fastq_r1_ch = guide_fastq_ch.map { pair -> pair[0] }
    guide_fastq_r2_ch = guide_fastq_ch.map { pair -> pair[1] }

    guide_fastq_r1_ch.view { "R1 file: $it" }
    guide_fastq_r2_ch.view { "R2 file: $it" }
    
    guide_seqSpecCheck = seqSpecCheck(
        guide_fastq_r1_ch,
        guide_fastq_r2_ch,
        file(params.guide_metadata),
        'guide'
    )

    emit:
    guide_seqSpecCheck_plots = guide_seqSpecCheck.seqSpecCheck_plots
    guide_position_table = guide_seqSpecCheck.position_table
}

workflow hashingWorkflow {

    Channel
        .of(params.fastq_files_hashing[0])
        .set { raw_input }

    raw_input
        .map { line -> 
            def parts = line.trim().split(/\s+/)  
            if (parts.size() >= 2) {
                return [parts[0], parts[1]]  
            } else {
                return []
            }
        }
        .filter { parts -> parts.size() == 2 }  
        .map { r1_path, r2_path -> [file(r1_path), file(r2_path)] }
        .set { hashing_fastq_ch }

    hashing_fastq_r1_ch = hashing_fastq_ch.map { pair -> pair[0] }
    hashing_fastq_r2_ch = hashing_fastq_ch.map { pair -> pair[1] }

    hashing_fastq_r1_ch.view { "R1 file: $it" }
    hashing_fastq_r2_ch.view { "R2 file: $it" }

    hashing_seqSpecCheck = seqSpecCheck(
        hashing_fastq_r1_ch,
        hashing_fastq_r2_ch,
        file(params.hashing_metadata),
        'hashing'
    )

    emit:
    hashing_seqSpecCheck_plots = hashing_seqSpecCheck.seqSpecCheck_plots
    hashing_position_table = hashing_seqSpecCheck.position_table
}

workflow seqSpecCheck_pipeline_HASHING {
    main:
    guide = guideWorkflow()
    hashing = hashingWorkflow()

    emit:
    guide_seqSpecCheck_plots = guide.guide_seqSpecCheck_plots
    guide_position_table = guide.guide_position_table
    hashing_seqSpecCheck_plots = hashing.hashing_seqSpecCheck_plots
    hashing_position_table = hashing.hashing_position_table
}
