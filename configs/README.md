# Pipeline Configuration Documentation

> **IMPORTANT**: All parameters defined in this configuration must have valid, non-empty values. The Nextflow pipeline requires every parameter to be properly set - empty or missing parameters will cause pipeline execution failures. Please ensure all parameters are specified before running the pipeline.

This document explains the configuration parameters for the single-cell RNA sequencing and CRISPR screening analysis pipeline.

## Input Data Parameters

### Metadata Paths
- `user_inference`: Path to CSV file containing user-defined cell pairs for testing
- `guide_metadata`: Path to TSV file containing guide RNA metadata
- `hashing_metadata`: Path to TSV file containing cell hashing metadata

### Dataset Configuration
- `DATASET_HASHING`: Boolean flag to enable/disable dataset hashing ('true'/'false')
- `transcriptome`: Species specification for transcriptome analysis ('human' or 'mouse')
- `seqspecs_directory`: Directory containing all YAML specification files and the cell barcode file(txt).
- `genome_download_path`: Download URL for human genome (hg38) reference
- `genome_local_path`: Local path for custom genome file
- `gtf_download_path`: Download URL for GENCODE gene annotations (v46)
- `gtf_local_path`:  Local path for custom GTF file

### Sequence Specification Files
- `scRNA_seqspec_yaml`: YAML file specifying single-cell RNA sequencing parameters
- `Guides_seqspec_yaml`: YAML file specifying guide RNA sequencing parameters
- `Hash_seqspec_yaml`: YAML file specifying cell hashing parameters

## Analysis Parameters

### scRNA Preprocessing Configuration
- `min_genes`: Minimum number of genes required per cell 
- `min_cells`: Minimum number of cells required
- `pct_mito`: Maximum percentage of mitochondrial genes allowed 

### Guide Assignment Configuration
- `assignment_method`: Method for guide-cell assignment ('sceptre' or 'cleanser')
- `capture_method`: Capture Method ('direct capture' or 'CROP-seq')
- `THRESHOLD`: Numerical threshold for assignments (default: 1)

### Perturbation Inference Configuration
- `inference_method`: Methods for statistical inference. 
   - 'sceptre,perturbo': select both methods
   - 'sceptre' or 'perturbo': Select either method individually

- `inference_option`: Type of inference analysis         
   - 'predefined_pairs': Uses a custom pairs_to_test.csv file from 'user_inference' for perturbation inference
   - 'by_distance': Generates target-guide pairs based on a defined distance threshold
   - 'all_by_all': Enables inference for all targets and all guides
- `distance_from_center`: Distance threshold for the 'by_distance' strategy to select targets for pairs_to_test (default: 1,000,000)

The following parameters apply to sceptre and can be adjusted accordingly:
- `moi`: Multiplicity of infection setting ('undecided', 'high' or 'low')
- `side`: Direction of statistical tests ('both')
- `grna_integration_strategy`: Strategy for combining guide information ('union')
- `resampling_approximation`: Distribution for resampling ('skew_normal')
- `control_group`: Control group specification ('default')
- `resampling_mechanism`: Mechanism for resampling ('default')
- `formula_object`: Statistical formula specification ('default')

### Evaluation Configuration

- `user_central_nodes`: User-defined central nodes 
   - 'undefined': Uses 'central_nodes_num' to determine central nodes.
   - Otherwise, specify the target names of interest.
- `central_nodes_num`: Number of central nodes to include in the network analysis (default: 2)

### FASTQ Files Configuration
The pipeline requires full paths to FASTQ files for three modalities:

- `fastq_files_rna`: Paired-end FASTQ files for scRNA sequencing.  
- `fastq_files_guide`: Paired-end FASTQ files for gRNA sequencing.  
- `fastq_files_hashing`: Paired-end FASTQ files for cell hashing.

### Test Files
Specific test files are configured for validation:
- Guide RNA test files (R1 and R2)
- Cell hashing test files (R1 and R2)

### Batch Configuration
- Batch identifiers: example 'batch_a', 'batch_b'
- Covariate list includes:
  - Batch information
  - Additional covariates (e.g., lane information)

## Usage Notes

1. All file paths should be relative to the project directory (`${projectDir}`)
2. FASTQ files should be properly paired and organized in the specified directory
3. The pipeline supports both local and remote genome references
4. Ensure all YAML specification files are present in the specified directory
5. Adjust threshold values and parameters according to your experimental design
6. Validate that all parameters have proper values - the pipeline will fail if any parameter is empty or missing
7. Use placeholder values like 'default', 'undefined', or 'undecided' when actual values are not yet determined, but never leave parameters empty

## Parameter Validation Requirements

- Every parameter in the configuration must be explicitly defined
- Use appropriate default values when actual values are not yet determined
- Ensure all file paths point to existing files or locations
- Boolean parameters must be explicitly set to 'true' or 'false'