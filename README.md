## This is a modified version of the [IGVF CRISPR pipeline](https://github.com/pinellolab/CRISPR_Pipeline), specifically adapted for TAPseq and used in the CRISPR Jamboree 2025. After cloning the repository, you will need to decompress the barcode whitelist using the following command:  

```bash
gunzip 9p21_data/yaml_files/3M-february-2018.txt.gz
```




# CRISPR Pipeline

A comprehensive pipeline for single-cell Perturb-Seq analysis that enables robust processing and analysis of CRISPR screening data at single-cell resolution.

## Prerequisites

The following dependencies must be installed before running the pipeline:

### Nextflow
Workflow manager for executing the pipeline:
```bash
conda install bioconda::nextflow
```

### Mamba
Fast package manager for dependency resolution:
```bash
conda install conda-forge::mamba
```

### Singularity
Container platform that must be available in your execution environment.

## Installation

To install the pipeline:

```bash
git clone https://github.com/pinellolab/CRISPR_Pipeline.git
```

## Input Requirements

### File Descriptions

#### FASTQ Files
- `{sample}_R1.fastq.gz`: Contains cell barcode and UMI sequences
- `{sample}_R2.fastq.gz`: Contains transcript sequences

#### YAML Configuration Files
- `rna_seqspec.yml`: Defines RNA sequencing structure and parameters
- `guide_seqspec.yml`: Specifies guide RNA detection parameters
- `hash_seqspec.yml`: Defines cell hashing structure (required if using cell hashing)
- `whitelist.txt`: List of valid cell barcodes

#### Metadata Files
- `guide_metadata.tsv`: Contains guide RNA information and annotations
- `pairs_to_test.csv`: Defines perturbation pairs for comparison analysis (required if testing predefined pairs)

This pipeline requires a specific data structure to function properly. Below is an overview of the required directory organization:

```
📁 example_data/
   │
   ├── 📁 fastq_files/
   │   ├── 📄 {sample}_R1.fastq.gz
   │   └── 📄 {sample}_R2.fastq.gz
   │
   ├── 📁 yaml_files/
   │   ├── 📄 rna_seqspec.yml
   │   ├── 📄 guide_seqspec.yml
   │   └── 📄 whitelist.txt
   │
   ├── 📄 guide_metadata.tsv
   └── 📄 pairs_to_test.csv
```

For detailed specifications, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

## Running the Pipeline 

### Configuration Steps

#### 1. Pipeline Settings

Make sure to specify your data paths and analysis parameters in `configs/pipeline.config`.

#### 2. Resource Configuration

Configure `input.config` to match your computing environment. For example:

```bash
withName:process_name {
   cpus = 24               # Number of CPU cores per mapping process (default: 4)
   memory = 64.GB         # RAM allocation per mapping process (default: 64GB)
}
```
💡 **Note:** Start with these default values and adjust based on your dataset size and system capabilities.

### Running the Pipeline

1. First, make the scripts executable:
   ```bash
   chmod +x bin/*
   ```

2. Launch the pipeline:
   ```bash
   nextflow run main.nf -c input.config
   ```

### Monitoring and Troubleshooting

#### During Execution
- Watch the terminal output for progress updates
- Check the `.nextflow.log` file for detailed execution logs

#### Common Issues and Solutions
- **Memory errors**: Increase the `memory` parameter in `input.config`
- **Missing files**: Double-check paths in `configs/pipeline.config` and actual files in `example_data`

Need help? Check our documentation or raise an issue on GitHub.

## Output Description

The output files will be generated in the `pipeline_outputs` and `pipeline_dashboard` directory.

### Generated Files

Within the `pipeline_outputs` directory, you will find:

- inference_mudata.h5mu - MuData format output
- per_element_output.tsv - Per-element analysis
- per_guide_output.tsv - Per-guide analysis

**Structure:**

```
📁 pipeline_outputs/
   ├── 📄 inference_mudata.h5mu    
   ├── 📄 per_element_output.tsv    
   └── 📄 per_guide_output.tsv     
```

For details, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

### Generated Figures

The pipeline produces several figures:

Within the `pipeline_dashboard` directory, you will find:

1. **Evaluation Output**:
   - `network_plot.png`: Gene interaction networks visualization.
   - `volcano_plot.png`: gRNA-gene pairs analysis.
   - IGV files (`.bedgraph` and `bedpe`): Genome browser visualization files.

2. **Analysis Figures**:
   - `knee_plot_scRNA.png`: Knee plot of UMI counts vs. barcode index.
   - `scatterplot_scrna.png`: Scatterplot of total counts vs. genes detected, colored by mitochondrial content.
   - `violin_plot.png`: Distribution of gene counts, total counts, and mitochondrial content.
   - `scRNA_barcodes_UMI_thresholds.png`: Number of scRNA barcodes using different Total UMI thresholds.
   - `guides_per_cell_histogram.png`: Histogram of guides per cell.
   - `cells_per_guide_histogram.png`: Histogram of cells per guide.
   - `guides_UMI_thresholds.png`: Simulating the final number of cells with assigned guides using different minimal number thresholds (at least one guide > threshold value). (Use it to inspect how many cells would have assigned guides. This can be used to check if the final number of cells with guides fit with your expected number of cells)
   - `guides_UMI_thresholds.png`: Histogram of the number of sgRNA represented per cell
   - `cells_per_htp_barplot.png`: Number of Cells across Different HTOs
   - `umap_hto.png`: UMAP Clustering of Cells Based on HTOs (The dimensions represent the distribution of HTOs in each cell)
   - `umap_hto_singlets.png`: UMAP Clustering of Cells Based on HTOs (multiplets removed)

3. **seqSpec Plots**:

   - `seqSpec_check_plots.png`: The frequency of each nucleotides along the Read 1 (Use to inspect the expected read parts with their expected signature) and Read 2 (Use to inspect the expected read parts with their expected signature).

**Structure:**
```
📁 pipeline_dashboard/
  ├── 📄 dashboard.html                         
  │
  ├── 📁 evaluation_output/                      
  │   ├── 🖼️ network_plot.png                   
  │   ├── 🖼️ volcano_plot.png                  
  │   ├── 📄 igv.bedgraph                     
  │   └── 📄 igv.bedpe                         
  │
  ├── 📁 figures/
  │   ├── 🖼️ knee_plot_scRNA.png                
  │   ├── 🖼️ scatterplot_scrna.png              
  │   ├── 🖼️ violin_plot.png                    
  │   ├── 🖼️ scRNA_barcodes_UMI_thresholds.png  
  │   ├── 🖼️ guides_per_cell_histogram.png      
  │   ├── 🖼️ cells_per_guide_histogram.png      
  │   ├── 🖼️ guides_UMI_thresholds.png          
  │   ├── 🖼️ cells_per_htp_barplot.png          
  │   ├── 🖼️ umap_hto.png                       
  │   └── 🖼️ umap_hto_singlets.png              
  │
  ├── 📁 guide_seqSpec_plots/
  │   └── 🖼️ seqSpec_check_plots.png            
  │
  └── 📁 hashing_seqSpec_plots/
      └── 🖼️ seqSpec_check_plots.png             
```

## Pipeline Testing Guide

To ensure proper pipeline functionality, we provide two extensively validated datasets for testing purposes.

### Available Test Datasets

#### 1. SMC DC TAPseq tiny Dataset (Quertermous Lab)

### Step-by-Step Testing Instructions

1. **Environment Setup**
   ```bash
   # Clone and enter the repository
   git clone https://github.com/zhaoshuoxp/CRISPR_Pipeline.git
   cd CRISPR_Pipeline
   gunzip 9p21_data/yaml_files/3M-february-2018.txt.gz
   ```
   
2. **Follow the Instructions:**

   ```bash
   # Run with default configuration
   nextflow run main.nf -c input.config
   ```
   
   3. Expected Outputs
      

The pipeline generates two directories upon completion:
- `pipeline_outputs`: Contains all analysis results
- `pipeline_dashboard`: Houses interactive visualization reports

### Troubleshooting
If you encounter any issues during testing:
1. Review log files and intermediate results in the `work/` directory
2. Verify that all input files meet the required format specifications

For additional support or questions, please open an issue on our GitHub repository.
