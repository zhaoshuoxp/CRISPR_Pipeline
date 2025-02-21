#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np

def export_output(mudata, inference_method):
    """
    Generate per-guide and per-element outputs from MuData and save as TSV files.
    """
    # Check if test_results exists in uns
    if 'test_results' not in mudata.uns:
        print("Warning: No test results found in mudata.uns")
        return

    # Get test results 
    test_results = pd.DataFrame(mudata.uns['test_results'])
    if test_results.empty:
        print("Warning: Test results dataframe is empty")
        return

    # Define column mappings
    col_map = {
        'sceptre': {'log2_fc': 'sceptre_log2_fc', 'p_value': 'sceptre_p_value'},
        'perturbo': {'log2_fc': 'perturbo_log2_fc', 'p_value': 'perturbo_p_value'}
    }
    
    if inference_method not in col_map:
        raise ValueError("Invalid inference_method. Must be 'sceptre' or 'perturbo'.")

    try:
        # merge with guide data
        per_guide_output = test_results.merge(
            mudata.mod['guide'].var.reset_index(),
            how='left',
            on=['intended_target_name', 'guide_id']
        )

        # Calculate average gene expression 
        try:  
            avg_expressions = []
            for _, row in per_guide_output.iterrows():  # Use per_guide_output instead of test_results
                target = row['intended_target_name']
                mask = mudata.mod['gene'].var['symbol'] == target
                if mask.any():
                    gene_data = mudata.mod['gene'].X[:, mask]
                    if isinstance(gene_data, np.ndarray):
                        avg_exp = np.mean(gene_data)
                    else:  # Sparse matrix
                        avg_exp = np.mean(gene_data.toarray())
                else:
                    avg_exp = np.nan
                avg_expressions.append(avg_exp)
            
            print(f"Length of avg_expressions: {len(avg_expressions)}")
            per_guide_output['avg_gene_expression'] = avg_expressions

        except Exception as e:
            print(f"Error calculating average expression: {str(e)}")
            print(f"Shape of gene expression matrix: {mudata.mod['gene'].X.shape}")
            raise

        # Add columns
        per_guide_output['cell_number'] = mudata.shape[0]
        per_guide_output['avg_gene_expression'] = avg_expressions

        # Rename columns based on inference method
        per_guide_output = per_guide_output.rename(columns=col_map[inference_method])

        # Add missing inference method columns
        for col in ['sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc', 'perturbo_p_value']:
            if col not in per_guide_output.columns:
                per_guide_output[col] = None

        # Reorder and rename columns
        final_cols = [
            'intended_target_name', 'guide_id', 'intended_target_chr',
            'intended_target_start', 'intended_target_end', 'gene_id',
            'sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc',
            'perturbo_p_value', 'cell_number', 'avg_gene_expression'
        ]

        per_guide_output = per_guide_output[final_cols].rename(columns={'guide_id': 'guide_id(s)'})

        # Generate per-element output
        per_element_output = (
            per_guide_output.groupby('gene_id', as_index=False)
            .agg({'guide_id(s)': lambda x: ','.join(x.dropna().astype(str))})
            .merge(
                per_guide_output.drop_duplicates('gene_id')
                .drop(columns=['guide_id(s)']),
                on='gene_id',
                how='left'
            )
        )

        column_order = [
            'intended_target_name',
            'guide_id(s)',
            'intended_target_chr',
            'intended_target_start',
            'intended_target_end',
            'gene_id',
            'sceptre_log2_fc',
            'sceptre_p_value',
            'perturbo_log2_fc',
            'perturbo_p_value',
            'cell_number',
            'avg_gene_expression'
        ]

        per_element_output = per_element_output[column_order]
        # Save outputs with error handling
        per_guide_output.to_csv('per_guide_output.tsv', sep='\t', index=False)
        per_element_output.to_csv('per_element_output.tsv', sep='\t', index=False)
        print("Exported output files successfully.")
        
    except Exception as e:
        print(f"Error during export: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export outputs from MuData.")
    parser.add_argument("--inference_method", type=str, required=True, help="Inference method: 'sceptre' or 'perturbo'.")
    parser.add_argument("--mudata", type=str, required=True, help="Path to input MuData file.")
    args = parser.parse_args()

    try:
        mudata = mu.read_h5mu(args.mudata)
        export_output(mudata, args.inference_method)
    except Exception as e:
        print(f"Error loading or processing MuData file: {str(e)}")
        raise