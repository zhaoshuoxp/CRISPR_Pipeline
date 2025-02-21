#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData
import matplotlib.pyplot as plt
import os
from typing import Literal, Optional

def volcano(mdata: MuData, method: Optional[Literal['sceptre', 'perturbo']], 
          log2_fc_thresh: float, p_value_thresh: float, ax: plt.Axes = None):
   if ax is None:
       ax = plt.gca()
   
   # Get test results from MuData
   test_results = pd.DataFrame({k: v for k, v in mdata.uns['test_results'].items()})
   
   # Determine which columns to use
   if not method:  # If method is None, try to use generic columns
       if 'log2_fc' in test_results.columns and 'p_value' in test_results.columns:
           log2_fc_col = 'log2_fc'
           p_value_col = 'p_value'
       else:
           raise ValueError("No generic log2_fc and p_value columns found")
   else:  # Use method-specific columns
       log2_fc_col = f"{method}_log2_fc"
       p_value_col = f"{method}_p_value"
   
   # Plot all data points (non-significant genes)
   ax.scatter(x=test_results[log2_fc_col], 
             y=test_results[p_value_col].apply(lambda x: -np.log10(x)), 
             s=5, color="green", label="Not significant")
   
   # Highlight down or up-regulated genes based on the thresholds
   down = test_results[(test_results[log2_fc_col] <= -log2_fc_thresh) & 
                      (test_results[p_value_col] <= p_value_thresh)]
   up = test_results[(test_results[log2_fc_col] >= log2_fc_thresh) & 
                     (test_results[p_value_col] <= p_value_thresh)]

   filtered_down = down[np.isfinite(down[log2_fc_col]) & 
                       np.isfinite(down[p_value_col])]

   method_name = method if method else "Analysis"
   print(f"\n{method_name.capitalize()} results:")
   print(f"Number of down-regulated genes based on the thresholds: {len(down)}")
   print(f"Number of down-regulated genes (finite values): {len(filtered_down)}")
   print(f"Number of up-regulated genes based on the thresholds: {len(up)}")

   down_copy = down.copy()
   down_copy[log2_fc_col] = down_copy[log2_fc_col].replace([np.inf, -np.inf], [1e10, -1e10])
   down_copy[p_value_col] = down_copy[p_value_col].replace([np.inf, -np.inf], [1e10, -1e10])

   top_5_down_log2fc = down_copy.sort_values(by=log2_fc_col, ascending=True).head(10)
   top_5_down_pval = down_copy.sort_values(by=p_value_col, ascending=True).head(10)
   print(f"\nTop genes ({method_name}):")
   print("Top down-regulated by log2fc:", top_5_down_log2fc['intended_target_name'].tolist())
   print("Top down-regulated by p-value:", top_5_down_pval['intended_target_name'].tolist())
   annotated = pd.concat([top_5_down_log2fc, top_5_down_pval])

   # Plot down-regulated genes
   ax.scatter(x=down[log2_fc_col], 
             y=down[p_value_col].apply(lambda x: -np.log10(x)), 
             s=10, label="Down-regulated", color="blue")
   
   # Plot up-regulated genes
   ax.scatter(x=up[log2_fc_col], 
             y=up[p_value_col].apply(lambda x: -np.log10(x)), 
             s=10, label="Up-regulated", color="red")
   
   # Annotate genes with offset labels
   for index, row in annotated.iterrows():
       x = row[log2_fc_col]
       y = -np.log10(row[p_value_col])
       label = f"{row['intended_target_name']} ({row['guide_id']})"
       
       # Add small offset to avoid overlapping with points
       offset_x = 0.1 if x >= 0 else -0.1
       offset_y = 0.1
       
       # Add annotation with a simple arrow
       ax.annotate(label, 
                  xy=(x, y),
                  xytext=(x + offset_x, y + offset_y),
                  fontsize=8,
                  arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

   # Set plot boundaries and labels
   low, high = ax.get_xlim()
   bound = max(abs(low), abs(high)) + 0.5
   
   ax.set_xlabel("log2 fold change")
   ax.set_ylabel("p-value (-log10)")
   
   # Draw threshold lines
   ax.axvline(-log2_fc_thresh, color="grey", linestyle="--")
   ax.axvline(log2_fc_thresh, color="grey", linestyle="--")
   ax.axhline(-np.log10(p_value_thresh), color="grey", linestyle="--")
   
   ax.set_xlim(-bound, bound)
   ax.legend(loc="upper right")
   
   method_title = method.capitalize() if method else "Differential Expression"
   ax.set_title(f"{method_title} Volcano Plot")

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Generate comparison volcano plots from MuData")
   parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
   parser.add_argument("--log2_fc", type=float, default=1,
                     help="log2 fold change threshold")
   parser.add_argument("--p_value", type=float, default=0.01,
                     help="p-value threshold")
   
   args = parser.parse_args()
   
   mdata = mu.read(args.mdata_path)
   
   # Check available methods
   results_df = pd.DataFrame(mdata.uns['test_results'])
   cols = results_df.columns
   
   # Check for generic columns first
   if 'log2_fc' in cols and 'p_value' in cols:
       print("Using generic log2_fc and p_value columns")
       fig, ax = plt.subplots(1, 1, figsize=(10, 10))
       volcano(mdata, None, args.log2_fc, args.p_value, ax)
       plt.suptitle("Differential Expression Results", y=1.02)
   else:
       # Check for method-specific columns
       available_methods = []
       if 'sceptre_log2_fc' in cols and not results_df['sceptre_log2_fc'].isna().all():
           available_methods.append('sceptre')
       if 'perturbo_log2_fc' in cols and not results_df['perturbo_log2_fc'].isna().all():
           available_methods.append('perturbo')
           
       print(f"Available methods: {available_methods}")
       
       if len(available_methods) == 1:
           fig, ax = plt.subplots(1, 1, figsize=(10, 10))
           volcano(mdata, available_methods[0], args.log2_fc, args.p_value, ax)
           plt.suptitle(f"{available_methods[0].title()} Results", y=1.02)
       elif len(available_methods) == 2:
           fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
           volcano(mdata, 'sceptre', args.log2_fc, args.p_value, ax1)
           volcano(mdata, 'perturbo', args.log2_fc, args.p_value, ax2)
           plt.suptitle("Sceptre vs Perturbo Results", y=1.02)
   
   # Save the plot
   output_dir = "evaluation_output"
   os.makedirs(output_dir, exist_ok=True)
   output_file = os.path.join(output_dir, "volcano_plot.png")

   print(f"Saving plot to {output_file}")
   plt.tight_layout()
   plt.savefig(output_file, dpi=300, bbox_inches='tight')
   plt.close()
   print("Plot saved successfully")