#!/usr/bin/env python3
"""
rMATS Results Aggregator

This script aggregates rMATS output files from multiple comparisons and 
alternative splicing types into a single Excel file with separate sheets 
for JC and JCEC results.

Usage: python aggregate_rmats_results.py
"""

import pandas as pd
import os
import glob
from pathlib import Path

def main():
    # Define paths and parameters
    base_output_dir = "/global/scratch/users/enricocalvane/minRNAseq/as/rmats_output"
    output_file = os.path.join(base_output_dir, "aggregated_rmats_results.xlsx")
    
    # Define comparison directories and their labels
    comparisons = {
        "wt_vs_siz1": "WTvssiz1",
        "wt_vs_siz1_oe1": "WTvsSIZ1_OE1", 
        "wt_vs_siz1_oe2": "WTvsSIZ1_OE2"
    }
    
    # Define alternative splicing types
    as_types = ["SE", "A5SS", "A3SS", "MXE", "RI"]
    
    # Define result types (JC and JCEC)
    result_types = ["JC", "JCEC"]
    
    print("Starting rMATS results aggregation...")
    print(f"Base directory: {base_output_dir}")
    print(f"Comparisons: {list(comparisons.keys())}")
    print(f"AS types: {as_types}")
    print("="*50)
    
    # Initialize dictionaries to store aggregated data
    aggregated_data = {"JC": [], "JCEC": []}
    
    # Process each comparison
    for comp_dir, comp_label in comparisons.items():
        comp_path = os.path.join(base_output_dir, comp_dir)
        
        if not os.path.exists(comp_path):
            print(f"Warning: Directory not found: {comp_path}")
            continue
            
        print(f"Processing comparison: {comp_label}")
        print(f"  Directory: {comp_path}")
        
        # Process each alternative splicing type
        for as_type in as_types:
            print(f"  Processing {as_type}...")
            
            # Process both JC and JCEC files
            for result_type in result_types:
                file_pattern = f"{as_type}.MATS.{result_type}.txt"
                file_path = os.path.join(comp_path, file_pattern)
                
                if os.path.exists(file_path):
                    try:
                        # Read the file
                        df = pd.read_csv(file_path, sep='\t')
                        
                        if len(df) > 0:
                            # Add metadata columns
                            df['AS_Type'] = as_type
                            df['Comparison'] = comp_label
                            df['Result_Type'] = result_type
                            
                            # Add to aggregated data
                            aggregated_data[result_type].append(df)
                            
                            print(f"    ✓ {result_type}: {len(df)} events")
                        else:
                            print(f"    - {result_type}: Empty file")
                            
                    except Exception as e:
                        print(f"    ✗ {result_type}: Error reading file - {e}")
                else:
                    print(f"    - {result_type}: File not found")
    
    print("="*50)
    print("Combining data and creating Excel file...")
    
    # Create Excel writer
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        
        for result_type in result_types:
            if aggregated_data[result_type]:
                # Combine all dataframes for this result type
                combined_df = pd.concat(aggregated_data[result_type], 
                                      ignore_index=True, sort=False)
                
                # Reorder columns to put metadata first
                metadata_cols = ['Comparison', 'AS_Type', 'Result_Type']
                other_cols = [col for col in combined_df.columns if col not in metadata_cols]
                final_cols = metadata_cols + other_cols
                combined_df = combined_df[final_cols]
                
                # Sort by comparison and AS type for better organization
                combined_df = combined_df.sort_values(['Comparison', 'AS_Type'])
                
                # Write to Excel sheet
                sheet_name = result_type
                combined_df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                print(f"Sheet '{sheet_name}': {len(combined_df)} total events")
                
                # Print summary by comparison and AS type
                summary = combined_df.groupby(['Comparison', 'AS_Type']).size()
                print("  Event counts by comparison and AS type:")
                for (comp, as_type), count in summary.items():
                    print(f"    {comp} - {as_type}: {count}")
                    
            else:
                print(f"No data found for {result_type}")
    
    print("="*50)
    print(f"✓ Aggregation complete!")
    print(f"Output file: {output_file}")
    
    # Final summary
    if os.path.exists(output_file):
        print(f"File size: {os.path.getsize(output_file) / 1024 / 1024:.2f} MB")
        
        # Quick verification of the Excel file
        try:
            with pd.ExcelFile(output_file) as xls:
                print(f"Excel sheets created: {xls.sheet_names}")
                for sheet in xls.sheet_names:
                    df_check = pd.read_excel(xls, sheet_name=sheet, nrows=0)
                    print(f"  {sheet} columns: {list(df_check.columns[:10])}{'...' if len(df_check.columns) > 10 else ''}")
        except Exception as e:
            print(f"Warning: Could not verify Excel file - {e}")
    else:
        print("✗ Output file was not created")

if __name__ == "__main__":
    main()
