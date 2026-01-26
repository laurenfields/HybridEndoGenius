"""
Peptide Fragment Analysis Script - Multi-Folder Batch Processing
=================================================================
This script processes peptide fragment data from MULTIPLE folders.

Each folder should contain:
1. final_results_EG_score.csv - contains target peptide IDs
2. fragment_matches.parquet - contains fragment match information

The script will:
- Find all folders in a parent directory
- Process each folder independently
- Combine all results into one master report
- Track which folder each result came from

Requirements:
- pandas
- pyarrow (for reading parquet files)

Instructions for Spyder:
1. Set your working directory to the PARENT folder containing all your data folders
2. Modify the parent_directory variable below if needed
3. Run the script (F5) or run sections with Ctrl+Enter
4. Check Variable Explorer to see all dataframes
"""

import pandas as pd
import os
import sys
import re
from pathlib import Path

#%% Configuration and Setup
print("="*70)
print("PEPTIDE FRAGMENT ANALYSIS - BATCH PROCESSING")
print("="*70)

# Parent directory containing all folders to process
# Options:
#   '.' = current working directory (default)
#   'path/to/parent/folder' = specific path
parent_directory = r"F:\postDefenseDataTransfer\Backup_D\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\crab_denovo_output_EGscore"


# File names to look for in each folder
csv_filename = 'final_results_EG_score.csv'
parquet_filename = 'fragment_matches.parquet'
output_file = r"F:\postDefenseDataTransfer\Backup_D\Manuscripts\2025_denovo_sequencing\by_analysis\peptide_ion_counts.csv"

# Filter thresholds
error_threshold = 0.02  # Fragment error (Da) must be <= this value
intensity_threshold = 1000  # Fragment actual intensity must be > this value

def get_peptide_length(peptide_seq):
    """Remove PTM annotations and return peptide length"""
    clean_seq = re.sub(r'\([^)]*\)', '', peptide_seq)
    return len(clean_seq)

def find_data_folders(parent_directory):
    """
    Find all data folders in the parent directory.
    Each sample folder should contain a subfolder with the data files.
    """
    parent_path = Path(parent_directory)
    if not parent_path.exists():
        print(f"   ERROR: Directory '{parent_directory}' not found!")
        return []
    
    folders_to_process = []
    
    for sample_folder in parent_path.iterdir():
        if sample_folder.is_dir():
            subfolders = [f for f in sample_folder.iterdir() if f.is_dir()]
            
            if len(subfolders) == 0:
                # Check if files are directly in this folder
                csv_path = sample_folder / csv_filename
                parquet_path = sample_folder / parquet_filename
                
                if csv_path.exists() and parquet_path.exists():
                    folders_to_process.append({
                        'sample_name': sample_folder.name,
                        'data_path': sample_folder
                    })
            else:
                # Take the first subfolder
                data_folder = subfolders[0]
                csv_path = data_folder / csv_filename
                parquet_path = data_folder / parquet_filename
                
                if csv_path.exists() and parquet_path.exists():
                    folders_to_process.append({
                        'sample_name': sample_folder.name,
                        'data_path': data_folder
                    })
    
    return folders_to_process

def process_analysis(analysis_name, parent_directory):
    """
    Process an entire analysis directory and return results dataframe.
    """
    print(f"\n{'='*70}")
    print(f"PROCESSING: {analysis_name}")
    print(f"{'='*70}")
    print(f"Directory: {parent_directory}")
    
    # Find folders
    print("\nFinding data folders...")
    folders = find_data_folders(parent_directory)
    
    if len(folders) == 0:
        print(f"   ERROR: No valid data folders found in {parent_directory}")
        return None
    
    print(f"   ✓ Found {len(folders)} sample folders")
    for folder_info in folders:
        print(f"      - {folder_info['sample_name']}")
    
    # Process each folder
    print("\nProcessing samples...")
    all_results = []
    
    for folder_info in folders:
        sample_name = folder_info['sample_name']
        data_path = folder_info['data_path']
        
        print(f"\n   Processing: {sample_name}")
        
        # Read CSV
        csv_path = data_path / csv_filename
        targets_df = pd.read_csv(csv_path)
        print(f"      ✓ Loaded {len(targets_df)} target peptides")
        
        # Read Parquet
        parquet_path = data_path / parquet_filename
        fragments_df = pd.read_parquet(parquet_path)
        print(f"      ✓ Loaded {len(fragments_df)} fragment matches")
        
        # Extract peptide-scan pairs
        target_pairs = targets_df[['Peptide', 'Scan']].drop_duplicates()
        
        # Filter fragments
        filtered_fragments = fragments_df[
            (fragments_df['Fragment error (Da)'] <= error_threshold) & 
            (fragments_df['Fragment actual intensity'] > intensity_threshold)
        ].copy()
        print(f"      ✓ Filtered to {len(filtered_fragments):,} fragments")
        
        # Process each peptide
        results = []
        for idx, row in target_pairs.iterrows():
            peptide = row['Peptide']
            scan = row['Scan']
            
            matches = filtered_fragments[
                (filtered_fragments['Sequence'] == peptide) & 
                (filtered_fragments['Scan'] == scan)
            ]
            
            if len(matches) == 0:
                b_ion_count = 0
                y_ion_count = 0
            else:
                ion_types = matches['ion'].value_counts()
                b_ion_count = sum(count for ion, count in ion_types.items() 
                                  if isinstance(ion, str) and ion.startswith('b'))
                y_ion_count = sum(count for ion, count in ion_types.items() 
                                  if isinstance(ion, str) and ion.startswith('y'))
            
            results.append({
                'Sample': sample_name,
                'Peptide': peptide,
                'Scan': scan,
                'Number of b-ions': b_ion_count,
                'Number of y-ions': y_ion_count
            })
        
        results_df = pd.DataFrame(results)
        all_results.append(results_df)
        print(f"      ✓ Processed {len(results_df)} peptides")
    
    # Combine all results
    print("\nCombining results...")
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Calculate coverage
    combined_df['Peptide_Length'] = combined_df['Peptide'].apply(get_peptide_length)
    combined_df['Expected_Ions_Per_Type'] = combined_df['Peptide_Length'] - 1
    combined_df['% Coverage b'] = (combined_df['Number of b-ions'] / combined_df['Expected_Ions_Per_Type']) * 100
    combined_df['% Coverage y'] = (combined_df['Number of y-ions'] / combined_df['Expected_Ions_Per_Type']) * 100
    combined_df['% Coverage b'] = combined_df['% Coverage b'].fillna(0)
    combined_df['% Coverage y'] = combined_df['% Coverage y'].fillna(0)
    combined_df['Total_Coverage'] = combined_df['% Coverage b'] + combined_df['% Coverage y']
    
    # Filter for high quality
    high_quality_df = combined_df[combined_df['Total_Coverage'] > coverage_threshold].copy()
    
    print(f"\n✓ Analysis complete:")
    print(f"   Total peptides: {len(combined_df)}")
    print(f"   High-quality peptides (>{coverage_threshold}% total coverage): {len(high_quality_df)} ({len(high_quality_df)/len(combined_df)*100:.1f}%)")
    
    return high_quality_df

#%% Process Both Analyses

# Process pre-de novo
pre_denovo_df = process_analysis("PRE-DE NOVO", pre_denovo_path)

# Process post-de novo
post_denovo_df = process_analysis("POST-DE NOVO", post_denovo_path)

# Check if both succeeded
if pre_denovo_df is None or post_denovo_df is None:
    print("\nERROR: One or both analyses failed to process!")
    sys.exit(1)

# Add analysis label
pre_denovo_df['Analysis'] = 'Pre-De Novo'
post_denovo_df['Analysis'] = 'Post-De Novo'

#%% Comparative Statistics

print("\n" + "="*70)
print("COMPARATIVE STATISTICS")
print("="*70)

print("\n1. SAMPLE SIZE:")
print(f"   Pre-De Novo:  {len(pre_denovo_df)} high-quality peptides")
print(f"   Post-De Novo: {len(post_denovo_df)} high-quality peptides")

print("\n2. B-ION COVERAGE:")
print(f"   Pre-De Novo:  Mean = {pre_denovo_df['% Coverage b'].mean():.2f}%, Median = {pre_denovo_df['% Coverage b'].median():.2f}%")
print(f"   Post-De Novo: Mean = {post_denovo_df['% Coverage b'].mean():.2f}%, Median = {post_denovo_df['% Coverage b'].median():.2f}%")
print(f"   Difference:   {post_denovo_df['% Coverage b'].mean() - pre_denovo_df['% Coverage b'].mean():.2f}%")

print("\n3. Y-ION COVERAGE:")
print(f"   Pre-De Novo:  Mean = {pre_denovo_df['% Coverage y'].mean():.2f}%, Median = {pre_denovo_df['% Coverage y'].median():.2f}%")
print(f"   Post-De Novo: Mean = {post_denovo_df['% Coverage y'].mean():.2f}%, Median = {post_denovo_df['% Coverage y'].median():.2f}%")
print(f"   Difference:   {post_denovo_df['% Coverage y'].mean() - pre_denovo_df['% Coverage y'].mean():.2f}%")

print("\n4. TOTAL COVERAGE:")
print(f"   Pre-De Novo:  Mean = {pre_denovo_df['Total_Coverage'].mean():.2f}%")
print(f"   Post-De Novo: Mean = {post_denovo_df['Total_Coverage'].mean():.2f}%")
print(f"   Difference:   {post_denovo_df['Total_Coverage'].mean() - pre_denovo_df['Total_Coverage'].mean():.2f}%")

print("\n5. B/Y RATIO:")
pre_b_y_ratio = (pre_denovo_df['% Coverage b'] / (pre_denovo_df['% Coverage y'] + 0.001)).median()
post_b_y_ratio = (post_denovo_df['% Coverage b'] / (post_denovo_df['% Coverage y'] + 0.001)).median()
print(f"   Pre-De Novo:  Median B/Y ratio = {pre_b_y_ratio:.2f}")
print(f"   Post-De Novo: Median B/Y ratio = {post_b_y_ratio:.2f}")

print("\n6. PEPTIDES FAVORING B-IONS:")
pre_favor_b = (pre_denovo_df['% Coverage b'] > pre_denovo_df['% Coverage y']).sum()
post_favor_b = (post_denovo_df['% Coverage b'] > post_denovo_df['% Coverage y']).sum()
print(f"   Pre-De Novo:  {pre_favor_b} / {len(pre_denovo_df)} ({pre_favor_b/len(pre_denovo_df)*100:.1f}%)")
print(f"   Post-De Novo: {post_favor_b} / {len(post_denovo_df)} ({post_favor_b/len(post_denovo_df)*100:.1f}%)")

# Statistical tests
print("\n7. STATISTICAL TESTS (Independent t-tests):")
t_stat_b, p_val_b = stats.ttest_ind(pre_denovo_df['% Coverage b'], post_denovo_df['% Coverage b'])
t_stat_y, p_val_y = stats.ttest_ind(pre_denovo_df['% Coverage y'], post_denovo_df['% Coverage y'])
t_stat_total, p_val_total = stats.ttest_ind(pre_denovo_df['Total_Coverage'], post_denovo_df['Total_Coverage'])

print(f"   B-ion coverage: t={t_stat_b:.3f}, p={p_val_b:.4e}")
print(f"   Y-ion coverage: t={t_stat_y:.3f}, p={p_val_y:.4e}")
print(f"   Total coverage: t={t_stat_total:.3f}, p={p_val_total:.4e}")

if p_val_b < 0.05:
    print(f"   ✓ B-ion coverage is SIGNIFICANTLY different between analyses (p < 0.05)")
if p_val_y < 0.05:
    print(f"   ✓ Y-ion coverage is SIGNIFICANTLY different between analyses (p < 0.05)")
if p_val_total < 0.05:
    print(f"   ✓ Total coverage is SIGNIFICANTLY different between analyses (p < 0.05)")


#%% Create Comparative Visualizations

print("\n" + "="*70)
print("CREATING COMPARATIVE PLOTS")
print("="*70)

# Calculate intensities for both analyses
print("\nCalculating ion intensities for both analyses...")

def calculate_ion_intensities(df, folders_list, parent_dir):
    """Calculate b-ion and y-ion intensities for a dataframe"""
    b_intensities = []
    y_intensities = []
    
    for folder_info in folders_list:
        sample_name = folder_info['sample_name']
        data_path = folder_info['data_path']
        
        parquet_path = data_path / parquet_filename
        fragments_df = pd.read_parquet(parquet_path)
        
        filtered_fragments = fragments_df[
            (fragments_df['Fragment error (Da)'] <= error_threshold) & 
            (fragments_df['Fragment actual intensity'] > intensity_threshold)
        ].copy()
        
        sample_results = df[df['Sample'] == sample_name]
        
        for idx, row in sample_results.iterrows():
            peptide = row['Peptide']
            scan = row['Scan']
            
            matches = filtered_fragments[
                (filtered_fragments['Sequence'] == peptide) & 
                (filtered_fragments['Scan'] == scan)
            ]
            
            if len(matches) > 0:
                b_ions = matches[matches['ion'].str.startswith('b', na=False)]
                y_ions = matches[matches['ion'].str.startswith('y', na=False)]
                
                b_avg_intensity = b_ions['Fragment actual intensity'].mean() if len(b_ions) > 0 else 0
                y_avg_intensity = y_ions['Fragment actual intensity'].mean() if len(y_ions) > 0 else 0
            else:
                b_avg_intensity = 0
                y_avg_intensity = 0
            
            b_intensities.append(b_avg_intensity)
            y_intensities.append(y_avg_intensity)
    
    return b_intensities, y_intensities

# Get folder lists for both analyses
pre_folders = find_data_folders(pre_denovo_path)
post_folders = find_data_folders(post_denovo_path)

# Calculate intensities
pre_b_int, pre_y_int = calculate_ion_intensities(pre_denovo_df, pre_folders, pre_denovo_path)
post_b_int, post_y_int = calculate_ion_intensities(post_denovo_df, post_folders, post_denovo_path)

pre_denovo_df['B_Avg_Intensity'] = pre_b_int
pre_denovo_df['Y_Avg_Intensity'] = pre_y_int
post_denovo_df['B_Avg_Intensity'] = post_b_int
post_denovo_df['Y_Avg_Intensity'] = post_y_int

print("   ✓ Calculated intensities for both analyses")

fig = plt.figure(figsize=(18, 12))

# ============================================================================
# PLOT 1: Combined scatter plot with simplified labels (2 regions only)
# ============================================================================
ax1 = fig.add_subplot(231)

# Pre-De Novo
ax1.scatter(pre_denovo_df['% Coverage y'], pre_denovo_df['% Coverage b'],
           color='orange', alpha=0.4, s=40, edgecolors='darkorange', linewidth=0.5,
           label=f'Pre-De Novo (n={len(pre_denovo_df)})')

# Post-De Novo
ax1.scatter(post_denovo_df['% Coverage y'], post_denovo_df['% Coverage b'],
           color='purple', alpha=0.4, s=40, edgecolors='indigo', linewidth=0.5,
           label=f'Post-De Novo (n={len(post_denovo_df)})')

max_axis = max(pre_denovo_df['% Coverage b'].max(), pre_denovo_df['% Coverage y'].max(),
               post_denovo_df['% Coverage b'].max(), post_denovo_df['% Coverage y'].max(), 100) + 10

ax1.plot([0, max_axis], [0, max_axis], 'k--', alpha=0.5, linewidth=2, label='Equal b and y')

# Calculate percentages for above/below diagonal only
# Above diagonal: more b-ions
pre_above = (pre_denovo_df['% Coverage b'] > pre_denovo_df['% Coverage y']).sum()
post_above = (post_denovo_df['% Coverage b'] > post_denovo_df['% Coverage y']).sum()

# Below diagonal: more y-ions
pre_below = (pre_denovo_df['% Coverage b'] < pre_denovo_df['% Coverage y']).sum()
post_below = (post_denovo_df['% Coverage b'] < post_denovo_df['% Coverage y']).sum()

# Add shaded regions
ax1.fill_between([0, max_axis], [0, max_axis], max_axis, alpha=0.1, color='red')
ax1.fill_between([0, max_axis], 0, [0, max_axis], alpha=0.1, color='blue')

# Add text labels - only 2 labels
# Above diagonal (more b-ions)
ax1.text(max_axis * 0.3, max_axis * 0.7, 
        f'More B-ions\n'
        f'Pre: {pre_above} ({pre_above/len(pre_denovo_df)*100:.1f}%)\n'
        f'Post: {post_above} ({post_above/len(post_denovo_df)*100:.1f}%)',
        fontsize=10, ha='center', va='center',
        bbox=dict(boxstyle='round', facecolor='salmon', alpha=0.7, edgecolor='darkred'))

# Below diagonal (more y-ions)
ax1.text(max_axis * 0.7, max_axis * 0.3,
        f'More Y-ions\n'
        f'Pre: {pre_below} ({pre_below/len(pre_denovo_df)*100:.1f}%)\n'
        f'Post: {post_below} ({post_below/len(post_denovo_df)*100:.1f}%)',
        fontsize=10, ha='center', va='center',
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7, edgecolor='darkblue'))

ax1.set_xlabel('% Coverage Y-ions', fontsize=11, fontweight='bold')
ax1.set_ylabel('% Coverage B-ions', fontsize=11, fontweight='bold')
ax1.set_title('B vs Y Coverage: Pre vs Post De Novo', fontsize=12, fontweight='bold')
ax1.set_xlim(-5, max_axis)
ax1.set_ylim(-5, max_axis)
ax1.grid(True, alpha=0.3)
ax1.set_aspect('equal')
ax1.legend(loc='upper left', fontsize=9)

# ============================================================================
# PLOT 2: Paired box plot - B-ions and Y-ions together
# ============================================================================
ax2 = fig.add_subplot(232)

# Prepare data: [Pre B, Pre Y, Post B, Post Y]
box_data_paired = [
    pre_denovo_df['% Coverage b'],
    pre_denovo_df['% Coverage y'],
    post_denovo_df['% Coverage b'],
    post_denovo_df['% Coverage y']
]

positions = [1, 2, 4, 5]
bp_paired = ax2.boxplot(box_data_paired, positions=positions, 
                        labels=['B-ions', 'Y-ions', 'B-ions', 'Y-ions'],
                        patch_artist=True, showmeans=True, meanline=True, widths=0.6)

# Color the boxes
bp_paired['boxes'][0].set_facecolor('orange')  # Pre B
bp_paired['boxes'][0].set_alpha(0.7)
bp_paired['boxes'][1].set_facecolor('orange')  # Pre Y
bp_paired['boxes'][1].set_alpha(0.7)
bp_paired['boxes'][2].set_facecolor('purple')  # Post B
bp_paired['boxes'][2].set_alpha(0.7)
bp_paired['boxes'][3].set_facecolor('purple')  # Post Y
bp_paired['boxes'][3].set_alpha(0.7)

# Add group labels
y_lim = ax2.get_ylim()
y_offset = y_lim[0] - (y_lim[1] - y_lim[0]) * 0.08
ax2.text(1.5, y_offset, 'Pre-De Novo', ha='center', fontsize=11, fontweight='bold')
ax2.text(4.5, y_offset, 'Post-De Novo', ha='center', fontsize=11, fontweight='bold')

# Add mean annotations
mean_pre_b = pre_denovo_df['% Coverage b'].mean()
mean_pre_y = pre_denovo_df['% Coverage y'].mean()
mean_post_b = post_denovo_df['% Coverage b'].mean()
mean_post_y = post_denovo_df['% Coverage y'].mean()

ax2.text(1, mean_pre_b, f'{mean_pre_b:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=8)
ax2.text(2, mean_pre_y, f'{mean_pre_y:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=8)
ax2.text(4, mean_post_b, f'{mean_post_b:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=8)
ax2.text(5, mean_post_y, f'{mean_post_y:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=8)

ax2.set_ylabel('% Coverage', fontsize=11, fontweight='bold')
ax2.set_title('Ion Coverage Comparison (Paired)', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_xlim(0, 6)

# ============================================================================
# PLOT 3: Total Coverage Comparison
# ============================================================================
ax3 = fig.add_subplot(233)

box_data_total = [pre_denovo_df['Total_Coverage'], post_denovo_df['Total_Coverage']]
bp_total = ax3.boxplot(box_data_total, labels=['Pre-De Novo', 'Post-De Novo'], 
                       patch_artist=True, showmeans=True, meanline=True)
bp_total['boxes'][0].set_facecolor('orange')
bp_total['boxes'][0].set_alpha(0.7)
bp_total['boxes'][1].set_facecolor('purple')
bp_total['boxes'][1].set_alpha(0.7)

ax3.set_ylabel('% Total Coverage (B + Y)', fontsize=11, fontweight='bold')
ax3.set_title('Total Coverage Comparison', fontsize=12, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='y')

mean_pre_total = pre_denovo_df['Total_Coverage'].mean()
mean_post_total = post_denovo_df['Total_Coverage'].mean()
ax3.text(1, mean_pre_total, f'{mean_pre_total:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=9)
ax3.text(2, mean_post_total, f'{mean_post_total:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=9)

# ============================================================================
# PLOT 4: Ion Count Comparison (Paired)
# ============================================================================
ax4 = fig.add_subplot(234)

# Ion count data
box_data_counts = [
    pre_denovo_df['Number of b-ions'],
    pre_denovo_df['Number of y-ions'],
    post_denovo_df['Number of b-ions'],
    post_denovo_df['Number of y-ions']
]

positions_counts = [1, 2, 4, 5]
bp_counts = ax4.boxplot(box_data_counts, positions=positions_counts,
                        labels=['B-ions', 'Y-ions', 'B-ions', 'Y-ions'],
                        patch_artist=True, showmeans=True, meanline=True, widths=0.6)

# Color the boxes
bp_counts['boxes'][0].set_facecolor('orange')
bp_counts['boxes'][0].set_alpha(0.7)
bp_counts['boxes'][1].set_facecolor('orange')
bp_counts['boxes'][1].set_alpha(0.7)
bp_counts['boxes'][2].set_facecolor('purple')
bp_counts['boxes'][2].set_alpha(0.7)
bp_counts['boxes'][3].set_facecolor('purple')
bp_counts['boxes'][3].set_alpha(0.7)

# Add group labels
y_lim4 = ax4.get_ylim()
y_offset4 = y_lim4[0] - (y_lim4[1] - y_lim4[0]) * 0.08
ax4.text(1.5, y_offset4, 'Pre-De Novo', ha='center', fontsize=11, fontweight='bold')
ax4.text(4.5, y_offset4, 'Post-De Novo', ha='center', fontsize=11, fontweight='bold')

ax4.set_ylabel('Number of Ions Detected', fontsize=11, fontweight='bold')
ax4.set_title('Ion Count Comparison (Paired)', fontsize=12, fontweight='bold')
ax4.grid(True, alpha=0.3, axis='y')
ax4.set_xlim(0, 6)

# ============================================================================
# PLOT 5: Average Intensity Comparison (LOG SCALE)
# ============================================================================
ax5 = fig.add_subplot(235)

# Calculate average intensities (b + y) / 2
pre_denovo_df['Avg_Total_Intensity'] = (pre_denovo_df['B_Avg_Intensity'] + pre_denovo_df['Y_Avg_Intensity']) / 2
post_denovo_df['Avg_Total_Intensity'] = (post_denovo_df['B_Avg_Intensity'] + post_denovo_df['Y_Avg_Intensity']) / 2

# Filter non-zero values for log scale
pre_nonzero_int = pre_denovo_df['Avg_Total_Intensity'][pre_denovo_df['Avg_Total_Intensity'] > 0]
post_nonzero_int = post_denovo_df['Avg_Total_Intensity'][post_denovo_df['Avg_Total_Intensity'] > 0]

box_data_intensity = [pre_nonzero_int, post_nonzero_int]

bp_int = ax5.boxplot(box_data_intensity, labels=['Pre-De Novo', 'Post-De Novo'], 
                     patch_artist=True, showmeans=True, meanline=True)
bp_int['boxes'][0].set_facecolor('orange')
bp_int['boxes'][0].set_alpha(0.7)
bp_int['boxes'][1].set_facecolor('purple')
bp_int['boxes'][1].set_alpha(0.7)

ax5.set_ylabel('Average Fragment Intensity (log scale)', fontsize=11, fontweight='bold')
ax5.set_title('Average Intensity Comparison', fontsize=12, fontweight='bold')
ax5.set_yscale('log')  # Set log scale
ax5.grid(True, alpha=0.3, axis='y', which='both')  # Show grid for both major and minor ticks

# ============================================================================
# PLOT 6: B/Y Ratio comparison
# ============================================================================
ax6 = fig.add_subplot(236)

# Calculate ratios
pre_ratios = pre_denovo_df['% Coverage b'] / (pre_denovo_df['% Coverage y'] + 0.001)
post_ratios = post_denovo_df['% Coverage b'] / (post_denovo_df['% Coverage y'] + 0.001)

# Filter extreme values
pre_ratios_filt = pre_ratios[pre_ratios < 5]
post_ratios_filt = post_ratios[post_ratios < 5]

ax6.hist(pre_ratios_filt, bins=30, alpha=0.6, color='orange', label='Pre-De Novo', edgecolor='darkorange')
ax6.hist(post_ratios_filt, bins=30, alpha=0.6, color='purple', label='Post-De Novo', edgecolor='indigo')
ax6.axvline(x=1, color='red', linestyle='--', linewidth=2, label='Equal (ratio=1)')

ax6.set_xlabel('B-ion / Y-ion Coverage Ratio', fontsize=11, fontweight='bold')
ax6.set_ylabel('Number of Peptides', fontsize=11, fontweight='bold')
ax6.set_title('B/Y Ratio Distribution Comparison', fontsize=12, fontweight='bold')
ax6.legend(fontsize=9)
ax6.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()

print(f"\n   ✓ Comparison plots created!")

# Save figure
fig.savefig(output_plots, dpi=300, bbox_inches='tight')
print(f"   ✓ Plots saved to: {output_plots}")

#%% Save Comparison Data

# Combine for export
combined_export = pd.concat([pre_denovo_df, post_denovo_df], ignore_index=True)
combined_export.to_csv(output_comparison_csv, index=False)
print(f"\n   ✓ Comparison data saved to: {output_comparison_csv}")

print("\n" + "="*70)
print("COMPARISON ANALYSIS COMPLETE!")
print("="*70)
print(f"\nFiles saved:")
print(f"  - {output_comparison_csv}")
print(f"  - {output_plots}")
print(f"\nKey variables in workspace:")
print(f"  - pre_denovo_df: Pre-de novo results ({len(pre_denovo_df)} peptides)")
print(f"  - post_denovo_df: Post-de novo results ({len(post_denovo_df)} peptides)")
print(f"  - combined_export: Combined data for export")
print("="*70)