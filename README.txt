The quality control (QC) of NGS data can be conducted by following the steps outlined in the provided R scripts. Below is a detailed description of each script and its purpose:

File1_ApplyinheritanceQC.r:

This script processes VCF files and calculates the Identity-by-State (IBS) between sample pairs.
It is recommended to split VCF files by chromosome to facilitate easier handling of large files and improve processing speed by working with SNPs in smaller chunks.
Adjust the 'num.cores' parameter in the script to define the number of cores to use for parallel processing.
The output includes files with counts of unfiltered and filtered SNPs in their names. IBS calculation begins with 1,000 unfiltered SNPs, doubling the number used iteratively until all SNPs are processed.
Use the final output file with the maximum number of unfiltered SNPs for your analysis. This script generates output files for both 'RELATED' pairs (IBS ≥ 0.06) and all sample pairs. For this method, only the file with 'RELATED' pairs is needed.
This script requires subroutines (hwe.r and annotate_SNPs_subroutines.r) to be correctly read from their paths. Output files are saved in the directory specified by 'path.save'.
File2_Population_stratification.r:

This script performs Principal Component Analysis (PCA) and includes detailed instructions for correcting population stratification.
Pay close attention to the comments in the script, which provide guidance on manipulating data files using tools such as PLINK, PLINK2, Shellfish, GATK, and Picard.
File3_Plot_PCA_and_infer_ethnicity_4d.r:

This script is used to plot PCA and infer the ethnicity of individuals.
Load the final 'RELATED' output file from File1_ApplyinheritanceQC.r into this script.
Using the first four principal components, it calculates the Euclidean distance between case samples and HapMap samples to infer sample ethnicity. The inferred ethnicities are saved in the output file.
Note: Execute all three scripts (File1_ApplyinheritanceQC.r, File2_Population_stratification.r, and File3_Plot_PCA_and_infer_ethnicity_4d.r) line by line to ensure no steps are missed during data processing, reading, or handling.
