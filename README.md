# Random Forest Analysis for Gene Expression Data

This project applies Random Forest analysis to gene expression data (e.g., RNA-Seq counts) to identify the most important genes related to a specified outcome variable. Additionally, it performs KEGG pathway enrichment analysis on the most important genes.

# Key Features:

Load and preprocess RNA-Seq counts data.

Perform Random Forest analysis to determine gene importance.

Perform KEGG pathway enrichment on the top genes.

Visualize the enriched KEGG pathways.

# Dependencies:

ranger: Random Forest implementation.

tidyverse: Data manipulation and visualization.

data.table: Fast data handling.

tidymodels: Machine learning workflows.

clusterProfiler: KEGG pathway enrichment.

org.Mm.eg.db: Mouse gene annotations.

future: Parallel processing for faster execution.

rstudioapi: For selecting files interactively.

# File Structure:

counts_data.csv: A CSV file containing RNA-Seq counts (genes as rows, samples as columns).

treatment_data.csv: A CSV file containing treatment group information (Sample and Outcome columns).

# Treatment Group File:

The treatment file must include two columns:

Sample: Unique sample identifiers (these should match the sample names in the counts data).

Outcome: The experimental group or treatment condition (e.g., Control, Treated).

Example:

Sample,Treatment
Sample1,Control
Sample2,Control
Sample3,Treated
Sample4,Treated

# How to Use:

Install Required Packages: Before running the script, ensure all the necessary R libraries are installed. You can install them using:

install.packages(c("ranger", "tidyverse", "data.table", "tidymodels", "clusterProfiler", "org.Mm.eg.db", "future", "rstudioapi"))

Run the Script: Open an R session and run the script. You will be prompted to:

Select your RNA-Seq counts file (CSV format).

Select the treatment group file (CSV format). This file should be a two column file with sample names that match the counts file in the first column and the treatment groups in the second column.

Enter the outcome variable name (the column name in your treatment group file representing the experimental groups, e.g., "Treatment").

Output:

The script will output the top 200 most important genes based on Random Forest feature importance.

It will also generate a bar plot showing the top 10 enriched KEGG pathways.

How the Code Works:
Data Preprocessing:

RNA-Seq counts are transposed and filtered based on the coefficient of variation of gene expression levels.

Treatment information is merged with the counts data by the sample identifier.

Random Forest Analysis:

The data is split into training and testing sets.

Random Forest models are trained using the outcome variable (e.g., treatment group) and the gene expression data.

Gene importance scores are computed using permutation importance.

KEGG Pathway Enrichment:

The top 200 most important genes are subjected to KEGG pathway enrichment analysis.

The results are visualized using a bar plot showing the top 10 enriched pathways.
