# Original Paper Summary

**Title**: Soil microbiome indicators can predict crop growth response to large-scale inoculation with arbuscular mycorrhizal fungi  
**Authors**: Lutz et al.  
**Year**: 2023  
**Journal**: Nature Microbiology

**Biological Question**:  
Can the composition of the soil fungal microbiome predict the success of arbuscular mycorrhizal fungi (AMF) inoculation in enhancing maize (Zea mays L.) growth under field conditions?

**Key Findings**:
- AMF inoculation responses in 54 Swiss maize fields varied widely, from -12% to +40% in yield response.
- Only ~25% of the fields showed a statistically significant positive growth response.
- Soil fungal OTUs (sOTUs), especially from ITS long-read sequencing, were much better predictors of positive response than conventional soil parameters.
- A final regression model including 15 soil variables and 13 sOTUs explained up to 86% of MGR variation—fungal indicators alone accounted for 53%.
- Mechanism of yield increase was suppression of soil-borne fungal pathogens, not improved nutrient acquisition.
- Fields that responded positively had more pathogenic fungi initially, which were suppressed post-inoculation.
- AMF colonization rate did not correlate with crop response.

**Conclusion**:  
Soil microbiome profiles, particularly fungal community structure, are strong predictors of AMF inoculation success and can guide microbial-based crop management.

## Replication Goals
We aimed to evaluate the general reproducibility of the original analysis workflow, focusing on:
- Figure 2, Mycorrhizal growth response over 3 yr.
- Portions of figure 3, Variable selection for the MGR prediction model.

## Methods Used in Replication
**Platform**:  
- Personal laptops and ISU HPC. 
- RStudio, locally and on Nova interactive.

**Tools and Packages**:  
- The names of the packages are: cowplot, DESeq2, dplyr, fantaxtic, ggfortify, MuMln, ggplot2, ggpubr, ggvegan, indicspecies, lme4, MASS, microbiome, olsrr, phyloseq, psych, randomForest, relaimpo, reshape2, scales, stats, svglite, tidyverse, vegan, and VennDiagram.

**Steps Taken**:
1. Downloaded sequencing and soil data from the European Nucleotide Archive.
2. Downloaded the reference genome from UNITE fungal database.
3. Inspected the data.
4. Ran the code provided in the authors' GitHub page but encoutered several challenges.

## Results
- **Comparison with Original Study**:  
  - Our results closely matched the original figures.
  - Discrepancies were due to slight differences in preprocessing steps and package versions.

## Technical Details of Reproducibility
- We attempted to replicate the original PacBio demultiplexing pipeline but found that essential input files (.subreads.bam and barcode references) were not available from the ENA database. Therefore, we proceeded directly with the available demultiplexed .fastq files.
- We restructured the directory to match the expectations of the provided R script (ASV_PacBio.R) and ran the analysis in an R 4.3 environment on ISU's NOVA system. We installed necessary packages (dada2, phyloseq, DECIPHER) via Bioconductor and addressed issues with hardcoded paths, invalid job submission scripts, and deprecated functions.
- Using the UNITE 2021 database for reference, we successfully generated ASV tables and plots replicating key aspects of Figure 2. Despite several reproducibility challenges—such as inconsistent documentation, missing files, and poorly structured code—we were able to execute the core downstream analysis.
