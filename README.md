# Conversion-efficiency-of-genomic-selection

# Simulation Code for Genomic Selection Strategies in Aquaculture Breeding

This repository contains the simulation code associated with the unpublished research paper:

**"Optimizing Genetic Gain and Diversity in Aquaculture Family-Based Breeding Programs for Indirectly Measurable Traits: A Simulation Study of Genomic Selection Strategies."**

## Introduction

This program simulates various genomic selection strategies in family-based aquaculture breeding programs, focusing on traits that are difficult to measure directly in candidates. The simulations aim to optimize genetic gain and maintain genetic diversity by exploring different panel densities, genotype imputation strategies, and reference population sizes.

## Requirements

Before using this code, please ensure that the following software is installed on your system:

- [AlphaMate](https://www.alphagenes.roslin.ed.ac.uk/alphamate/)
- [GCTA](http://cnsgenomics.com/software/gcta/)
- [PLINK](https://www.cog-genomics.org/plink/)
- [HiBLUP](https://github.com/xiaolei-lab/HiBLUP)

Additionally, install the R packages listed in `package.R`.

**Note:** The `ocs.R` script is sourced from the AlphaLearn software.

## Usage

### Directory Naming Convention

The program automatically extracts parameters from the directory name. Please name the directory using the following format:

```
[model][panel_density][I/F][reference_size]
```

**Example:**

```
pop114I10
```

This represents:

- **`pop`**: Using the "pop" model. ("fam" model is unnecessary in this study)
- **`114`**: Panel density of **114 SNPs per chromosome**.
- **`IF`** or **`IT`**: Indicates whether to use genotype imputation.
  - **`IT`**: Using genotype imputation.
  - **`TF`**: Not using genotype imputation.
- **`10`**: Reference population size of 10.

**Additional Examples:**

- **`pop114IF10`**: Using the "pop" model, panel density of 114 SNPs per chromosome, **not** using genotype imputation, reference size of 10.
- **`pop23IT30`**: Using the "pop" model, panel density of **23 SNPs per chromosome**, using genotype imputation, reference size of 30.

### Scripts

- **`createddirandpbs.sh`**: This script automatically creates the directories for all parameter combinations used in the study.

### Setup Instructions

1. **Download All Scripts**: Ensure all scripts are downloaded to the appropriate working directory on your computer.

2. **Set Working Paths**: Update the working paths in all scripts to match the directory structure on your system.

   - Open each script and modify the path variables to reflect your working directory.
   - Example:

     ```sh
     WORK_PATH="/path/to/your/working/directory"
     ```

3. **Run `createddirandpbs.sh`**: Execute this script to generate the necessary directories and scripts for the simulations.

   ```sh
   sh createddirandpbs.sh
   ```

### Notes on Job Submission Scripts

- The `.pbs` files are job submission scripts intended for use on high-performance computing clusters utilizing the Portable Batch System (PBS).

- **If you are not using a PBS-based cluster**, these files may not be necessary, and you may need to adapt the scripts to your scheduling system or run the simulations locally.

## Additional Information

- **Parameter Extraction**: The program reads parameters directly from the directory names, automating the configuration process for different simulation scenarios.

- **Genotype Imputation**: Genotype imputation can significantly impact the results, especially at lower panel densities. Indicate whether to use imputation by including `IT` (imputation) or `IF` (no imputation) in the directory name.

- **Panel Density**: The panel density (`[panel_density]`) should be specified as the number of SNPs per chromosome. This allows the simulation to adjust marker densities accordingly.

## Contact Information

For questions or further assistance, please contact:

- **[Ziyi Kang]**
- **Email**: [kangziyi1998@163.com]
- **Institution**: [Ocean university of China]

- - **[Sheng Luan]**
- **Email**: [luansheng@ysfri.ac.cn]
- **Institution**: [Yellow Sea Fisheries Research Institute, Chinese Academy of Fishery Sciences]

---

By following these instructions, you should be able to replicate the simulations conducted in the study and explore the effects of different genomic selection strategies on genetic gain and diversity in aquaculture breeding programs.
