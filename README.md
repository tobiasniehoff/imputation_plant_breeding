# Low-density marker imputation in plant breeding populations

This is the repository for scripts used in our low-density marker imputation manuscript. If more details are needed, please contact [Tobias](tobias.niehoff@wur.nl).

## data
The folder `data` contains a simple simulated data set to test the workflows. The data set does not contain genetic data as used in our study. It merely serves to provide input data of the same format as used in our study. The data set was simulated with MoBPS (https://github.com/tpook92/MoBPS). 
In brief, simulation was done as follows:
 - 63000 SNPs were simulated over 9 chromosomes each with length of 0.8 Morgan
 - 2700 SNPs were assigned a QTL effect sampled from a Gamma distribution with shape parameter 0.4
 - 1000 founders a simulated
 - 100 distinct generations of random mating starting with the founders were simulated with a linearly decreasing population size to 100 after 100 generations
 - next, the population expanded rapidly linearly from 100 to 1000 plants in 4 generations
 - next, 10 generations of breeding activities are simulated:
   - based on phenotypic selection
   - with heritability of 0.5
   - selection of the 30 best plants out of 2000
   - selection based on F3 plants
   - random mating
 - data set generation:
   - 20,000 markers with a minor allele frequency over 0.01 were randomly sampled and used for the 20k array
   - 3,000 out of the 20,000 SNPs were randomly sampled for the 3k array
   - 2,000 out of the 3,000 SNPs  were randomly sampled for the 2k array
   - 1,000 out of the 2,000 SNPs  were randomly sampled for the 1k array
   - 30 plants of the last generation were randomly sampled and mated to obtain F1 plants
   - these were selfed twice to get F3 plants (these are used as parents later)
   - the workflows simulate offspring based on these parents anew but we also provide 50 simulated F3 offspring lines per family in the data set

**make_simulated_data.R** is an R script to produce the simulated data mimicking sugar beet with the R package MoBPS.

**simulated_data.vcf.gz** is an vcf file with true full and phased marker data of 36 families each with 50 offspring lines.

**simulated_data_ped.txt** is a text file that denotes the pedgree of the individuals.

**simulated_map.txt** is a text file with the genetic map of all 20000 markers on all 9 chromosomes.

**simulated_map3.txt** is a text file with the genetic map of 3000 markers that are on the 3k low-density panel.

**simulated_map2.txt** is a text file with the genetic map of 2000 markers that are on the 2k low-density panel.

**simulated_map1.txt** is a text file with the genetic map of 1000 markers that are on the 1k low-density panel.

## imputation workflows
The folder `imputation workflows` contains R scripts with which imputation can be done from within R. The scripts are reading in a parameter table that is stored in the folder `input parameter tables`.The scripts were written to work on a cluster with the lines in the input parameter table executed in parallel. Users without access to a cluster can test the scripts by modifying the R scripts at the positions indicated with **MODIFY_HERE** in the code.
The imputation workflows will create a new folder for every replicate (row in excel input parameter tables). After a run is finished, its directory is deleted.
The output of the workflows are `.txt` and `.xlsx` files (read section **Quality measures**).

### R scripts
The R scripts will create directories. Results are stored in a directory with the same name as the input parameter file and a time stamp.

**API2_workflow.R** is the R script to use if imputing with AlphaPlantImpute2.

**BEAGLE_workflow.R** is the R script to use if imputing with Beagle version 5.1.

**BEAGLE52_workflow.R** is the R script to use if imputing with Beagle version 5.2.

**BEAGLE53_workflow.R** is the R script to use if imputing with Beagle version 5.3.

**single-genotype_BEAGLE_workflow.R** is the R script to use if using the single-genotype imputation approach with Beagle version 5.1.

**BEAGLE_API2_workflow.R** is the R script to use for the Beagle+API2 approach. Beagle version 5.1 is used.


### Finding software
The R scripts using AlphaPlantImpute2 expect that AlphaPlantImpute2 is added to the `$PATH` variable. AlphaPlantImpute2 can be downloaded [here](https://github.com/AlphaGenes/AlphaPlantImpute2). Make sure you have at least python version 3.7 and recent versions of NumPy and Numba installed. Instructions how to add a directory to the `$PATH` variable can be found [here](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/).

The R scripts using Beagle will search for the Beagle file in the directory of the script where it is executed and in the parent directory. Most straightforward is to store the Beagle files in the folder `imputation workflow`. Beagle 5.1 (18May20.d20) can be downloaded [here](https://faculty.washington.edu/browning/beagle/b5_1.html), Beagle 5.2 (28Jun21.220) [here](https://faculty.washington.edu/browning/beagle/old.beagle.html) and Beagle 5.3 (08Feb22.fa4) can be downloaded [here](https://faculty.washington.edu/browning/beagle/beagle.html).

### SLURM jobs
The `.sh` scripts contain the job description for SLURM jobs. The parameter table to use (stored in folder `input parameter tables` has to be specified in the shell scripts.

**API2_workflow.sh** is the SLURM job description for AlphaPlantImpute2 imputation.

**BEAGLE_workflow.sh** is the SLURM job description for Beagle version 5.1 imputation.

**BEAGLE52_workflow.sh** is the SLURM job description for Beagle version 5.2 imputation.

**BEAGLE53_workflow.sh** is the SLURM job description for Beagle version 5.3 imputation.

**single-genotype_BEAGLE_workflow.sh** is the SLURM job description for the single-genotype imputation approach with Beagle.

**BEAGLE_API2_workflow.sh** is the SLURM job description for the Beagle+API2 approach.

## input parameter tables
Input parameter tables are stored in the folder with the same name. The procedure and the input parameters can be specified with the input parameter tables in excel.  Names of the tables are self-explanatory. Functionality of parameters is explained in the excel files. 
The parameters can in principle specified to perform tests not presented in our manuscript. 
For instance, imputation of low-coverage data sets can be simulated by using homozygous parent and offspring population types (e.g., DH or F\[high number\]), genotyping all offspring and parents at high-density and setting the fraction of randomly missing calls per individual to the desired fraction (e.g., 0.5).
Simple quality control can be used by specifying minimum values for the minor allele frequency (MAF) and/or minor allele count (MAC). These filteres are applied to all lines before any phasing or imputation.

## R functions
This folder stored R scripts in which functions are defined that are sourced by the workflows.

**analysis.R** contains functions that analyze data sets.

**API2_functions.R** contains functions specific for AlphaPlantImpute2 workflows.

**BEAGLE_functions.R** contains functions that are specific for Beagle workflows and for reading and writing .vcf files. The imputation functions are written to only work on Windows or Linux systems. If you are on a different operation system, modify the function `impute_BEAGLE[5.1|5.3]()`.

**general_functions.R** contains functions that are used by all scripts.

**simulation_functions.R** contains functions that are used for simulation of data in the workflows using the R package `Meiosis`.

## Quality measures
The file **parameter_finding3.xlsx** is empty and only contains the names for the quality parameters. The accuracy reported in the manuscript is written in the column `GenotypeMeanAccuracy_whole_pop_all_loci` of the output tables created by the workflows. The output tables are written in `.txt` format and `.xlsx` in a new directory created by the workflows with the same name as the input parameter table. The new directory will be created in the directory named `imputation workflows`.
