#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
rm(list = ls())
#setwd("~/Studium Agrarwissenschaften In Goettingen/Master/Masterarbeit/R_thesis_tobias")


# load packages
library(data.table)
library(microbenchmark)
library(stringr)

setwd("./R_thesis_tobias/")
odir <- getwd()

setwd("Niehoff_Beissinger_Gholami/Niehoff_Beissinger_Gholami/")

# read in data
marker_table <- fread("Sugar-beet-data1_Marker_012.txt")
pedigree <- fread("Sugar-beet-data1_ped.txt")
pedigree[is.na(pedigree)] <- 0

map <- read.table("Sugar-beet-data1_map.txt", header = TRUE)
map_LD1 <- read.table("Sugar-beet-data1_map-LD1.txt", header = TRUE)
map_LD2 <- read.table("Sugar-beet-data1_map-LD2.txt", header = TRUE)
map_LD3 <- read.table("Sugar-beet-data1_map-LD3.txt", header = TRUE)

# load functions developed for other software into environment
setwd(odir)
setwd("R AlphaPlantImpute/")
source("API_functions.R")

setwd(odir)
setwd("R AlphaPlantImpute2/")
source("API2_functions.R")

setwd(odir)
setwd("R BEAGLE/")
source("analysis.R")
source("BEAGLE_functions.R")
setwd(odir)

setwd("R simulation/")
source("simulation_functions.R")
source("chromoMap_functions.R")
setwd(odir)

##################################################################
# some basic analysis of the data first
ei <- analyze_mendel_nerr_ind(marker_table = marker_table, pedigree = pedigree)
em <- analyze_mendel_nerr_marker(marker_table = marker_table, pedigree = pedigree)
ai <- analyze_inds_marker_table(marker_table)
am <- analyze_marker_marker_table(marker_table)

##################################################################
# Imputation for one family for one chromosome

mt <- extract_family_marker_table(pedigree, marker_table, "Fam08")
ped <- extract_individuals_pedigree(pedigree, ID_names = colnames(mt)[-1])

hd <- extract_parents_pedigree(ped)[[1]]
mt <- order_marker_table_by_mapo(mapo = map, marker_table = mt)
mt <- order_marker_table_by_pedigree(pedigree = ped, marker_table = mt)

mt <- mt[mt[[1]] %in% map[map[[2]] != 99,1], ]

# extract only specific chromosome (here chromosome 2)
mt <- mt[mt[[1]] %in% map[map[[2]] == 2, 1], ]

# cover/set to missing all the marker of non hd individuals taht are not on the 
# LD panel
masked_mt <- make_masked_marker_table(mt, map_LD2, HD_indiv = hd)
mapo <- extract_marker_mapo(mapo = map, marker_IDs = masked_mt[[1]])
##################################################################

##### create a subdirectory for this analysis
dir.create("temporary")
setwd("./temporary")
######################

write_marker_table_to_vcf(marker_table = masked_mt, pedigree = ped, mapo = mapo, 
                          outfile_name = "chrom_02_fam08.vcf", gzipped = FALSE)

map_PLINK <- convert_mapo_to_PLINK_map(mapo = mapo)
write_PLINK_map(PLINK_map_object = map_PLINK, file_prefix = "test_map")

#arguments <- read_arguments("already_existing_BEAGLE_arguments.txt")
arguments <- make_arguments_BEAGLE(gt = "chrom_02_fam08.vcf", out = "out_fam08.gt", 
                                   write_args_to_txt = F, map = "test_map.map")

# 'impute_BEAGLE' doesn't actually produce any output than just the time needed
# the imputed files are written to a file
run_time <- impute_BEAGLE(arguments = arguments, memory_limit_GB = 8, verbose = TRUE)

################################
#### read in results
imputed_mt <- read_vcf_to_marker_table("out_fam08.gt.vcf.gz", set_missing_gt_to = 9)
imputed_phase <- read_vcf_to_phase("out_fam08.gt.vcf.gz", set_missing_gt_to = 9)

#### after reading in the results, the subdirectory and its files can be deleted
setwd("..")
unlink("./temporary", recursive = TRUE)

###########
# analysis
imputed_mt <- order_marker_table_by_pedigree(marker_table = imputed_mt, pedigree = colnames(mt)[-1])

analyze_inds_imputation(mt, imputed_mt)
analyze_pop_imputation(mt, imputed_mt)

analyze_marker_marker_table(mt)
analyze_marker_marker_table(imputed_mt)

analyze_inds_marker_table(mt)
analyze_inds_marker_table(imputed_mt)

# note that BEAGLE does mistakes. So quality filtering and supplying phased data 
# is recommended

##################################################
# for whole population
#extract_family_parents_pedigree(pedigree, per_parent = FALSE)
mt <- marker_table

ped <- extract_individuals_pedigree(pedigree, ID_names = colnames(mt)[-1])

hd <- extract_parents_pedigree(ped)[[1]]
mt <- order_marker_table_by_mapo(mapo = map, marker_table = mt)
mt <- order_marker_table_by_pedigree(pedigree = ped, marker_table = mt)

mt <- mt[mt[[1]] %in% map[map[[2]] != 99,1], ]

# cover/set to missing all the marker of non hd individuals that are not on the 
# LD panel
masked_mt <- make_masked_marker_table(mt, map_LD2, HD_indiv = hd)
mapo <- extract_marker_mapo(mapo = map, marker_IDs = masked_mt[[1]])
#############

##### create a subdirectory for this analysis
dir.create("temporary")
setwd("./temporary")
######################

write_marker_table_to_vcf(marker_table = masked_mt, pedigree = ped, mapo = mapo, 
                          outfile_name = "tot_pop.vcf", gzipped = TRUE)

map_PLINK <- convert_mapo_to_PLINK_map(mapo = map)
write_PLINK_map(PLINK_map_object = map_PLINK, file_prefix = "test_map")

#arguments <- read_arguments("already_existing_BEAGLE_arguments.txt")
arguments <- make_arguments_BEAGLE(gt = "tot_pop.vcf.gz", out = "out_totpop.gt", 
                                   write_args_to_txt = F, map = "test_map.map")

run_time <- impute_BEAGLE(arguments = arguments, memory_limit_GB = 2, verbose = TRUE)

################################
#### read in results
imputed_mt <- read_vcf_to_marker_table("out_totpop.gt.vcf.gz", set_missing_gt_to = 9)
imputed_phase <- read_vcf_to_phase("out_totpop.gt.vcf.gz", set_missing_gt_to = 9)

#### after reading in the results, the subdirectory and its files can be deleted
setwd("..")
unlink("./temporary", recursive = TRUE)

###########
# analysis
imputed_mt <- order_marker_table_by_pedigree(marker_table = imputed_mt, pedigree = colnames(mt)[-1])


ai <- analyze_inds_imputation(mt, imputed_mt)
ap <- analyze_pop_imputation(mt, imputed_mt)

analyze_marker_marker_table(mt)
analyze_marker_marker_table(imputed_mt)

analyze_inds_marker_table(mt)
analyze_inds_marker_table(imputed_mt)

# per marker
a <- analyze_marker_imputation(true_marker_table = mt, 
                          imputed_marker_table = imputed_mt)
b <- analyze_marker_panel_imputation(true_marker_table = mt, 
                                imputed_marker_table = imputed_mt)


##############################################
# some visualization

# this will extract the markers that are informative.
infor <- extract_informative_marker_mapo(mapo, p1_genotype = mt[[2]], p2_genotype = mt[[3]])
# this extracts markers for which the individuals in the marker table are not all the same
# (so it will return only markers for which not all individuals aer 0, 1 or 2 (missing data
# is filtered out in advanced by the function))
infor <- extract_informative_marker_mapo(mapo, marker_table = mt[,c("marker", "P1", "P2")])

# this extracts all the informative markers out of a mapo (e.g. a low density
# SNP panel)
LD_infor <- extract_marker_mapo(mapo = map_LD1, marker_IDs = infor[[1]])

visualize_mapo_chromoMap(HD_mapo = mapo, LD_mapo = mapo, ploidy = 1)
visualize_mapo_chromoMap(HD_mapo = mapo, LD_mapo = map_LD1, ploidy = 1)
visualize_mapo_chromoMap(HD_mapo = mapo, LD_mapo = map_LD2, ploidy = 1)
visualize_mapo_chromoMap(HD_mapo = mapo, LD_mapo = map_LD3, ploidy = 1)
visualize_mapo_chromoMap(HD_mapo = mapo, LD_mapo = LD_infor, ploidy = 1)

# ideally done if you have true none missing genotype data
visualize_imp_errors_mt_chromoMap(mapo = mapo, true_genotype_vector = mt[[5]]
                                  ,imputed_genotype_vector = imputed_mt[[5]])
# this can only be done if you have true phase data for parents
visualize_origin_chomoMap(parents_phase = true_phase[1:4,]
                          , offspring_phase = imputed_phase[7:8,], mapo = mapo
                          ,direct_offspring = FALSE)

#############################################
# you can extract the phase data and use it as a reference for the next imputation round
parents <- extract_parents_pedigree(ped)
parents_phase <- extract_individuals_phase(imputed_phase, ID_inds = hd)

#############################################
# familywise imputation

mt <- extract_family_marker_table(pedigree, marker_table, c("Fam08", "Fam10"))
ped <- extract_individuals_pedigree(pedigree, ID_names = colnames(mt)[-1])

hd <- extract_parents_pedigree(ped)[[1]]
mt <- order_marker_table_by_mapo(mapo = map, marker_table = mt)
mt <- order_marker_table_by_pedigree(pedigree = ped, marker_table = mt)

mt <- mt[mt[[1]] %in% map[map[[2]] != 99,1], ]

# cover/set to missing all the marker of non hd individuals that are not on the 
# LD panel
masked_mt <- make_masked_marker_table(mt, map_LD2, HD_indiv = hd)
mapo <- extract_marker_mapo(mapo = map, marker_IDs = masked_mt[[1]])
#############
setwd(odir)

# this will return a phase object for every family
IN <- impute_familywise_BEAGLE(arguments = arguments, marker_table = masked_mt, pedigree = ped
                              , memory_limit_GB = 8, mapo = mapo, verbose = TRUE)
IN_phase <- IN[[1]]
(run_time <- IN[[2]])
IN_mt <- lapply(IN_phase, convert_phase_to_marker_table)

# IN_mt contains a marker table for every family
# sometimes, imputation results differ per family for the same individual (a parent)
# This function here will return a marker table with only information at those 
# marker-individual combination that were predicted always the same
# others are missing (9)

# this shows you the individual predictions per individual
order_ind_imps_lsmt(IN_mt)
unq_mt <- extract_unique_imp_lsmt(IN_mt)

unq_mt <- order_marker_table_by_mapo(marker_table = unq_mt, mapo = colnames(mt)[-1])

extract_family_parents_pedigree(pedigree, per_parent = TRUE)
