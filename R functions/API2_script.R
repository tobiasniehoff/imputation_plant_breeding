rm(list = ls())

library(data.table)
library(microbenchmark)
library(stringr)

setwd("~/Studium Agrarwissenschaften In Goettingen/Master/Masterarbeit/R_thesis_tobias")
odir <- getwd()

setwd("Niehoff_Beissinger_Gholami/Niehoff_Beissinger_Gholami/")
getwd()

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

###########
# test imputation for one family
mt <- extract_family_marker_table(pedigree, marker_table, "Fam08")
ped <- extract_individuals_pedigree(pedigree, ID_names = colnames(mt)[-1])

hd <- extract_parents_pedigree(ped)[[1]]
mt <- order_marker_table_by_mapo(mapo = map, marker_table = mt)
mt <- order_marker_table_by_pedigree(pedigree = ped, marker_table = mt)

mt <- mt[mt[[1]] %in% map[map[[2]] != 99,1], ] # remove markers with unknown location

# extract only specific chromosome (here chromosome 1)
mt <- mt[mt[[1]] %in% map[map[[2]] == 1, 1], ]

# cover/set to missing all the marker of non hd individuals that are not on the 
# LD panel
masked_mt <- make_masked_marker_table(mt, map_LD2, HD_indiv = hd)
mapo <- extract_marker_mapo(mapo = map, marker_IDs = masked_mt[[1]])

founders_obj <- make_founders_obj(pedigree = ped) 
genotypes <- convert_marker_table_to_GenotypeFile(marker_table = masked_mt)

##### create a haplotype library for API2
arguments_library <- make_arguments_createlib_API2(out = "lib_fam08"
                                                   , genotypes = "masked.genotypes", 
                                                      hd_threshold = 0.9, n_sample_rounds = 2,
                                                      write_args_to_txt = TRUE
                                                   , n_haplotypes = 10)
gmap <- convert_mapo_to_GeneticMapFileAPI(mapo)

# API2 does not use a map but I need it to loop over the chromosomes
IN_createlib <- createlib_API2(arguments_createlib = arguments_library, genotypes = genotypes
                   ,GeneticMapFile = gmap, verbose = TRUE)
libphase <- IN_createlib[[1]]
run_time <- IN_createlib[[2]]
# btw a library is nothing else than a phase file so you could take the phase 
# information from wherever (BEAGLE) and use it as library

##### after creating the library, you can impute
arguments_impute <- make_arguments_impute_API2(out = "imputed_fam08", libphase = "lib_fam08.phase",
                                               genotypes = "masked.genotypes", 
                                               write_args_to_txt = TRUE, 
                                               founders = "founders_file.txt",
                                               overwrite = TRUE, seed = 48)
# you can also read in an existing arguments file
# arguments_impute <- read_arguments("arguments_impute_API2.txt")
IN_impute <- impute_API2(arguments_impute = arguments_impute, genotypes = genotypes
                     ,founders = founders_obj, GeneticMapFile = gmap, libphase = libphase)
imputed_mt <- IN_impute[[1]]
(run_time <- IN_impute[[2]])

imputed_mt <- order_marker_table_by_pedigree(marker_table = imputed_mt, pedigree = colnames(mt)[-1])

# needed because the parents are not imputed (actually true but I include them 
# in the impute function )
# actually, for a fair comparison, HD individuals (so parents) should be excluded 
# from the analysis
true <- extract_individuals_marker_table(marker_table = mt, ID_inds = colnames(imputed_mt)[-1]) 
true <- order_marker_table_by_pedigree(marker_table = true, pedigree = colnames(imputed_mt)[-1])

anal_inds <- analyze_inds_imputation(true, imputed_mt)
anal_pop <- analyze_pop_imputation(true, imputed_mt)

anal_pop
View(anal_inds)

#####################################################
# imputation with a whole population as library


############### has been changed according to new functions but has not been tested yet
######### DO NOT EXECUTE. TAKES WAY TOO MUCH TIME

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

masked_genotypes <- convert_marker_table_to_GenotypeFile(masked_mt)

#arguments_library <- read_arguments("already_existing_arguments.txt")
arguments_library <- make_arguments_createlib_API2(out = "lib_totpop", genotypes = "masked.genotypes", 
                                                      hd_threshold = 0.9, n_sample_rounds = 5)

gmap <- convert_mapo_to_GeneticMapFileAPI(mapo)

# API2 does not use a map but I need it to loop over the chromosomes
libphase <- createlib_API2(arguments_createlib = arguments_library, genotypes = masked_genotypes
                           ,GeneticMapFile = gmap, verbose = TRUE)

arguments_impute <- make_arguments_impute_API2(out = "imputed", libphase = "lib_totpop.phase",
                                               genotypes = "masked.genotypes")

# you can also read in an existing arguments file
# arguments_impute <- read_arguments("arguments_impute_API2.txt")
imputed_mt <- impute_API2(arguments_impute = arguments_impute, genotypes = masked_genotypes
                          , GeneticMapFile = gmap, libphase = libphase)


imputed_mt <- order_marker_table_by_pedigree(marker_table = imputed_mt, pedigree = colnames(mt)[-1])

# this is done because the parents weren't actually imputed. Only the offspring
parents <- extract_parents_pedigree(pedigree)[[1]]
true <- extract_individuals_marker_table(marker_table = mt, ID_inds = parents, negate = TRUE) 
imputed_mt <- extract_individuals_marker_table(marker_table = imputed_mt, ID_inds = parents, negate = TRUE) 

true <- order_marker_table_by_pedigree(marker_table = true, pedigree = pedigree)
imputed_mt <- order_marker_table_by_pedigree(marker_table = imputed_mt, pedigree = pedigree)

anal_inds <- analyze_inds_imputation(true, imputed_mt)
anal_pop <- analyze_pop_imputation(true, imputed_mt)

#################
# for fun, compare if imputation per family or with the whole population is better