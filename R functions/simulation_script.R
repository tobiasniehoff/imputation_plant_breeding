#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
rm(list = ls())

# load packages
library(data.table)
library(microbenchmark)
library(stringr)
library("Meiosis") # the 'Meiosis package from Dominik Mueller will be used
library(rlist)

setwd("~/Studium Agrarwissenschaften In Goettingen/Master/Masterarbeit/R_thesis_tobias")
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

setwd("R simulation/")
source("simulation_functions.R")
setwd(odir)

setwd(odir)
setwd("R BEAGLE/")
source("analysis.R")
source("BEAGLE_functions.R")
setwd(odir)

########
vignette('Introduction', package = 'Meiosis')
## Use Meiosis::seed_rng to seed the RNG of the C++ routines (Mersenne Twister)

##########
tempdir()
#########

set.seed(123L) ## Seed R's rng
Meiosis::seed_rng(seed = 123L) ## Seed rng used by Meiosis
######################

#############################################################################
# apparently you can have NAs or any other allele coding than 0,1, e.g. 2 or 3 as long as
# it is integer
#parents[[1]]$maternal[[10]] <- as.integer(rep(2, length(parents[[1]]$maternal[[10]])))
#parents[[1]]$paternal[[10]] <- as.integer(rep(3, length(parents[[1]]$paternal[[10]])))
#############################################################################

#marker_table <- marker_table[,1:20]
n_chr <- unique(map[[2]])

xoparam <- convert_mapo_to_xoparam(map)
map <- extract_marker_mapo(mapo = map, marker_IDs = marker_table[[1]])
positions <- extract_positions_mapo(mapo = map)
check_positions(positions)
marker_table <- order_marker_table_by_mapo(mapo = map, marker_table = marker_table)
geno <- convert_marker_table_to_geno(marker_table)
phase <- convert_geno_to_phase(geno)
phase[phase == 9] <- 0
true_pop <- convert_phase_to_meiosis_pop(phase, mapo = map)
new_phase <- convert_meiosis_pop_to_phase(meiosis_pop = true_pop)
colnames(new_phase)[-1] <- map[[1]]
identical(new_phase, phase)
#####################################
pop <- true_pop
F1_pop <- breed_cross_inds(ind1 = pop[[1]], ind2 = pop[[2]], n_offspring = 100, 
                         positions = positions, xoparam = xoparam)
S3_pop <- breed_self_pop(F1_pop, n_selfings = 10, positions = positions, xoparam = xoparam)
DH_pop <- breed_make_dh_pop(F1_pop, positions = positions, xoparam = xoparam)
rndmat_pop <- breed_random_mating_pop(F1_pop, n_random_mating = 10, positions = positions, xoparam = xoparam)

analyze_marker_marker_table(convert_meiosis_pop_to_marker_table(pop))
analyze_marker_marker_table(convert_meiosis_pop_to_marker_table(F1_pop))
analyze_marker_marker_table(convert_meiosis_pop_to_marker_table(S3_pop))
analyze_marker_marker_table(convert_meiosis_pop_to_marker_table(DH_pop))
analyze_marker_marker_table(convert_meiosis_pop_to_marker_table(rndmat_pop))

analyze_inds_marker_table(convert_meiosis_pop_to_marker_table(pop))
analyze_inds_marker_table(convert_meiosis_pop_to_marker_table(F1_pop))
analyze_inds_marker_table(convert_meiosis_pop_to_marker_table(S3_pop))
analyze_inds_marker_table(convert_meiosis_pop_to_marker_table(DH_pop))
analyze_inds_marker_table(convert_meiosis_pop_to_marker_table(rndmat_pop))

# for getting a homozygousity parameter for teh whole population
mean(analyze_inds_marker_table(convert_meiosis_pop_to_marker_table(rndmat_pop))[[3]])


set.seed(123L) ## Seed R's rng
Meiosis::seed_rng(seed = 123L) ## Seed rng used by Meiosis