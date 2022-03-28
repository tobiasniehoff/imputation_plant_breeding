rm(list = ls())

# load packages
library(data.table)
library(microbenchmark)
library(stringr)
library(rlist)
suppressMessages(library("Meiosis"))
library(readxl)


# this loop is searching for a folder with the ending 'imputation workflow' in 
# the current working directory and up to two parent directories
for(i in 1:3) {
  ld <- list.dirs(recursive = TRUE, full.names = TRUE)
  workflow_dir <- ld[str_detect(ld, pattern = "/imputation workflow$")]
  if(length(workflow_dir) == 0) {setwd("..")} 
  else {break}
}
setwd(workflow_dir)
(odir <- getwd())

# load in functions
setwd("../R functions")
source("analysis.R")
source("general_functions.R")
source("simulation_functions.R")
source("BEAGLE_functions.R")
setwd(odir)


args_script <- commandArgs(trailingOnly = TRUE)
index_for_scenario <- as.numeric(args_script[1]) 
i <- index_for_scenario
# name_parameter_table <- "BEAGLE53_example_file" ### MODIFY_HERE: uncomment the line
name_parameter_table <- args_script[2]
n_cores <- as.numeric(Sys.getenv(c("SLURM_CPUS_PER_TASK"))) ### MODIFY_HERE: specify number of cores
memory_limit_GB <- as.numeric(Sys.getenv(c("SLURM_MEM_PER_NODE")))/1024### MODIFY_HERE: specify memory in GB

# read in the parameter table
setwd("../input parameter tables")

suppressMessages(
  col.names <- data.frame(readxl::read_excel(paste0(name_parameter_table, '.xlsx'), sheet = 1, skip = 0))
)
# there are two columns with this name, one for phasing, one for imputation
# I want to detect the first one
# I don't want to hard code it cause the table layout may change later
start.imputation.parameters <- min(which(stringr::str_detect(col.names[1,], pattern = 'per_family_imputation')))

col.names <- c(col.names[1,1:(start.imputation.parameters-1)], 
               paste0(col.names[1,start.imputation.parameters:(start.imputation.parameters+15)], '_imputing'), 
               paste0(
                 col.names[1,(start.imputation.parameters+15+1):(start.imputation.parameters+15+1+15)], '_phasing'))
col.names <- stringr::str_replace_all(col.names, pattern = '[/-]', '.')

suppressMessages(
  para_tab <- data.frame(readxl::read_excel(paste0(name_parameter_table, '.xlsx'), sheet = 1, skip = 1))
)

colnames(para_tab) <- col.names

para_tab$random_offspring <- as.logical(para_tab$random_offspring)
para_tab$parents_HD <- as.logical(para_tab$parents_HD)
para_tab$provide_parents_as_ref <- as.logical(para_tab$provide_parents_as_ref)
para_tab$provide_HD_inds_as_ref <- as.logical(para_tab$provide_HD_inds_as_ref)
para_tab$provide_HD_fam_members_as_ref <- as.logical(para_tab$provide_HD_fam_members_as_ref)
para_tab$provide_parents_as_ref <- as.logical(para_tab$provide_parents_as_ref)
para_tab$HD_inds_in_gt_set <- as.logical(para_tab$HD_inds_in_gt_set)

para_tab$per_family_imputation_imputing <- as.logical(para_tab$per_family_imputation_imputing)
para_tab$geneticMap_imputing <- as.logical(para_tab$geneticMap_imputing)
para_tab$BEAGLE_ref_phased_imputing <- as.logical(para_tab$BEAGLE_ref_phased_imputing)

para_tab$per_family_imputation_phasing <- as.logical(para_tab$per_family_imputation_phasing)
para_tab$geneticMap_phasing <- as.logical(para_tab$geneticMap_phasing)
para_tab$`BEAGLE_ref_phased(no_effect)_phasing` <- as.logical(para_tab$`BEAGLE_ref_phased(no_effect)_phasing`)

para_tab$BEAGLE_em_imputing <- str_to_lower(para_tab$BEAGLE_em_imputing)
para_tab$BEAGLE_em_phasing <- str_to_lower(para_tab$BEAGLE_em_phasing)
################################################################################
### MODIFY_HERE: for (i in 1:nrow(para_tab)){
setwd(odir)
setwd("../data/")

map <- read.table('simulated_map.txt', header = TRUE)
pedigree <- fread('simulated_data_ped.txt')
imputed_phase <- read_vcf_to_phase("simulated_data.vcf.gz", set_missing_gt_to = 9)

map_LD1 <- read.table('simulated_map1.txt', header = TRUE)
map_LD2 <- read.table('simulated_map2.txt', header = TRUE)
map_LD3 <- read.table('simulated_map3.txt', header = TRUE)

offspring_type <- para_tab$offspring_type[i]
parent_type <- para_tab$parent_type[i]
n_off.family <- para_tab$offspring.family[i]

sim_data <- make_simulated_data(offspring_type = offspring_type, 
                                parent_type = parent_type,
                                n_off.family = n_off.family, 
                                phase = imputed_phase, 
                                pedigree = pedigree, 
                                map = map, 
                                seed = para_tab$seed_sampling[i])

pedigree <- sim_data$sim_ped
phase <- sim_data$sim_phase
marker_table <- convert_phase_to_marker_table(phase)

rm(sim_data); rm(imputed_phase)
################################################################################
setwd('..')

suppressMessages(
  performance_nam <- readxl::read_excel('parameter_finding3.xlsx', sheet = 1, skip = 1
                                        , col_names = TRUE)
)
performance_nam <- str_remove_all(colnames(performance_nam), pattern = '\\.*[0-9]*$')
performance_nam[2:7] <- paste0(performance_nam[2:7], '_whole_pop_var_loci')
performance_nam[8:13] <- paste0(performance_nam[8:13], '_per_fam_var_loci')
performance_nam[c(14:24, 39:41)] <- paste0(performance_nam[c(14:24, 39:41)], '_whole_pop_all_loci')
performance_nam[c(25:35, 42:44)] <- paste0(performance_nam[c(25:35, 42:44)], '_whole_pop_var_loci')

performance_df <- data.frame(matrix(NA, nrow = nrow(para_tab), ncol = length(performance_nam)))
colnames(performance_df) <- performance_nam
performance_df[[1]] <- para_tab[[1]]

################################################################################

marker_table <- order_marker_table_by_pedigree(pedigree = pedigree, 
                                               marker_table = marker_table)
pedigree <- extract_individuals_pedigree(pedigree = pedigree, 
                                         ID_names = colnames(marker_table))
marker_table <- extract_individuals_marker_table(marker_table = marker_table, 
                                                 ID_inds = pedigree[[1]])
mapo <- map[map[[2]] != 99,]
mapo <- extract_marker_mapo(mapo = mapo, marker_IDs = marker_table[[1]])
marker_table <- extract_marker_marker_table(marker_table = marker_table, 
                                            marker_IDs = mapo[[1]])
marker_table <- order_marker_table_by_mapo(mapo = mapo, marker_table = marker_table)
################################################################################

# The functions used in this script produce new directories and delete them all the time.
# This makes running two scripts in parallel difficult as they will delete each 
# others directories.
# Therefore, I create a unique directory (theoretically problems could still 
# occur but with an extremely low probability)

new_dir_name <- 
  paste0('temp_',name_parameter_table,'_run',i,'_')
for(number in 1:10000){
  if(dir.exists(new_dir_name)){new_dir_name <- paste0(new_dir_name, number)}
  else{
    dir.create(new_dir_name)
    break
  }
}
cat(paste0('Will work in directory ', new_dir_name,'.\n'))
setwd(paste0('./',new_dir_name))
################################################################################


# make arguments for BEAGLE phasing
arguments_phasing <- make_arguments_BEAGLE5.3(gt = "tot_pop.vcf.gz", out = "out_totpop.gt", 
                                              write_args_to_txt = F, map = "test_map.map"
                                              , burnin = para_tab$BEAGLE_burnin_phasing[i]
                                              , seed = para_tab$seed_phasing[i]
                                              , phase_states = para_tab$BEAGLE_phase.states_phasing[i]
                                              , imp_states = para_tab$BEAGLE_imp.states_phasing[i]
                                              , imp_segment = para_tab$BEAGLE_imp.segment_phasing[i]
                                              , imp_step = para_tab$BEAGLE_imp.step_phasing[i]
                                              , imp_nsteps = para_tab$BEAGLE_imp.nsteps_phasing[i]
                                              , cluster = para_tab$BEAGLE_cluster_phasing[i]
                                              , ne = para_tab$BEAGLE_ne_phasing[i]
                                              , window = para_tab$BEAGLE_window_phasing[i]
                                              , iterations = para_tab$BEAGLE_iterations_phasing[i]
                                              , overlap = para_tab$BEAGLE_overlap_phasing[i]
                                              , em = para_tab$BEAGLE_em_phasing[i]
)


# make arguments for BEAGLE imputation
arguments_imputing <- make_arguments_BEAGLE5.3(gt = "tot_pop.vcf.gz", out = "out_totpop.gt", 
                                               write_args_to_txt = F, map = "test_map.map"
                                               , ref = 'reference.vcf.gz'
                                               , burnin = para_tab$BEAGLE_burnin_imputing[i]
                                               , seed = para_tab$seed_imputing[i]
                                               , phase_states = para_tab$BEAGLE_phase.states_imputing[i]
                                               , imp_states = para_tab$BEAGLE_imp.states_imputing[i]
                                               , imp_segment = para_tab$BEAGLE_imp.segment_imputing[i]
                                               , imp_step = para_tab$BEAGLE_imp.step_imputing[i]
                                               , imp_nsteps = para_tab$BEAGLE_imp.nsteps_imputing[i]
                                               , cluster = para_tab$BEAGLE_cluster_imputing[i]
                                               , ne = para_tab$BEAGLE_ne_imputing[i]
                                               , window = para_tab$BEAGLE_window_imputing[i]
                                               , iterations = para_tab$BEAGLE_iterations_imputing[i]
                                               , overlap = para_tab$BEAGLE_overlap_imputing[i]
                                               , em = para_tab$BEAGLE_em_imputing[i]
)

if(para_tab$BEAGLE_ref_phased_imputing[i] == FALSE){
  arguments_imputing <- stringr::str_remove(arguments_imputing, pattern = ' ref=.+gz')
}

################################################################################
if(para_tab$map_LD[i] == 1){LD_map <- map_LD1}
if(para_tab$map_LD[i] == 2){LD_map <- map_LD2}
if(para_tab$map_LD[i] == 3){LD_map <- map_LD3}

# the below is to extract the correct number of individuals
fam_names <- unique(pedigree[[2]])
n_off <- para_tab$offspring.family[i]

# this line finds out where there first individual for a certain family is
# this assumes that all members of a family are right below each other
# don't worry about the parents for now
first_ind <- match(fam_names, pedigree[[2]])

# there are a couple of samplings involved in the script. Here, I am setting a seed
set.seed(para_tab$seed_sampling[i])

# testing if offspring should be drawn randomly or not
if(para_tab$random_offspring[i]){
  # -1 cause of 1 off error
  # drawing random offspring
  pos <- c()
  for(j in 1:length(fam_names)){
    if(j == 1){next} ################################### this here is done to exclude 
    # parents from sampling. Sampling parents is actually not a problem since they are 
    # added afterwards anyway but if the sample size is larger than the number of
    # parents present, the sample() function complains
    #
    # be careful if the to tested population doesn't have parents present
    pos <- c(pos, sample(
      seq(from = first_ind[j], 
          to = first_ind[j]+sum(pedigree[[2]] == fam_names[j])-1), 
      size = n_off, 
      replace = FALSE)
    )
  }
} else{
  # if offspring should not be drawn randomly, take the first 5, 10, whatever
  first_ind <- first_ind[-1] # this line is included to exclude all parents (who are the first family)
  pos <- unlist(lapply(first_ind, function(j){seq(from = j, to = (j+n_off-1))}))
}

parents <- extract_parents_pedigree(pedigree = pedigree)

# this is kinda hardcoded as the parents of the OP population are unknown but
# I know which ones where in the original mixture
# so what this condition asks for is if there is an 'OP' for open pollination 
# in the directory name BUT NOT an 'perfam' meaning 'per family open pollination'.
if(str_detect(para_tab$offspring_type[i], 
              pattern = 'OP') & str_detect(para_tab$offspring_type[i], 
                                           pattern = 'perfam', negate = TRUE)){
  parents <- pedigree[1:30,]
}

# so the parents will be included in the final pedigree
ped <- rbind(parents, pedigree[pos,])

# interestingly unique works, so no chaining with 'duplicated' is needed
ped <- unique(ped) 

# not sure if it is needed to store the pedigree (i.e. whether the variable ped
# is changed in later lines), but just to be sure
ped_store <- ped
mt <- extract_individuals_marker_table(marker_table = marker_table, 
                                       ID_inds = ped[[1]])
################################################################################
# this bit was not needed in previous scripts because the provided marker table
# was always the true one. But now since we are inducing errors, this is not the 
# case anymore. mt_true is needed for the analysis in the end
mt_true <- mt
################################################################################
##
# now inducing some genotyping errors
# note that the below function produces errors independent of the state of the 
# truth. I.e. a '0' is as likely as a '1' if the truth is a '2'. In reality,
# it would probably be more likely to get a '1' if the truth is a '2'
mt <- make_genotype_marker_errors_marker_table(marker_table = mt, 
                                               fraction_incorrect = para_tab$frac_incorrect_marker.individual[i])

##
# now inducing some missing values
mt <- make_genotype_marker_missing_marker_table(marker_table = mt, 
                                                fraction_missing = para_tab$frac_NA_marker.individual[i])

##
HD_inds <- c()
# now choosing the individuals that should be given as HD individuals
if(para_tab$parents_HD[i]){
  HD_inds <- extract_parents_pedigree(ped)[[1]]
  
  if(str_detect(para_tab$offspring_type[i], 
                pattern = 'OP') & str_detect(para_tab$offspring_type[i], 
                                             pattern = 'perfam', negate = TRUE)){
    parents <- pedigree[1:30,]
    HD_inds <- parents[[1]]
  }
}
# here I randomly draw offspring from the families that should be HD
# btw, has no affect if the parameter is 0
if(str_detect(para_tab$offspring_type[i], 
              pattern = 'OP') & str_detect(para_tab$offspring_type[i], 
                                           pattern = 'perfam', negate = TRUE)){
  ped_off_temp <- ped[31:nrow(ped),]
  HD_inds <- c(HD_inds, 
               extract_random_HD_offspring_per_fam(pedigree = ped_off_temp, 
                                                   fraction_HD_offspring_per_family = para_tab$frac_offspring_HD_per_family[i])
  )
} else{
  HD_inds <- 
    c(HD_inds, 
      extract_random_HD_offspring_per_fam(pedigree = ped, 
                                          fraction_HD_offspring_per_family = para_tab$frac_offspring_HD_per_family[i])
    )
}

##
HD_inds <- unique(HD_inds)

masked_mt <- make_masked_marker_table(HD_marker_table = mt, LD_map = LD_map, HD_indiv = HD_inds)
################################################################################
## masked_mt is the information that one would get to analyze
# we can also apply filters on it, like filtering for MAF, MAC, or NA

## here, I am throwing out all markers that have too many missing values
# temp_mt enthaelt nur HD individuals. This is based only on the HD individuals
# because the LD individuals already have very many missing information
# since they are genotyped with a much smaller marker panel
temp_mt <- extract_individuals_marker_table(marker_table = masked_mt, ID_inds = HD_inds)
temp_analysis <- analyze_marker_marker_table(marker_table = temp_mt)
# here I am keeping all the markers that are below the threshold for missing markers
masked_mt <- masked_mt[temp_analysis$freq_NA < para_tab$NA_threshold_HD_inds[i],]
rm(temp_mt)
rm(temp_analysis)

mapo <- extract_marker_mapo(mapo = mapo, marker_IDs = masked_mt[[1]])
LD_map <- extract_marker_mapo(mapo = LD_map, marker_IDs = mapo[[1]])

# here filtering for MAF is done.
# Markers with a low MAF are not thrown out but the individuals that have
# such an allele are set to missing
analysis <- analyze_marker_marker_table(masked_mt)

m <- masked_mt
m[m == 9] <- NA
# this is checking which allele the minor allele is
minor_allele <- abs(round(apply(m[,-1], MARGIN = 1, function(d){mean(d, na.rm = TRUE)})/2)-1)

# here filtering is done for MAF
# the way it works is like this:
# it is checked if at a locus the MAF is below the threshold
# if yes, it is checked whether the minor allele is 0,
# if yes, all individuals that have at least one 0 allele meaning their genotype
# is either 0 or 1, are set to missing (9) at this locus
# if the minor allele is 1, then the same as explained about is done for all
# individuals that are either 1 or 2
# for this short while, masked_mt must be changed to data.frame
masked_mt <- as.data.frame(masked_mt)
for(x in 1:nrow(masked_mt)){
  # is TRUE if the MAF is below the threshold
  if(analysis$MAF[x] < para_tab$MAF_threshold_whole_pop[i]){
    # if the minor allele is 0, it is present in 0 and 1. These should be set to 
    # 9 (=missing) then
    if(minor_allele[x] == 0){
      masked_mt[x, masked_mt[x,] %in% c(0,1)] <- 9
    } else{
      # else, so the other case, is when the minor allele is 1
      masked_mt[x, masked_mt[x,] %in% c(1,2)] <- 9
    }
  }
}

# here filtering is done for MAC
for(x in 1:nrow(masked_mt)){
  # is TRUE if the MAF is below the threshold
  if(analysis$MAC[x] < para_tab$MAC_threshold_whole_pop[i]){
    # if the minor allele is 0, it is present in 0 and 1. These should be set to 
    # 9 (=missing) then
    if(minor_allele[x] == 0){
      masked_mt[x, masked_mt[x,] %in% c(0,1)] <- 9
    } else{
      # else, so the other case, is when the minor allele is 1
      masked_mt[x, masked_mt[x,] %in% c(1,2)] <- 9
    }
  }
}
# masked_mt is converted back to data.table because that is the format 
# all other functions require it to be 
masked_mt <- as.data.table(masked_mt)

# who knows what this is good
masked_mt_store <- masked_mt
################################################################################
# whether phasing shall be done per family 
runtime <- NULL
# the line below is complicated but I didn't want to change the whole code cause 
# it grew like this
if(para_tab$per_family_imputation_phasing[i] & para_tab$BEAGLE_ref_phased_imputing[i]){
  # this creates a list of marker tables
  ls_masked_mt <- 
    extract_family_marker_table(pedigree = ped, marker_table = masked_mt)
  
  # now phasing (and therefore also imputation) is done for every family individually
  runtime <- NULL
  for(k in 1:length(ls_masked_mt)){
    ped_fam <- extract_individuals_pedigree(pedigree = ped, 
                                            ID_names = colnames(ls_masked_mt[[k]]))
    # the 'parents' population is confusing
    # they are not related so API cannot impute them, still I don't want to 
    # filter them out by name but use a general rule
    
    # here I check which individuals are parents to others
    ped_fam_p <- extract_parents_pedigree(pedigree = ped_fam)
    # who is an offspring of these parents?
    ped_fam_o <- extract_offspring_pedigree(pedigree = ped_fam, 
                                            parents_ID = ped_fam_p[[1]])
    
    # combine them
    ped_fam <- rbind(ped_fam_p,ped_fam_o)
    # because an individual could be parent AND offspring, I need to use 'unique'
    # (not sure if API actually would accept such individuals)
    ped_fam <- unique(ped_fam)
    
    # this bit is important: if less than 3 individuals are present in the pedigree,
    # so not even a trio, skip this round of the for loop
    # of course it is possible that only one parent is known and API wouldn't 
    # accept it but it works for the structure we have
    if(nrow(ped_fam) <= 2) {
      cat(paste0('Skipped family ', k,' in trial ',i,'.\n'))
      next
    }
    # this is to set the whole pedigree to the same family
    ped_fam[[2]] <- ped_fam[[2]][1]
    
    ########
    # maybe you want to do phasing for a different number of individuals than
    # for the imputation step
    # this block here takes care of that
    # 'extract_random_HD_offspring_per_fam' can be used here even though it 
    # has the 'HD' in its name
    o_id <- extract_random_HD_offspring_per_fam(pedigree = ped_fam, 
                                                fraction_HD_offspring_per_family = para_tab$fraction_offspring_provided_for_phasing[i])
    o_id <- o_id[match(ped_fam_o[[1]], o_id)]
    o_id <- o_id[!is.na(o_id)]
    rand_p <- sample(ped_fam_p[[1]], replace = FALSE,
                     size = ceiling(
                       para_tab$fraction_parents_provided_for_phasing[i]*length(ped_fam_p[[1]])
                     ))
    
    rand_p <- rand_p[match(ped_fam_p[[1]], rand_p)]
    rand_p <- rand_p[!is.na(rand_p)]
    ped_fam_phasing <- extract_individuals_pedigree(pedigree = ped_fam, 
                                                    ID_names = unique(c(rand_p, o_id)))
    masked_mt_for_phasing <- extract_individuals_marker_table(marker_table = ls_masked_mt[[k]], 
                                                              ID_inds = ped_fam_phasing[[1]])
    ########
    
    # here, it is just phased
    # it is done only for one family
    IN <- impute_familywise_BEAGLE5.3(arguments = arguments_phasing
                                   , marker_table = masked_mt_for_phasing
                                   , pedigree = ped_fam
                                   , memory_limit_GB = memory_limit_GB
                                   , mapo = mapo
                                   , verbose = FALSE)
    run_time <- IN$run_time
    ref_phase <- IN$out_phase[[1]]
    cat(
      paste0('\tFinished phasing for family ',k,' of ',
             length(ls_masked_mt),' families in trial ',i,'.\n')
    )
    
    # here I am checking which individuals should be given as reference
    ref_inds <- c()
    if(para_tab$provide_parents_as_ref[i]){
      ref_inds <- ped_fam_p[[1]]
    } else{
      ref_phase <- extract_individuals_phase(phase = ref_phase, 
                                             ID_inds = ped_fam_p[[1]], 
                                             negate = TRUE)
    }
    # it is important to distinguish between ref_inds and libphase here
    # ref_inds will be the individuals that will be given as reference
    # If para_tab$provide_parents_as_ref is FALSE, I don't want to provide parents
    # even if they are HD. That's why I need to remove the parents first
    
    # this finds the HD family members (effectively offspring) and adds them to 'ref_inds'
    # or removes them from the library phase
    #
    # NOTE: this, at this position, is not needed. You could also manipulate with
    # para_tab$provide_HD_inds_as_ref[i]. However, for consistency sake, I insert
    # it here to be the same with other positions in this script
    if(para_tab$provide_HD_fam_members_as_ref[i]){
      ref_inds <- unique(c(ref_inds, ped_fam_o[[1]][ped_fam_o[[1]] %in% HD_inds]))
    } else{
      ref_phase <- extract_individuals_phase(phase = ref_phase, 
                                            ID_inds = ped_fam_o[[1]][ped_fam_o[[1]] %in% HD_inds], 
                                            negate = TRUE)
    }
    
    if(para_tab$provide_HD_inds_as_ref[i]){
      ref_inds <- unique(c(ref_inds, unlist(ref_phase[ref_phase[[1]] %in% HD_inds,1])))
      ref_phase <- extract_individuals_phase(phase = ref_phase, ID_inds = ref_inds)
    } else{
      ref_phase <- extract_individuals_phase(phase = ref_phase, ID_inds = ref_inds)
    }
    
    # here, HD individuals are removed, if desired
    if(!para_tab$HD_inds_in_gt_set[i]){
      ls_masked_mt[[k]] <- extract_individuals_marker_table(marker_table = ls_masked_mt[[k]], 
                                                            ID_inds = HD_inds, negate = T)
    }
    # here it is imputed
    IN <- impute_familywise_BEAGLE5.3(arguments = arguments_imputing
                                   , marker_table = ls_masked_mt[[k]]
                                   , pedigree = ped_fam
                                   , memory_limit_GB = memory_limit_GB
                                   , mapo = mapo
                                   , verbose = FALSE
                                   , reference_phase = ref_phase)
    run_time <- run_time + IN$run_time
    
    
    imputed_mt <- convert_phase_to_marker_table(IN$out_phase[[1]])  
    
    # in the end I only want to end up with one marker table to analyze
    # here I add the results from every round of imputation
    # of course parents will be imputed and/or phased a couple of times
    # but I won't consider them later in the analysis
    if(is.null(runtime)){
      out_mt <- imputed_mt
      runtime <- run_time
    } else{
      out_mt <- cbind(out_mt, imputed_mt[,-1])
      runtime <- runtime + run_time
    }
    cat(
      paste0('\tFinished family ',k,' of ',length(ls_masked_mt),' families in trial ',i,'.\n')
    )
  }
} else{
  # the below is for the case that no familywise phasing shall be performed
  ########
  # here it is controlled how many individuals are provided for phasing
  if(str_detect(para_tab$offspring_type[i], 
                pattern = 'OP') & str_detect(para_tab$offspring_type[i], 
                                             pattern = 'perfam', negate = TRUE)){
    
    o_id <- sample(ped[[1]][31:nrow(ped)], replace = FALSE, size = ceiling(
      para_tab$fraction_offspring_provided_for_phasing[i]*length(ped[[1]][31:nrow(ped)])
    ))
    o_id <- o_id[match(ped[[1]][31:nrow(ped)], o_id)]
    o_id <- o_id[!is.na(o_id)]
    
    ped_p <- ped[1:30,]
    rand_p <- sample(ped_p[[1]], replace = FALSE,
                     size = ceiling(
                       para_tab$fraction_parents_provided_for_phasing[i]*length(ped_p[[1]])
                     ))
    rand_p <- rand_p[match(ped_p[[1]], rand_p)]
    rand_p <- rand_p[!is.na(rand_p)]
    
    ped_phasing <- extract_individuals_pedigree(pedigree = ped, 
                                                ID_names = unique(c(rand_p, o_id)))
    masked_mt_for_phasing <- extract_individuals_marker_table(marker_table = masked_mt, 
                                                              ID_inds = ped_phasing[[1]])
    
  } else{
    o_id <- extract_random_HD_offspring_per_fam(pedigree = ped, 
                                                fraction_HD_offspring_per_family = para_tab$fraction_offspring_provided_for_phasing[i])
    ped_o <- extract_offspring_pedigree(pedigree = ped)
    o_id <- o_id[match(ped_o[[1]], o_id)]
    o_id <- o_id[!is.na(o_id)]
    
    ped_p <- extract_parents_pedigree(pedigree = ped)
    rand_p <- sample(ped_p[[1]], replace = FALSE,
                     size = ceiling(
                       para_tab$fraction_parents_provided_for_phasing[i]*length(ped_p[[1]])
                     ))
    rand_p <- rand_p[match(ped_p[[1]], rand_p)]
    rand_p <- rand_p[!is.na(rand_p)]
    
    ped_phasing <- extract_individuals_pedigree(pedigree = ped, 
                                                ID_names = unique(c(rand_p, o_id)))
    masked_mt_for_phasing <- extract_individuals_marker_table(marker_table = masked_mt, 
                                                              ID_inds = ped_phasing[[1]])
  }
  ########
  
  # this is included to keep the code the same
  # ped_fam for imputation is freshly created later in the code and thus
  # won't be affected by this here
  ped_fam <- ped_phasing
  ped_fam[[2]] <- ped_phasing[[2]][[1]]
  # so the below is for the case that a phasing step shall be included
  # (just not familywise phasing)
  runtime <- NULL
  if(para_tab$BEAGLE_ref_phased_imputing[i]) {
    # here, it is just phased
    IN <- impute_familywise_BEAGLE5.3(arguments = arguments_phasing
                                   , marker_table = masked_mt_for_phasing
                                   , pedigree = ped_fam
                                   , memory_limit_GB = memory_limit_GB
                                   , mapo = mapo
                                   , verbose = TRUE)
    ref_phase <- IN$out_phase[[1]]
    run_time <- IN$run_time
    cat(paste0('\tFinished producing a reference phase in trial ',i,'.\n'))
  }
  
  ###########################################################
  if(para_tab$per_family_imputation_imputing[i]) {
    ls_masked_mt <- 
      extract_family_marker_table(pedigree = ped, marker_table = masked_mt)
    
    # now imputation is done for every family individually
    for(w in 1:length(ls_masked_mt)){
      ped_fam <- extract_individuals_pedigree(pedigree = ped, 
                                              ID_names = colnames(ls_masked_mt[[w]]))
      # the 'parents' population is confusing
      # they are not related so API cannot impute them, still I don't want to 
      # filter them out by name but use a general rule
      
      # here I check which individuals are parents to others
      ped_fam_p <- extract_parents_pedigree(pedigree = ped_fam)
      # who is an offspring of these parents?
      ped_fam_o <- extract_offspring_pedigree(pedigree = ped_fam, 
                                              parents_ID = ped_fam_p[[1]])
      # combine them
      ped_fam <- rbind(ped_fam_p,ped_fam_o)
      # because an individual could be parent AND offspring, I need to use 'unique'
      # (not sure if API actually would accept such individuals)
      ped_fam <- unique(ped_fam)
      
      # this bit is important: if less than 3 individuals are present in the pedigree,
      # so not even a trio, skip this round of the for loop
      # of course it is possible that only one parent is known and API wouldn't 
      # accept it but it works for the structure we have
      if(nrow(ped_fam) <= 2) {
        cat(paste0('Skipped family ', w,' in trial ',i,'.\n'))
        next
      }
      # this is to set the whole pedigree to the same family
      ped_fam[[2]] <- ped_fam[[2]][1]
      
      # maybe this is unnecessary cause impute_familywise_BEAGLE will handle it
      # according to the arguments provided but this option is safe
      if(para_tab$BEAGLE_ref_phased_imputing[i]){
        ref_phase_fam <- ref_phase
        
        # here I am checking which individuals should be given as reference
        ref_inds <- c()
        if(para_tab$provide_parents_as_ref[i]){
          ref_inds <- ped_fam_p[[1]]
        } else{
          ref_phase_fam <- extract_individuals_phase(phase = ref_phase_fam, 
                                                     ID_inds = ped_fam_p[[1]], 
                                                     negate = TRUE)
        }
        
        # this finds the HD family members (effectively offspring) and adds them to 'ref_inds'
        # or removes them from the library phase
        if(para_tab$provide_HD_fam_members_as_ref[i]){
          ref_inds <- unique(c(ref_inds, ped_fam_o[[1]][ped_fam_o[[1]] %in% HD_inds]))
        } else{
          ref_phase_fam <- extract_individuals_phase(phase = ref_phase_fam, 
                                                ID_inds = ped_fam_o[[1]][ped_fam_o[[1]] %in% HD_inds], 
                                                negate = TRUE)
        }
        
        if(para_tab$provide_HD_inds_as_ref[i]){
          ref_inds <- unique(c(ref_inds, unlist(ref_phase_fam[ref_phase_fam[[1]] %in% HD_inds,1])))
          ref_phase_fam <- extract_individuals_phase(phase = ref_phase_fam, ID_inds = ref_inds)
        } else{
          ref_phase_fam <- extract_individuals_phase(phase = ref_phase_fam, ID_inds = ref_inds)
        }
        
        # here, HD individuals are removed, if desired
        if(!para_tab$HD_inds_in_gt_set[i]){
          ls_masked_mt[[w]] <- extract_individuals_marker_table(marker_table = ls_masked_mt[[w]], 
                                                                ID_inds = HD_inds, negate = T)
        }
        
        IN <- impute_familywise_BEAGLE5.3(arguments = arguments_imputing
                                       , marker_table = ls_masked_mt[[w]]
                                       , pedigree = ped_fam
                                       , memory_limit_GB = memory_limit_GB
                                       , mapo = mapo
                                       , verbose = FALSE
                                       , reference_phase = ref_phase_fam)
        
        imputed_mt <- convert_phase_to_marker_table(IN$out_phase[[1]])  
        
        if(is.null(runtime)){
          out_mt <- imputed_mt
          runtime <- run_time + IN$run_time
        }  else{
          out_mt <- cbind(out_mt, imputed_mt[,-1])
          runtime <- runtime + IN$run_time
        }
      } else{
        IN <- impute_familywise_BEAGLE5.3(arguments = arguments_imputing
                                       , marker_table = ls_masked_mt[[w]]
                                       , pedigree = ped_fam
                                       , memory_limit_GB = memory_limit_GB
                                       , mapo = mapo
                                       , verbose = FALSE)
        
        imputed_mt <- convert_phase_to_marker_table(IN$out_phase[[1]])  
        
        if(is.null(runtime)){
          out_mt <- imputed_mt
          runtime <- IN$run_time
        } else{
          out_mt <- cbind(out_mt, imputed_mt[,-1])
          runtime <- runtime + IN$run_time
        }
      }
    }
  } else{
    # the below is for the case that imputation shall not be done familywise
    # but the reference phase from the whole population shall be used
    
    # I want to use the function 'impute_familywise_BEAGLE' even if I don't actually
    # impute familywise. I am using this function because it requires less 
    # preparation (by me) and does everything by itself.
    # If everyone is from the same family, 'impute_familywise_BEAGLE' is producing
    # the same result as 'impute_BEAGLE'
    #
    # so, I only need to set everyone to the same family, and then it works
    ped_fam <- ped
    ped_fam[[2]] <- ped_fam[[2]][[1]]
    
    if(para_tab$BEAGLE_ref_phased_imputing[i]){
      
      # here I am checking which individuals should be given as reference
      # again, be aware of the difference between ref_inds and ref_phase
      ped_fam_p <- extract_parents_pedigree(ped_fam)
      ped_fam_o <- extract_offspring_pedigree(ped_fam)
      # this is kinda hardcoded in here to enable extracting parents for the OP10n100 population
      if(str_detect(para_tab$offspring_type[i], 
                    pattern = 'OP') & str_detect(para_tab$offspring_type[i], 
                                                 pattern = 'perfam', negate = TRUE)){
        ped_fam_p <- ped_fam[1:30,]
        ped_fam_o <- ped_fam[31:nrow(ped_fam),]
      }
      
      ref_inds <- c()
      if(para_tab$provide_parents_as_ref[i]){
        ref_inds <- ped_fam_p[[1]]
      } else{
        # if parents should not provided as reference, this here is executed (removes the parents)
        ref_phase <- extract_individuals_phase(phase = ref_phase, 
                                               ID_inds = ped_fam_p[[1]], negate = TRUE)
      }
      
      # this finds the HD family members (effectively offspring) and adds them to 'ref_inds'
      # or removes them from the library phase
      if(para_tab$provide_HD_fam_members_as_ref[i]){
        ref_inds <- unique(c(ref_inds, ped_fam_o[[1]][ped_fam_o[[1]] %in% HD_inds]))
      } else{
        ref_phase <- extract_individuals_phase(phase = ref_phase, 
                                              ID_inds = ped_fam_o[[1]][ped_fam_o[[1]] %in% HD_inds], 
                                              negate = TRUE)
      }
      
      if(para_tab$provide_HD_inds_as_ref[i]){
        ref_inds <- unique(c(ref_inds, unlist(ref_phase[ref_phase[[1]] %in% HD_inds,1])))
        ref_phase <- extract_individuals_phase(phase = ref_phase, ID_inds = ref_inds)
      } else {
        ref_phase <- extract_individuals_phase(phase = ref_phase, ID_inds = ref_inds)
      }
      
      # impute_familywise_BEAGLE separates the families based on the family 
      # identifier, not based on the parents
      IN <- impute_familywise_BEAGLE5.3(arguments = arguments_imputing
                                     , marker_table = masked_mt
                                     , pedigree = ped_fam
                                     , memory_limit_GB = memory_limit_GB
                                     , mapo = mapo
                                     , verbose = FALSE
                                     , reference_phase = ref_phase)
      
      out_mt <- convert_phase_to_marker_table(IN$out_phase[[1]])
      runtime <- run_time + IN$run_time
    } else {
      
      # here, HD individuals are removed, if desired
      if(!para_tab$HD_inds_in_gt_set[i]){
        masked_mt <- extract_individuals_marker_table(marker_table = masked_mt, 
                                                      ID_inds = HD_inds, negate = T)
      }
      
      IN <- impute_familywise_BEAGLE5.3(arguments = arguments_imputing
                                     , marker_table = masked_mt
                                     , pedigree = ped_fam
                                     , memory_limit_GB = memory_limit_GB
                                     , mapo = mapo
                                     , verbose = FALSE)
      
      out_mt <- convert_phase_to_marker_table(IN$out_phase[[1]])
      runtime <- IN$run_time
    }
  }
}

################################################################################
# the marker table that was fed to the software is stored in mt_for_imputation
# this may be useful to have later
mt_for_imputation <- masked_mt_store
# later on, comparisons are made to the true marker table which is not necessarily
# the one that was provided to the software
mt <- mt_true

# out_mt should be stored in case it is needed somewhere
out_mt_store <- out_mt
################################################################################
########### analysis
performance_df$RuntimeOfSoftware[i] <- runtime
################################################################################

# in previous versions of this script, the parents were always excluded 
# from the analysis
# Now, I want this to depend on the input parameter table
# if none of the below conditions is right, the analysis is based on the whole 
# population. But I am not sure how that will work out if if the phasing step is
# done per family
if(para_tab$who_to_analyze[i] == 'only_HD_inds'){
  out_mt <- extract_individuals_marker_table(marker_table = out_mt, ID_inds = HD_inds)
} else if(para_tab$who_to_analyze[i] == 'only_LD_inds') {
  out_mt <- extract_individuals_marker_table(marker_table = out_mt, ID_inds = HD_inds, negate = TRUE)
} else if(para_tab$who_to_analyze[i] == 'non_parents_inds'){
  # plan is to extract all offspring 
  # and extract all parents and kick those out that are overlapping
  # check documentation of the functions of unsure how they are working
  offspring_IDs <- extract_offspring_pedigree(pedigree = ped_store)[[1]]
  offspring_IDs <- offspring_IDs[!(offspring_IDs %in% extract_parents_pedigree(pedigree = ped_store)[[1]])]
  out_mt <- extract_individuals_marker_table(marker_table = out_mt, ID_inds = offspring_IDs)
  
  # this will take care of the parents in case of open pollination
  # the second str_detect is negated
  if(str_detect(para_tab$offspring_type[i], 
                pattern = 'OP') & str_detect(para_tab$offspring_type[i], 
                                             pattern = 'perfam', negate = TRUE)){
    # this line is needed because out_mt doesn't contain
    # individuals cause the above lines are executed
    out_mt <- out_mt_store 
    out_mt <- extract_individuals_marker_table(marker_table = out_mt, 
                                               ID_inds = ped_store[[1]][1:30], negate = TRUE)
  }
} else if(para_tab$who_to_analyze[i] == 'only_parents_inds'){
  parents_id <- extract_parents_pedigree(pedigree = ped_store)[[1]]
  out_mt <- extract_individuals_marker_table(marker_table = out_mt, ID_inds = parents_id, negate = FALSE)
  
  if(str_detect(para_tab$offspring_type[i], 
                pattern = 'OP') & str_detect(para_tab$offspring_type[i], 
                                             pattern = 'perfam', negate = TRUE)){
    out_mt <- out_mt_store
    out_mt <- extract_individuals_marker_table(marker_table = out_mt, 
                                               ID_inds = ped_store[[1]][1:30], negate = FALSE)
  }
}

################################################################################
mt <- 
  extract_individuals_marker_table(marker_table = mt, 
                                   ID_inds = colnames(out_mt)[-1])
mt <- extract_marker_marker_table(marker_table = mt, 
                                  marker_IDs = out_mt[[1]])

out_mt <- order_marker_table_by_pedigree(marker_table = out_mt, 
                                         pedigree = colnames(mt)[-1])

masked_mt <- extract_individuals_marker_table(marker_table = masked_mt,
                                              ID_inds = colnames(out_mt)[-1])
masked_mt <- extract_marker_marker_table(marker_table = masked_mt, 
                                         marker_IDs = out_mt[[1]])
masked_mt <- order_marker_table_by_pedigree(marker_table = masked_mt, 
                                            pedigree = colnames(out_mt)[-1])

# this line is checking how many markers have been imputed
# it looks quite complicated but that because a software may leave some markers
# as missing after imputation so I can't just use the 
# number of markers (total) - number of markers (LD_map)
# Also, later we will experiment with simulating missing data so the number
# of loci/marker at which an imputation was done may exceed the number of 
# markers not on the LD map. So now, if in the whole marker table there is 
# at least one missing information at a locus and the software proposes a 
# solution (imputation result) this loci counts as imputed regardless of how many
# individuals had missing information
# vice versa, if there is one missing genotype at a locus left in the 
# output, the locus counts as not imputed
performance_df[i,37] <-  
  (sum(apply(masked_mt[,-1], MARGIN = 1, function(row){any(row == 9)}))
   - sum(apply(out_mt[,-1], MARGIN = 1, function(row){any(row == 9)})))

# 'parents_exclude' is chosen to be FALSE because the user can specify
# this in the input excel table
# Genotype analysis
anal_wp <- analyze_pop_imputation(true_marker_table = mt, 
                                  imputed_marker_table = out_mt, 
                                  parents_exclude = FALSE, 
                                  only_variable_loci = TRUE, 
                                  per_family = FALSE, 
                                  pedigree = ped)

# wp = whole population, pf = per family 
anal_pf <- analyze_pop_imputation(true_marker_table = mt, 
                                  imputed_marker_table = out_mt, 
                                  parents_exclude = FALSE, 
                                  only_variable_loci = TRUE, 
                                  per_family = TRUE,
                                  pedigree = ped)

performance_df[i,c(2:7, 36)] <- anal_wp[c(2,3,4,5,6,7,1)]
performance_df[i,8:13] <- anal_pf[c(2,3,4,5,6,7)]

# Marker analysis
# analysis for whole population at once
performance_df[i, c(14:21,24)] <- 
  analyze_marker_panel_imputation(true_marker_table = mt, 
                                  imputed_marker_table = out_mt)[c(2,3,4,5,6,7,8,9,1)]

# here I count how many loci are invariable
an <- analyze_marker_marker_table(mt)
performance_df[i,22] <- sum(an$MAF %in% 0)/nrow(an)

an <- analyze_marker_marker_table(out_mt)
performance_df[i,23] <- sum(an$MAF %in% 0)/nrow(an)

# calculates the number of fields with information given
# for the whole marker table (except HD individuals)
performance_df[i, 39] <- sum(mt[,-1] != 9) # for the true marker table
performance_df[i, 40] <- sum(masked_mt[,-1] != 9) # for the masked marker table
performance_df[i, 41] <- sum(out_mt[,-1] != 9) # for the imputed marker table

# Marker analysis for only variable markers
# only extract variable markers for analysis

# only extract variable markers for analysis
# important that this is done here because otherwise the parents will be excluded
# and then, by chance, different loci will be variable depending on the random draw
# of offspring
# it is a difficult question to solve: if the analysis should be based on only
# variable marker positions, from where do I get these positions?
# This was simple when I was supplying true data but with data that contains
# genoytping errors, I don't know where these errors are before the analysis
# also, the function below counts '9' as variation
# 
# I have decided to provide the true marker table to get the variable positions
# the interpretation of results based on variable positions is harder
# but this solution is the most meaningful for me

var_mt <- extract_variable_loci_marker_table(marker_table = mt)

# this step is included because it is possible that through NA filtering
#  some markers got lost
var_mt <- extract_marker_marker_table(marker_table = var_mt, 
                                      marker_IDs = out_mt[[1]])

var_mt <- extract_individuals_marker_table(marker_table = var_mt,
                                           ID_inds = colnames(out_mt)[-1])
var_mt <- order_marker_table_by_pedigree(marker_table = var_mt, 
                                         pedigree = colnames(out_mt)[-1])

var_out_mt <- extract_marker_marker_table(marker_table = out_mt, 
                                          marker_IDs = var_mt[[1]])
performance_df[i, c(25:32,35)] <- 
  analyze_marker_panel_imputation(true_marker_table = var_mt, 
                                  imputed_marker_table = var_out_mt)[c(2,3,4,5,6,7,8,9,1)]

# here I count how many loci are invariable
an <- analyze_marker_marker_table(var_mt)
performance_df[i,33] <- sum(an$MAF %in% 0)/nrow(an)

an <- analyze_marker_marker_table(var_out_mt)
performance_df[i,34] <- sum(an$MAF %in% 0)/nrow(an)


# calculates the number of fields with information given
# for only the variable positions of the marker table
performance_df[i, 42] <- sum(var_mt[,-1] != 9) # for the true marker table
var_masked_mt <- make_masked_marker_table(HD_marker_table = var_mt, 
                                          LD_map = LD_map, 
                                          HD_indiv = HD_inds, 
                                          set_NA_to = 9)
performance_df[i, 43] <- sum(var_masked_mt[,-1] != 9) # for the masked marker table
performance_df[i, 44] <- sum(var_out_mt[,-1] != 9) # for the imputed marker table

################################################################################
# the following part is not modified, but added to the advanced scripts like this one
# if you remember, 'mt_for_imputation' is the marker table contains the information
# that is provided to the software
performance_df[i, 45] <- sum(mt_for_imputation[,-1] != 9)

performance_df[i, 46] <- 
  mean(1-calc_nonmissing_marker_ind_marker_table(marker_table = mt_for_imputation)[[2]])

anal <- analyze_marker_marker_table(marker_table = mt_for_imputation)
performance_df[i, 47] <- mean(anal[[4]])
performance_df[i, 48] <- sum(anal[[3]] == 0)/nrow(anal)
performance_df[i, 49] <- mean(anal[[3]])

# this is also the number of markers after filtering
performance_df[i, 50] <- nrow(anal)

# the following is done for the very true data before any filtering
# or inducing errors
performance_df[i, 51] <- sum(mt[,-1] != 9)

performance_df[i, 52] <- 
  mean(1-calc_nonmissing_marker_ind_marker_table(marker_table = mt)[[2]])

anal <- analyze_marker_marker_table(marker_table = mt)
performance_df[i, 53] <- mean(anal[[4]])
performance_df[i, 54] <- sum(anal[[3]] == 0)/nrow(anal)
performance_df[i, 55] <- mean(anal[[3]])

# this is also the number of markers before filtering
performance_df[i, 56] <- nrow(anal)

##########
## n marker before and after filtering. Just to have it listed clearly
performance_df[i, 57] <- performance_df[i, 56]
performance_df[i, 58] <- performance_df[i, 50]
##########
# we only have genotype analyses for only the variable loci so far
# I want to include an analysis with all loci
# this is done in the following

anal_wp_all_loci <- analyze_pop_imputation(true_marker_table = mt, 
                                           imputed_marker_table = out_mt, 
                                           parents_exclude = FALSE, 
                                           only_variable_loci = FALSE, 
                                           per_family = FALSE, 
                                           pedigree = ped)

performance_df[i,c(59:64)] <- anal_wp_all_loci[c(2,3,4,5,6,7)]
##########

# and at last, I am recording the number of cores and amount of memory 
# that was allowing for this script8
performance_df[i, 65] <- memory_limit_GB
performance_df[i, 66] <- n_cores

################################################################################

cat(
  paste0('Finished trial ',i,' of ', nrow(para_tab),' trials.\n')
)
################################################################################
##############################################################################
# I don't want to wait until the whole script has finished to see some results
# Even though it takes some time, I will print the results after every run
# so that you have something in case the script crashes at a later run
setwd(odir)
# this line is included because I don't want to have a document for every run
# that's why I need to delete the previous one. This shall only be executed 
# after the first file was already produced (so starting with row 2)

out_table <- cbind(para_tab, performance_df)
time <- paste0("_",Sys.time())
time <- stringr::str_replace_all(time, pattern = '[-:]', "")
time <- stringr::str_replace_all(time, pattern = ' ', replacement = "_")
filename <- paste0(name_parameter_table,'_run',i,'_time-',time,'.txt')

# saving in a folder with the name of the parameter table
if(!dir.exists(name_parameter_table)){
  dir.create(name_parameter_table)
}

setwd(odir)
setwd(paste0('./',name_parameter_table))
data.table::fwrite(x = out_table, file = filename)

# saving the workspace
#saveRDS(object = out_table, file = paste0(name_parameter_table,'_run',i,'_time-',time,".RData"))
#save.image(file = paste0(name_parameter_table,'_run',i,'_time-',time,".RData"))
################################################################################
# this bit here is for merging the data so that I don't have to do it anymore
n_files_expected <- nrow(para_tab)

ls_files <- list.files()
txt_files <- stringr::str_detect(ls_files, pattern = '\\.txt$')
txt_files <- ls_files[txt_files]

# numbers <- 0 is included to keep the remaining part executable even if no
# file was found
numbers <- 0
numbers <- stringr::str_extract_all(txt_files, pattern = "run[[:digit:]]+")
numbers <- stringr::str_extract_all(unlist(numbers), pattern = '[[:digit:]]+', simplify = FALSE)
numbers <- as.numeric(unlist(numbers))

numbers <- unique(sort(numbers))

# the following should only become TRUE if all files are already written out
condition <- sum(numbers %in% 1:n_files_expected) == n_files_expected

if(condition){
  for(i in 1:n_files_expected){
    file <- txt_files[stringr::str_detect(txt_files, pattern = paste0('run',i))][1]
    IN <- data.table::fread(file)
    if(i == 1){out <- IN[i,]; next}
    
    out <- rbind(out, IN[i,])
  }
  # the same time as for the output above is used
  data.table::fwrite(out, file = paste0(name_parameter_table, '_all_time-',time,'.txt'))
  
  # here I am converting logical values to character
  ncol <- which(unlist(lapply(out, is.logical)))
  out <- as.data.frame(out)
  out[,ncol] <- lapply(out[,ncol], as.character)
  writexl::write_xlsx(out, paste0(name_parameter_table, '_all_time-',time,'.xlsx'))
}
################################################################################

################################################################################
# delete the directory in which this script was working
setwd(odir)
unlink(paste0('./',new_dir_name), recursive = TRUE)
### MODIFY_HERE: }