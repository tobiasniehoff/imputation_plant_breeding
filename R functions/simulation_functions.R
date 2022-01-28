#devtools::install_github("DominikMueller64/Meiosis", build_vignettes = TRUE)
library("Meiosis")
source("general_functions.R")

#options(digits = 7) # this is the default number of digits that will be printed
options(digits = 12) # this is needed (for visual inspection) for the make_unique_positions
# function. The elements are stored internally in R the same just, that if a lower
# digits number is chosen not all of them will be displayed


# convert functions

# this converts a mapo object to xoparam needed by the Meiosis package
convert_mapo_to_xoparam <- function(mapo) {
  n_chr <- unique(mapo[[2]])
  L <- unlist(lapply(n_chr, function(i) {max(mapo[mapo[[2]] == i, 3])}))
  xoparam <- Meiosis::create_xoparam(L)
  return(xoparam)
}
#xoparam <- convert_mapo_to_xoparam(map)


# this function converts a phase object to a pop object as it is produced by
# functions of the Meiosis package#
# use_names can be set to FALSE or TRUE according to whether one like sthe genotypes
# in the pop object to be named or not
convert_phase_to_meiosis_pop <- function(phase, mapo, use_names = FALSE) {
  pos <- extract_positions_mapo(mapo)
  n_loci <- unlist(lapply(pos, length))
  
  pos_hap1 <- seq(from = 1, to = nrow(phase), by = 2)
  pop <- list()
  
  # this function extracts the markers of each chromosomes and puts them in a list per chromosome
  #vector<- c("a", "b", "c", "d", "e", "f", "g", "h") # this is just an example
  list_markers <- function(vector, n_loci, marker_names = NULL) {
    l <- NULL
    start <- 0
      for(i in 1:length(n_loci)) {
        l <- c(l, list(vector[(start+1):(start + n_loci[i])]))
        start <- start + n_loci[i]
      }
    
    return(l)
  }
  #list_markers(vector, n_loci = c(2,3,2,1))
  
  for(i in pos_hap1) {
    hap1 <- list_markers(vector = as.integer(phase[i, -1]), n_loci = n_loci)
    hap2 <- list_markers(vector = as.integer(phase[i+1, -1]), n_loci = n_loci)
    ind <- list(paternal = hap1, maternal = hap2)
    pop <- list.append(pop, ind)
  }
  if(use_names) {names(pop) <- phase[[1]][pos_hap1]}
  return(pop)
}
# p <- convert_phase_to_meiosis_pop(phase, mapo = map)

# converts a pop object as produced my the meiosis functions into a phase object
convert_meiosis_pop_to_phase <- function(meiosis_pop, marker_names = NULL) {
  out <- lapply(meiosis_pop, function(ind) {
    lapply(ind, unlist)
  })
  out <- list.flatten(out)
  out <- list.rbind(out)
  pos_names <- sort(rep(seq(1,length(meiosis_pop),1), 2))
  if(!is.null(names(meiosis_pop))) {
    namse <- names(meiosis_pop)
  } else {namse <- paste0("ID_", seq(1, length(meiosis_pop)))}
  out <- data.table(Genotype = namse[pos_names], out)
  if(!is.null(marker_names)) {colnames(out)[-1] <- marker_names}
  return(out)
}
#phase <- convert_meiosis_pop_to_phase(meiosis_pop = pop)
#phase <- convert_meiosis_pop_to_phase(meiosis_pop = pop, marker_names = mapo[[1]])
# OR: #colnames(phase)[-1] <- mapo[[1]]

# wrapper function. To make it easier while programming
convert_meiosis_pop_to_marker_table <- function(meiosis_pop, marker_names = NULL){
  mt <- convert_phase_to_marker_table(
    convert_meiosis_pop_to_phase(meiosis_pop = meiosis_pop, marker_names = marker_names)
    )
  return(mt)
}


##################################
# make functions

# this function checks if any position occurs more than once and if so, increases
# the other identical positions by a tiny amount to make them acceptable by the Meiosis package
# the small error that is induced here is thought to not cause huge trouble later
# this is a recursive function
#
# only works if input is a vector
#
make_unique_positions <- function(positions, offset = 1/1000000000) {
  p <- positions
  # the second position that is identical will be increased by one million of a centi Morgan
  offset <- offset
  # to compare if any two positions are the same I verschiebe den Vektor um eins. 
  # Wenn der Nachbar gleich ist, müsste der Nachbar als TRUE widergegeben werden
  same <- c(p, -1) == c(-1, p)
  
  if(!any(same)) {return(p)}
  else {
    p[same] <- p[same] + offset
    p <- p[-length(p)]
    make_unique_positions(positions = p, offset = offset) # checks again(because maybe there are 3 identical positions)
  }
}
#p <- c(1,2,3,3,3,4,5,6,6)
#make_unique_positions(positions = p)


# this function takes an existing pedigree and makes a new one with the same
# pairings of parents but with desired number of offspring
make_simulated_pedigree <- function(pedigree, n_offspring) {
  a <- unique(ped[,3:4])
  d <- a[rowSums(a == 0) == 0,]
  
  df <- cbind(seq(1,nrow(d)*n_offspring)
              ,rep(seq(1, nrow(d)), each = n_offspring)
              ,rep(d[[1]], each = n_offspring)
              ,rep(d[[2]], each = n_offspring))
  df <- as.data.frame(df)
  df[stringr::str_length(df[,2]) == 1,2] <- paste0("0",df[str_length(df[,2]) == 1,2])
  df[,2] <- paste0("Fam", df[,2])
  df[,1] <- paste0("ID", df[,1], "_sim")
  parents <- extract_individuals_pedigree(pedigree = pedigree
                                          , ID_names = unique(c(pedigree[[3]], pedigree[[4]])))
  colnames(df) <- colnames(pedigree)
  df <- rbind(parents, df)
  
  return(as.data.table(df)) 
}
#make_simulated_pedigree(pedigree, n_offspring = 2)


# the below function is custom made for the imputation scripts. If you provide
# a pedigree, it is looked up which parents were crossed for a specific family
# and these crosses are done again (simulated). But the type of population for 
# parents and the offspring can be specified. For parent_type, 'real', 'DH' and 
# 'F[number]' are available. 'real' extracts the real parents from the provided
# phase object and uses them. If any of the other options is chosen, parents are
# simulated. For that, F1 individuals have to be generated. The F1s are produced 
# with female gametic control, i.e. a real parent is taken and pollinated with a
# pollen cloud of all parents (selfing is possible). Female gametic control is 
# used because I do want F1/DH/F... but these should be as genetically similar
# to the real parents from a diversity perspective so that I don't have an excess 
# in drift.
# 
# For 'offspring_type', 'DH', 'F[number]', 'OP[number]perfamn[number]' and 
# 'OP[number]n[number]' are available. 'OP' means 'open pollination' and can 
# either be done with all individuals or within a bi-parental family.
# The number of offspring can be controlled with 'n_off.family'.
# The number after OP refers to the number of generations of random mating
# AFTER the F1. E.g. to produce an F2, it takes 2 generations (-> F1 -> F2). To 
# produce an OP2 it takes 3 generations (->F1 -> OP1 -> OP2).
# The number after 'n' refers to the population size (number of individuals)
# in any of the generations
# of open pollination, except for the last one. The number of individuals in the
# last generation of open pollination can be controlled 'with n_off.family', 
# even if open pollination is not done per family.
#
# Note that to get the same number of individuals for 'OP2n10' as for 
# 'OP2perfamn10', for the latter case you need to specify 'n_off.family = x' 
# but for the first case it has to be 'n_off.family = x*[number families]'.
# For our data set, the number of families is 36.
make_simulated_data <- function(offspring_type, parent_type, n_off.family, seed,
                                phase, pedigree, map){
  
  if(missing(seed)){seed <- sample(1:10000,1)}
  ##############################################################################
  pedigree_store <- pedigree
  set.seed(seed) ## Seed R's rng
  Meiosis::seed_rng(seed = seed) ## Seed rng used by Meiosis
  
  # this extracts the parents from the pedigree based on whether they occur as parents
  # NOT based on their population name
  parents <- extract_parents_pedigree(pedigree)[[1]]
  parents_phase <- extract_individuals_phase(phase = phase, ID_inds = parents)
  
  mapo <- map[map[[2]] != 99,]
  mapo <- extract_marker_mapo(mapo = mapo, marker_IDs = colnames(phase))
  
  # this recreates the pedigree 'pedigree' in the way that the same crosses are 
  # made as in 'pedigree'.
  sim_ped <- make_simulated_pedigree(pedigree, n_offspring = n_off.family)
  
  # only executed for OP populations. The following is then overwriting 
  # 'sim_ped'
  # only execute if OP
  if(str_detect(offspring_type, pattern = 'OP')){
    pedigree_store <- pedigree
    # only execute if non familywise
    if(str_detect(offspring_type, pattern = 'perfam', negate = TRUE)){
      # what this is doing is to set the family name of all indidivuals that are not
      # founders (here: have a father) to the same family name
      pedigree[(pedigree[[3]] != '0'), 2] <- unique(pedigree[[2]])[2]
      pedigree[,3:4] <- 0
    }
    # n_off differs from n_off.family in that n_off means 'number of individuals
    # in every round of random mating'. 'n_off.family' means the number of 
    # individuals used for imputation.
    # For e.g. an F3 population, these two numbers are identical. Otherwise,
    # two or more individuals in F3 are derived from the same F2 plant or
    # not all F2 plants contribute to the F3 generation
    n_off <- as.numeric(str_extract(offspring_type, pattern = "[:digit:]+$"))
    sim_ped <- make_simulated_pedigree(pedigree, n_offspring = n_off)
    pedigree <- pedigree_store
  }
  
  # positions and xoparam are somthing that is required by the 'Meiosis' package
  # to do Meiosis simulations
  positions <- extract_positions_mapo(mapo)
  xoparam <- convert_mapo_to_xoparam(mapo)
  
  # normally, you would use the function 'breed_inds_from_pedigree' but that 
  # function doesn't allow DH, F1, F2,... or different types of inbreeding in the 
  # offspring
  parents_pop <- convert_phase_to_meiosis_pop(phase = parents_phase, mapo = mapo, use_names = TRUE)
  ##############################################################################
  # this bit is for simulating the parents
  # 'l' is just a placeholder for 'list' with no deeper meaning
  # the way parents are simulated if they should be the real ones, is that 
  # I first cross them randomly but with female gametic control. I.e. imagine a 
  # pollen cloud to which all parents contributed and with this pollen cloud you 
  # pollinate all of the original parents. The point of this female gametic control
  # is to get DH, F1, F2 parents that are from a diversity standpoint somewhat similar
  # to the real ones
  
  # the condition is TRUE if NOT the real parents should be used
  if(!str_detect(parent_type, pattern = 'real')){
    l <- list()
    # this is just making the F1
    for(i in 1:length(parents)){
      father_pos <- sample(seq(1:length(parents)),1)
      ind <- Meiosis::cross_geno(father = parents_pop[[father_pos]], mother = parents_pop[[i]], 
                                 positions = positions, xoparam = xoparam)
      l <- c(l, list(ind)) 
    }
    names(l) <- names(parents_pop)
    parents_pop <- l
    
    # the following is for selfing the parents or making DH
    if(parent_type == "DH"){
      parents_pop <- breed_make_dh_pop(pop = parents_pop, positions = positions, xoparam = xoparam)
    } else if(str_detect(parent_type, pattern = 'F[:digit:]+')){
      # this here is executed if the parents should be selfed
      
      # finding out how often the parents should be selfed
      # -1 because a F2 is only selfed once
      # If parents should be F1, n_selfings is 0. Then, 'breed_self_pop' returns
      # the input population unchanged
      n_selfings <- (as.numeric(str_extract(parent_type, pattern = '[:digit:]+')) - 1)
      
      parents_pop <- breed_self_pop(pop = parents_pop, n_selfings = n_selfings, 
                                    positions = positions, xoparam = xoparam)
    } else{
      stop(paste0("The parent type you specified cannot be recognized./n
                     Please name either 'real', 'DH', or 'F[number]'."))
    }
  }
  parents_phase <- convert_meiosis_pop_to_phase(meiosis_pop = parents_pop, marker_names = mapo[[1]])
  
  # this line could also be deleted
  parents_pop <- convert_phase_to_meiosis_pop(phase = parents_phase, mapo = mapo, use_names = TRUE) 
  # finished simulating parents
  ################################################################################
  # simulating the offspring
  
  if(str_detect(offspring_type, pattern = "F[:digit:]+")){
    # -1 because F2 means 1 selfing
    n_selfings <- as.numeric(str_extract(offspring_type, pattern = '[:digit:]+'))-1
    sim_phase <- breed_inds_from_pedigree(phase = parents_phase, pedigree = sim_ped, 
                                          mapo = mapo, n_selfings = n_selfings)
  } else if(offspring_type == 'DH'){
    l <- list()
    sim_ped_store <- sim_ped
    # extracts only individuals that have both parents known (i.e. no parent
    # is denoted as 0)
    sim_ped <- sim_ped[!apply(MARGIN = 1, sim_ped[, 3:4] == 0, any),]
    for(i in 1:nrow(sim_ped)) {
      # this monstrosity for extracting the parents is like "check in the pedigree
      # which parent to take", and then 'extract this parent from the meiosis pop
      # object'
      ind <- Meiosis::cross_geno(father = parents_pop[[sim_ped[[3]][i]]]
                                 , mother = parents_pop[[sim_ped[[4]][i]]],
                                 , positions = positions, xoparam = xoparam)
      l <- c(l, list(ind))
    }
    names(l) <- sim_ped[[1]]
    
    l <- breed_make_dh_pop(pop = l, positions = positions, xoparam = xoparam)
    out_phase <- convert_meiosis_pop_to_phase(l, marker_names = mapo[[1]])
    out_phase <- extract_individuals_phase(out_phase, ID_inds = parents, negate = TRUE)
    out_phase <- rbind(parents_phase, out_phase)
    
    sim_phase <- out_phase
    sim_ped <- sim_ped_store
  } else if(
    str_detect(offspring_type, 
               pattern = 'OP') & str_detect(offspring_type, 
                                            pattern = 'perfam', negate = TRUE)){
    # the above (horribly written) condition is TRUE for random mating but not 
    # per family random mating
    
    # this extracts the number of generations of random mating
    # -1 is used here because the number of individuals from the one before the last
    # to the last generation may change (if someone is interested in it)
    n_generations <- as.numeric(
      str_extract(
        str_extract(offspring_type, pattern = 'OP[:digit:]+'), 
        pattern = '[:digit:]+'))-1
    
    # 'n_off' should here better be understood as 'population size'
    n_off <- as.numeric(str_extract(offspring_type, pattern = "[:digit:]+$"))
    
    sim_pop <- breed_random_mating_pop(pop = parents_pop, n_random_mating = n_generations, 
                                       positions = positions, xoparam = xoparam, 
                                       n_offspring = n_off, allow_selfing = TRUE)
    # the parameter 'n_off.family' here means total pop size as there are no families
    # anymore
    sim_pop <- breed_random_mating_pop(pop = sim_pop, n_random_mating = 1, 
                                       positions = positions, xoparam = xoparam, 
                                       n_offspring = n_off.family, allow_selfing = TRUE)
    
    out_phase <- convert_meiosis_pop_to_phase(meiosis_pop = sim_pop, marker_names = mapo[[1]])
    
    # the below is generating a new simulated pedigree because the one used above
    # for generating the F1s can be different than the one needed for imputation
    pedigree <- pedigree_store
    
    # these many lines below are for makeing a new simulated pedigree
    sim_ped <- make_simulated_pedigree(pedigree = pedigree, n_offspring = n_off.family)
    parents_ped <- extract_parents_pedigree(sim_ped)
    
    offspring_ped <- extract_individuals_pedigree(pedigree = sim_ped, 
                                                  ID_names = parents_ped[[1]], 
                                                  negate = TRUE)
    offspring_ped <- offspring_ped[1:n_off.family,]
    offspring_ped[[2]] <- offspring_ped[[2]][1]
    
    sim_ped <- rbind(parents_ped, offspring_ped)
    sim_ped[,3:4] <- 0
    
    sim_phase <- rbind(parents_phase, out_phase)
    # assigning the right names to the phase object. Rep 2 because every ind occurs
    # twice in a phase object
    sim_phase[[1]] <- rep(sim_ped[[1]], each = 2)
  } else if(str_detect(offspring_type, pattern = 'perfam')){
    # the above if condition is only TRUE if per family random mating shall be done 
    l <- list()
    # here I am simulating F1 individuals
    sim_phase <- breed_inds_from_pedigree(phase = parents_phase, pedigree = sim_ped, mapo = mapo
                                          , n_selfings = 0)
    
    # this extracts the number of generations of random mating
    # -1 is used here because the number of individuals from the one before the last
    # to the last generation may change (if someone is interested in it)
    n_generations <- as.numeric(
      str_extract(
        str_extract(offspring_type, pattern = 'OP[:digit:]+'), 
        pattern = '[:digit:]+'))-1
    
    # 'n_off' should here better be understood as 'population size'
    n_off <- as.numeric(str_extract(offspring_type, pattern = "[:digit:]+$"))
    # I do this '-1' because the first family are the parents which I don't need
    # to simulate later
    ls_fam_peds <- extract_family_pedigree(pedigree = sim_ped)[-1]
    
    out_phase <- c()
    for(i in 1:length(ls_fam_peds)){
      fam_ids <- ls_fam_peds[[i]][[1]]
      # sim_phase includes the parents but the parents IDs are not in fam_ids
      fam_phase <- extract_individuals_phase(phase = sim_phase, ID_inds = fam_ids)
      fam_pop <- convert_phase_to_meiosis_pop(phase = fam_phase, mapo = mapo, use_names = TRUE)
      fam_pop <- breed_random_mating_pop(pop = fam_pop, 
                                         n_random_mating = n_generations, ########### important
                                         positions = positions, 
                                         xoparam = xoparam, n_offspring = n_off)
      
      fam_pop <- breed_random_mating_pop(pop = fam_pop, positions = positions,
                                         xoparam = xoparam, n_random_mating = 1, 
                                         n_offspring = n_off.family)
      fam_phase <- convert_meiosis_pop_to_phase(meiosis_pop = fam_pop, marker_names = mapo[[1]])
      
      if(is.null(out_phase)){out_phase <- fam_phase} 
      else {out_phase <- rbind(out_phase, fam_phase)}
    }
    
    # sim_ped was already created above
    # but the sim_ped above can be different than the one here, say because you 
    # want to have per family random mating with population size 5 but the final analysis
    # demands 100 individuals per family.
    sim_ped <- make_simulated_pedigree(pedigree = pedigree_store, n_offspring = n_off.family)
    
    sim_phase <- rbind(parents_phase, out_phase)
    sim_phase[[1]] <- rep(sim_ped[[1]], each = 2)
  } else{
    stop (paste0('Could not recognize the offspring type. Please use "DH", "F3", 
  "OP10perfamn5" or "OP10n100" (just as an example)./n'))
  }
  
  if((nrow(sim_phase)/2) != nrow(sim_ped)) stop('sim_ped and sim_phase contain a
                                              different number of individuals')
  
  return(list(sim_phase = sim_phase, sim_ped = sim_ped, mapo = mapo, seed = seed))
}

# in the example below, the parents are simulated to be F1. The crosses of the
# families in 'pedigree' are recreated and in this case, the F1 offspring is 
# randomly mating for 4(!) generations with population size of 10 within each 
# family. The population size per family in the 5th generation however
# is 50.
# IN_data <- make_simulated_data(offspring_type = 'OP5perfamn10', parent_type = 'F1',
#                                n_off.family = 50, phase = phase, 
#                                pedigree = pedigree, map = map, seed = 1)
# str(IN_data)



################################################################################
# extract functions

# extracts the positions of a mapo object and returns the positions in a format
# needed by Meiosis functions. Also, if any two neighboring markers have the same 
#position (which makes them non acceptable by Meiosis functions) a tiny error 
#(1/1000000000) is added to make the positions unique
extract_positions_mapo <- function(mapo) {
  # orders the mapo object first by positions and then by chromosome. Needed
  # in case the mapo object is messed up in which the checking of neighbouring 
  # positions won't work
  mapo <- as.data.frame(mapo)
  mapo <- order_mapo_by_position(mapo) 
  n_chr <- unique(mapo[[2]])
  # extracts the positions per chromosome and puts them in a separate list
  positions <- lapply(n_chr, function(i) {mapo[mapo[[2]] == i, 3]}) 
  
  # checks if neighboring positions are identical and if so adds a small offset
  positions <- lapply(positions, make_unique_positions)
  Meiosis::check_positions(positions) # a Meiosis function
  return(positions)
}
#extract_positions_mapo(map)
# not possible to leave values as missing, i.e. 9 for simulation

# breed functions
# all functions with the prefix breed refer to simulation breeding activities


# creates a number of genotypes from the cross between ind1 and ind2
# if ind2 is not specified, ind1 is selfed
breed_cross_inds <- function(ind1, ind2 = NULL, n_offspring = 1,
                             positions= positions, xoparam = xoparam) {
  if(class(ind2) == "NULL") {ind2 <- ind1}
  Meiosis::check_geno_individual(ind1)
  Meiosis::check_geno_individual(ind2)
  Meiosis::check_positions(positions)
  
  pop <- replicate(n_offspring, Meiosis::cross_geno(ind1, ind2, positions, xoparam), simplify = FALSE)
  
  return(pop)
}
#breed_cross_inds(ind1 = pop[[3]], ind2 = pop[[10]], n_offspring = 10, positions = positions, xoparam = xoparam)

# this function selfs all individuals in a population (default = 1)
breed_self_pop <- function(pop, n_selfings = 1
                           , positions = positions, xoparam = xoparam) {
  Meiosis::check_positions(positions)
  Meiosis::check_geno_individual(pop[[1]])
  # this if is kinda dumb but maybe important later when used automatically
  if(n_selfings > 0){
    for(s in 1:n_selfings) {
      for(i in 1:length(pop)) {
        pop[[i]] <- Meiosis::self_geno(individual = pop[[i]], positions = positions, 
                                       xoparam = xoparam)
      }
    }
  }
  return(pop)
}
#breed_self_pop1(pop, n_selfings = 5, positions = positions, xoparam = xoparam)

# takes a 
breed_make_dh_pop <- function(pop, positions = positions, xoparam = xoparam){
  Meiosis::check_geno_individual(pop[[1]])
  Meiosis::check_positions(positions)
  return(lapply(pop, function(i) {Meiosis::dh_geno(i, positions = positions, xoparam = xoparam)}))
}
# DH_pop <- breed_make_dh_pop(pop, positions = positions, xoparam = xoparam)
# DH_pop[[1]][[1]][[1]] == DH_pop[[1]][[2]][[1]] # this compares the maternal and paternal chromosome


# takes a pop object and allows the population to mate at random (so selfing is possible)
# the number of generations (n_random_mating) (default = 1) can be set

# Update 28.01.2021: It is possible now to breed a different number of individuals
# than the number of individuals present in the provided pop object
# Also, it is possible to control selfing. Parents are drawn randomly from pop.
# With 'allow_selfing' you can control if this could be the same parent by chance
# The newly built in functionalities should produce the same results in older scripts
# that are using this function
breed_random_mating_pop <- function(pop, n_random_mating = 1,
                                    positions = positions, xoparam = xoparam, 
                                    n_offspring, allow_selfing = TRUE){
  new_pop <- pop
  if(missing(n_offspring)) {n_offspring <- length(pop)}
  
  for(i in 1:n_random_mating) {
    new_pop <- lapply(seq(1,n_offspring), function(i) {
      parents <- sample(pop, size = 2, replace = allow_selfing)
      Meiosis::cross_geno(father = parents[[1]], 
                          mother = parents[[2]], 
                          positions = positions, 
                          xoparam = xoparam)
    })
    pop <- new_pop
  }
  return(pop)
}
#rndmat_pop <- breed_random_mating_pop(pop, n_random_mating = 10, positions = positions, xoparam = xoparam)


# this function takes a pedigree object and a phase object
# it will simulate all individuals for which their parents are known (so none 
# of them is 0). It will look up the parents from the supplied phase object
# it will return a phase object with the parents and the simulated offspring
# combined
#
# say an individual has both parents known but is also a parent of other
# individuals, then this individual will not be simulated but it's offspring
# will be simulated based on its parents real data
breed_inds_from_pedigree <- function(phase, pedigree, mapo, n_selfings = 0){
  parents <- extract_parents_pedigree(pedigree)[[1]]
  phase <- extract_individuals_phase(phase = phase, ID_inds = parents)
  
  pop <- convert_phase_to_meiosis_pop(phase = phase, mapo = mapo, use_names = TRUE)
  
  positions <- extract_positions_mapo(mapo)
  xoparam <- convert_mapo_to_xoparam(mapo)
  
  l <- list()
  # extracts only individuals that have both parents known (i.e. no parent
  # is denoted as 0)
  pedigree <- pedigree[!apply(MARGIN = 1, pedigree[, 3:4] == 0, any),]
  for(i in 1:nrow(pedigree)) {
    # this montrosity for extracting the parents is like "check in the pedigree
    # wh9ich parent to take", and then 'extract this parent from the meiosis pop
    # object'
    ind <- Meiosis::cross_geno(father = pop[[pedigree[[3]][i]]]
                               , mother = pop[[pedigree[[4]][i]]],
                               , positions = positions, xoparam = xoparam)
    l <- c(l, list(ind))
  }
  names(l) <- pedigree[[1]]
  l <- breed_self_pop(pop = l, n_selfings = n_selfings, positions = positions, xoparam = xoparam)
  out_phase <- convert_meiosis_pop_to_phase(l, marker_names = mapo[[1]])
  out_phase <- extract_individuals_phase(out_phase, ID_inds = parents, negate = TRUE)
  out_phase <- rbind(phase, out_phase)
  return(out_phase)
}
# breed_inds_from_pedigree(phase = IN, pedigree = sim_ped, mapo = mapo)
#
# the following is showing how the function works if an individual is offspring 
# and parent. In this example, G1 is a parent from P16 and P34 from the real dataset
# sp <- pedigree[c(16,34),]
# sp <- rbind(sp,list('G1', "f1", "P16", "P34"))
# sp <- rbind(sp, list('G2', "f2", "G1", "P16"))
# ph <- extract_individuals_phase(IN, ID_inds = c("P16", "P34", "G1"))
# d <- breed_inds_from_pedigree(phase = ph, pedigree = sp, mapo = mapo)
# d[, 1:10]