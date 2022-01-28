library(data.table)
library(stringr)
library(R.utils)
library(rlist)
library(readxl)
library(writexl)

`%notin%` <- Negate(`%in%`)

# dt <- as.data.table(matrix(sample(c(0,1,1,2), replace = T, size = 16), nrow = 4))
# dt <- cbind(c("ID1", "ID1", "ID2", "ID2"),dt)
# colnames(dt) <- c("Genotype", "C1M1","C1M2","C1M3","C2M1")
# phase <- dt

###### order functions
# orders the individuals in a marker_table by their occurrence in a pedigree object
# it is also possible to just supply a character vector. Then only the individuals 
# in the vector are used for ordering. All unmentioned individuals are added in the order
# of their occurrence after the mentioned ones
#
# this function was changed on 06.01.2021 and is now more than 4 times faster
order_marker_table_by_pedigree <- function(pedigree, marker_table) {
  if(any(class(pedigree) == "character")) {pedigree <- list(pedigree)}
  
  # extracting the individuals that are present in the marker table from the
  # pedigree (chopping off all individuals not mentioned in the marker table)
  pedigree <- pedigree[[1]][pedigree[[1]] %in% colnames(marker_table[,-1])]
  
  # getting the positions of individuals that were found in the marker table
  pos <- match(pedigree, colnames(marker_table)) 
  
  # adding 1 for the marker name column and all the individuals that were not 
  # mentioned by the pedigree
  pos <- c(1,pos,which(!colnames(marker_table) %in% pedigree)[-1])
  
  return(marker_table[, ..pos])
}

#orders the markers in a marker_table object by the occurence in a mapo object 
# (marker, chrom, pos_start). It is also possible to just give in the names of 
# the markers that should be ordered. Then the marker_table is ordered by them
# and all markers that were not specified but are still in the marker_table
# are added to the output in the order in which they occurred in the marker_table
order_marker_table_by_mapo <- function(mapo, marker_table) {
  if(any(class(mapo) == "character")) {mapo <- list(mapo)}
  pos <- lapply(mapo[[1]], function(i) {which(marker_table[[1]] == i)})
  out <- marker_table[unlist(pos),]
  out <- rbind(out, marker_table[marker_table[[1]] %notin% out[[1]],])
  
  return(out)
}

#order_marker_table_by_map(c("C1M1", "C1M2"), marker_table)[1:10, 1:10]
#order_marker_table_by_map(map, marker_table)[1:10, 1:10]

# orders the rows by their position in descending order within the chromosomes
# chromosomes are also ordered in descending order
order_mapo_by_position <- function(mapo) {
  # the conversion to numeric is necessary because you can do huge mistakes if 
  # the positions or chromosomes are saved as character by coincidence
  mapo[[2]] <- as.numeric(mapo[[2]])
  mapo[[3]] <- as.numeric(mapo[[3]])
  
  mapo <- mapo[order(mapo[[3]]),]
  mapo <- mapo[order(mapo[[2]]),]
  return(mapo)
}


# orders the markers (columns) of a phase object according to their occurence in a mapo
# object or a vector with marker names. All the markers that are not mentioned in the 
# vector or the map but are in the phase object are added to the ordered object in the 
# end
order_phase_by_mapo <- function(mapo, phase) {
  if(any(class(mapo) == "character")) {mapo <- list(mapo)}
  pos <- lapply(mapo[[1]], function(i) {which(colnames(phase) == i)})
  out <- subset(phase, select = unlist(pos))
  out <- data.table(Genotype = unname(phase[,1]), 
                    out, 
                    subset(phase, select = colnames(phase)[-1] %notin% colnames(out))
  )
  return(out)
}

# this function orders a phase object by a provided pedigree.
# It is possible to specify a character vector containing the names of individuals
# instead of a pedigree object. If the phase object contains more individuals
# than the pedigree, phase is ordered by pedigree and the remaining individuals
# are added to the ordered phase object in order of their occurrence in the 
# provided phase object. It works similar as order_marker_table_by_pedigree
#
# also works for GenotypeFile objects (output of API and API2)
order_phase_by_pedigree <- function(phase, pedigree) {
  if(any(class(pedigree) == "character")) {pedigree <- list(pedigree)} 
  # because it is also allowed to supply a character vector containing the names
  pos <- unlist(lapply(pedigree[[1]], function(i) {which(phase[[1]] == i)}))
  out <- phase[pos, ]
  
  out <- rbind(out, phase[phase[[1]] %notin% pedigree[[1]],])
  return(out)
}
#order_phase_by_pedigree(phase, pedigree)
#order_phase_by_pedigree(phase, c("P1", "G12", "P2"))
#GenotypeFile <- fread("ImputedDescendantGenotypes.txt")
#order_phase_by_pedigree(GenotypeFile, pedigree[[1]])

# works same as order_phase_by_pedigree,m just or geno objects
order_geno_by_pedigree <- function(geno, pedigree) {
  if(any(class(pedigree) == "character")) {pedigree <- list(pedigree)} 
  pos <- unlist(lapply(pedigree[[1]], function(i) {which(rownames(geno) == i)}))
  out <- geno[pos, ]
  
  out <- rbind(out, geno[rownames(geno) %notin% pedigree[[1]],])
  return(out)
}
#order_geno_by_pedigree(geno, pedigree)

# check functions

# this function checks if a pedigree is resolvable, i.e. if all individuals
# listed as parents can be traced back to founders or are founders.
# It first checks if all individuals listed as parents are actually also
# listed as individuals.
# Secondly, it check if there are loops in the pedigree.
# See the example below.
# This function will always return the same format:
# a list for indicating whether the pedigree is resolvable or not (TRUE/FALSE),
# then showing the names of parents that are missing and individuals that are affected
# and third the individuals that are involved in a loop
check_resolvable_pedigree <- function(pedigree){
  # first check if any listed parents are missing as individuals
  old <- pedigree
  l_p <- unique(unlist(old[,3:4]))
  l_p <- l_p[!(l_p %in% 0)]
  # gives a vector of names of individuals with missing parents
  missing_parent <- l_p[!(l_p %in% unique(old[[1]]))]
  # gives a pedigree object involving the individuals with problematic parents
  wrong_ind <- old[old[[3]] %in% missing_parent | old[[4]] %in% missing_parent,]
  
  # this will check if there are any loops in the pedigree
  continue <- TRUE
  while(continue) {
    new <- extract_parents_pedigree(old)
    continue <- !identical(new, old)
    if(continue) {old <- new}
  }
  loop_ped <- new
  
  out <- list(Resolvalbe = FALSE, 
              MissingParent = list(MissParentID = missing_parent, OffMissParent = wrong_ind),
              IndsLoop = loop_ped)
  
  if(length(missing_parent) >0){
    warning('Some individuals have parents listed that do not occur as individuals
            in the pedigree.')
  }
  
  if(nrow(loop_ped) > 0){
    warning("Not every individual can be traced back to founders. Maybe loops are in the pedigree.
            See individuals that cannot be traced back in the output")
  } 
  
  if(length(missing_parent)  == 0 & nrow(loop_ped) == 0) {
    cat("Good. The pedigree seems fine.\n")
    out <- list(Resolvalbe = TRUE, 
                MissingParent = list(MissParentID = missing_parent, OffMissParent = wrong_ind),
                IndsLoop = loop_ped)
  }
  
  return(out)
}
# p <- pedigree[1:13,]
# p[3,c(1,3,4)] <- list('O1',"P1","P2")
# p[4,c(1,3,4)] <- list('O2',"P1","P1")
# p[5,c(1,3,4)] <- list('P3',0,0)
# p[6,c(1,3,4)] <- list('O3',"O1","O2")
# p[7,c(1,3,4)] <- list('O4',"O2",0)
# p[8,c(1,3,4)] <- list('O5',"Miss1","O2")
# p[9,c(1,3,4)] <- list('O6',"O5","P2")
# p[10,c(1,3,4)] <- list('O7','O8',0)
# p[11,c(1,3,4)] <- list('O8','O7',0)
# p[12,c(1,3,4)] <- list('O9','O10',0)
# p[13,c(1,3,4)] <- list('O10','O9',0)
# check_resolvable_pedigree(p)
# check_resolvable_pedigree(pedigree)


###### convert functions

# makes the genotype object needed for imputation with individuals on rows and markers on columns
# also works on data.table
convert_marker_table_to_GenotypeFile <- function(marker_table) {
  transposed_mt <- transpose(marker_table[,-1])
  colnames(transposed_mt) <- unlist(marker_table[,1])
  transposed_mt <- cbind(Genotypes = colnames(marker_table)[-1], transposed_mt)
  return(transposed_mt)
}

# takes a pedigree object (Genotype ID, Population ID, Father ID, Mother ID) and
# returns a PedigreeObject as it is required as input for API per chromosome imputation
convert_pedigree_to_PedigreeFile <- function(pedigree) {
  if(ncol(pedigree) != 4) {
    stop("The pedigree input is expected to have 4 columns (ID, pop, father, mother)")
  } else {
    PedigreeFile <- pedigree[,c(1,3,4)]
  }
  return(PedigreeFile)
}

# takes a geno object (columns as SNP1_1 SNP1_2) and produces a marker_table object
# based on it
convert_geno_to_marker_table <- function(geno) {
  
  maternal <- str_detect(colnames(geno), pattern = "_1$")
  if(any(unique(maternal == c(T,F)) == FALSE))
  {stop("Something is odd about the order of columns.")}
  
  paternal <- str_detect(colnames(geno), pattern = "_2$")
  if(any(unique(paternal == c(F,T)) == FALSE))
  {stop("Something is odd about the order of columns.")}
  
  
  marker_names <- str_remove(colnames(geno), pattern = "_[1-2]$")
  if(any(marker_names[maternal] == marker_names[paternal]) == FALSE) {
    stop("Seems as if at least one marker haplotype is not next to his other haplotype.")
  }
  
  obj <- geno[, maternal] + geno[,paternal]
  obj[obj > 2] <- 9 # 1+1 = 2 -> fine, 1 + 9 = 10 >2 -> 9 in mt
  obj <- t(obj)
  obj <- apply(obj, MARGIN = 2, as.integer)
  obj <- data.table(marker = marker_names[maternal], obj)
  
  return(obj)
}

# converts from phase to geno objects. Phase object are objects as the yare read 
# like 
# fread("ImputedDescendantPhase.txt"). The geno object is a data.frame because 
# it needs rownames
convert_phase_to_geno <- function(phase_obj, marker_IDs = NULL) {
  cols <- (ncol(phase_obj)-1)*2
  rows <- nrow(phase_obj)/2
  m <- as.data.table(matrix(nrow = rows, ncol = cols))
  maternal <- seq(from = 1, to = cols, by = 2)
  paternal <- seq(from = 2, to = cols, by = 2)
  m[,maternal] <- phase_obj[seq(1, rows*2, 2), -1]
  m[,paternal] <- phase_obj[seq(2, rows*2, 2), -1]
  rownames(m) <- phase_obj[[1]][seq(1, rows*2, 2)]
  
  if(length(marker_IDs) == 0) {
    marker_IDs <- colnames(phase_obj)[-1]
  }
  
  if(ncol(m)/2 != length(marker_IDs)) {
    warning("Number of marker IDs does not match number of markers in phase object.")
  }
  
  colnames(m)[seq(from = 1, to = cols, by = 2)] <- paste(marker_IDs, 1, sep = "_")
  colnames(m)[seq(from = 2, to = cols, by = 2)] <- paste(marker_IDs, 2, sep = "_")
  
  geno <- as.data.frame(m)
  rownames(geno) <- rownames(m)
  return(geno)
}

# converts phase format to marker_table format
convert_phase_to_marker_table <- function(phase) {
  maternal <- seq(1, nrow(phase), 2)
  paternal <- maternal + 1
  
  ind_names <- phase[[1]][maternal]
  marker_names <- colnames(phase)[-1]
  phase <- phase[, -1]
  obj <- phase[maternal,] + phase[paternal,]
  obj[obj > 2] <- as.integer(9)
  obj <- transpose(obj)
  obj <- cbind(marker = marker_names, obj)
  colnames(obj)[-1] <- ind_names
  
  return(obj)
}
#convert_phase_to_marker_table(phase)

# convert_phase_to_geno(phase_desc)[1:10, 1:10]
# convert_phase_to_geno(phase_desc, marker_IDs = map[map[[2]] == 1, 1])[1:10, 1:10]
# if marker ID are procided like this. make sure the map object and the phase object 
# contain the same markers

# converts geno object to phase objects as required by API
convert_geno_to_phase <- function(geno) {
  maternal <- seq(1, ncol(geno), 2)
  paternal <- seq(2, ncol(geno), 2)
  
  haplo1 <- subset(geno, select = maternal)
  haplo2 <- subset(geno, select = paternal)
  
  phase <- rbind(as.matrix(haplo1), as.matrix(haplo2))
  n <- nrow(phase)/2
  pos <- seq(1, n, 1)
  pos <- sort(c(pos, pos))
  hap2 <- seq(2, length(pos), 2)
  pos[hap2] <- pos[hap2] + n
  
  namse <- rownames(phase[pos,])
  phase <- data.table(namse, apply(phase[pos,], MARGIN = 2, as.integer))
  # with pos, I reorder the row so that the right ones are underneath each other
  
  colnames(phase) <- c("Genotype",
                       stringr::str_extract(colnames(haplo1), pattern = "[:alnum:]+"))
  
  return(phase)
}


# Converts what you read in from the genotype file to marker table
convert_GenotypeFile_to_marker_table <- function(GenotypeFile) {
  mt <- transpose(GenotypeFile[,-1])
  mt <- cbind(marker = colnames(GenotypeFile)[-1], mt)
  colnames(mt)[-1] <- GenotypeFile[[1]]
  
  return(mt)
}
# in <- fread("ImputedDescendantGenotypes.txt")
# convert_GenotypeFile_to_marker_table(in)

# creates a map object in PLINK format
# PLINK wants to have the physical position in the fourth column
# if no fourth column is present in the mapo object, this function just numerates
# from 1 to number(marker/chrom) in that column
convert_mapo_to_PLINK_map <- function(mapo) {
  if(ncol(mapo) == 3){
    pos <- unlist(lapply(table(mapo[[2]]), function(i) {seq(1, i, 1)}))
  } else{
    if(ncol(mapo) != 4) {
      stop("Not sure what those columns are. Please provide
                              3 or 4 columns only.")
    }
    pos <- mapo[[4]]
  }
  map_object <- cbind(mapo[,c(2, 1, 3)], pos_bp = pos)
  
  return(map_object)
}
# pm <- convert_mapo_to_PLINK_map(mapo = map)

# converts a mapo object to a genetic map object that is required by API
# and impute_API
convert_mapo_to_GeneticMapFileAPI <- function(mapo) {
  return(mapo[,c(2,1,3)])
}
# convert_mapo_to_GeneticMapFileAPI(map)


# converts marker_table to geno object. Note that the geno object is NOT phased
# object needed for imputation with individuals on rows and markers on columns
# this function was written for API imputation as a whole chromosome at once
# Geno could then be converted to PLINK ped file. But since API splits up the 
# chromosomes internally anyway, we could also impute every chromosome and then
# merge the results to get the whole chromosome. I thought API whole genome
# imputation would allow to impute markers for which it it unknown on which 
# chromosome they are
convert_marker_table_to_geno <- function(marker_table) {
  transposed_mt <- as.matrix(t(marker_table[,-1]))
  
  hap_1 <- transposed_mt
  hap_2 <- transposed_mt
  
  pos <- hap_1 == 2 # assigned to pos cause it can be reused and saves time
  hap_1[pos] <- 1
  hap_2[hap_2 == 1] <- 0
  hap_2[pos] <- 1
  
  colnames(hap_1) <- paste(marker_table[[1]], 1, sep = "_")
  colnames(hap_2) <- paste(marker_table[[1]], 2, sep = "_")
  
  geno <- matrix(nrow = nrow(transposed_mt), ncol = ncol(transposed_mt)*2)
  
  geno[,seq(1,ncol(hap_1)*2, 2)] <- hap_1
  geno[,seq(2,ncol(hap_2)*2, 2)] <- hap_2
  
  geno <- as.data.frame(geno)
  colnames(geno)[seq(1,ncol(hap_1)*2, 2)] <- colnames(hap_1)
  colnames(geno)[seq(2,ncol(hap_2)*2, 2)] <- colnames(hap_2)
  
  rownames(geno) <- colnames(marker_table)[-1]
  
  return(geno)
}
#convert_marker_table_to_geno(mt)


####### extract functions

# A function to extract individuals specified in ID_names from a pedigree object
extract_individuals_pedigree <- function(pedigree, ID_names, negate = FALSE) {
  if(negate) {out <- pedigree[pedigree[[1]] %notin% ID_names,]}
  else{out <- pedigree[pedigree[[1]] %in% ID_names,]}
  
  return(out)
}

# A function that looks up all the parents in the father and mother column of a pedigree
# object and returns only those rows of the pedigree object that contain parents information
# if no individual name is specified, all individuals are reported that are at 
# least a parent once
extract_parents_pedigree <- function(pedigree, ID_names, unknown_parent = 0) {
  pedigree_inds <- pedigree
  if(!missing(ID_names)) {pedigree_inds <- pedigree[pedigree[[1]] %in% ID_names,]}
    parents <- unique(unname(unlist(pedigree_inds[, 3:4])))
  parents <- parents[parents %notin% unknown_parent]
  
  ped_parents <- extract_individuals_pedigree(pedigree, ID_names = parents)
  if(length(parents) != nrow(ped_parents)) {
    warning("Not all parents are listed as individuals in the pedigree object.")
  }
  
  return(ped_parents)
}
# extract_parents_pedigree(pedigree)
# extract_parents_pedigree(pedigree = pedigree, ID_names = c("G12", "G65"))

# extracts all the ancestors of any individual by tracing back all ancestors
# in the pedigree object
# 'unknown_parent' tells how unknown parents are coded (e.g. = or NA)
# if no ID names are specified, the function looks which individuals are ancestors
# (in other words are a parent of at least one individual)
# basically same output as 'extract_parents_pedigree(pedigree)', but slower
extract_ancestors_pedigree <- function(pedigree, ID_names, unknown_parent = 0) {
  if(missing(ID_names)) {
    ID_names <- pedigree[pedigree[[3]] != unknown_parent & 
                          pedigree[[4]] != unknown_parent,1][[1]]
  }
  
  ancestors <- extract_individuals_pedigree(pedigree, ID_names)[, 3:4]
  ancestors <- unique(unname(unlist(ancestors)))
  ancestors <- ancestors[ancestors %notin% unknown_parent]
  
  if(length(ancestors) != 0) {
    for(i in 1:length(ancestors)) {
      ancestors <- c(ancestors, extract_ancestors_pedigree(pedigree, ancestors[i], 
                                                 unknown_parent = unknown_parent))
    }
    return(unique(ancestors))
  }
  
  else {
    return(unique(ancestors))
  }
}
# extract_ancestors_pedigree(pedigree = pedigree, ID_names = c("G13", "G57"))
# extract_ancestors_pedigree(pedigree = pedigree, unknown_parent = NA)

# this function extracts all the offspring of specified parents
# if no parents_ID is supplied, it will return all individuals that are offspring
# of someone (i.e. the pedigree apart from founders)
extract_offspring_pedigree <- function(pedigree, parents_ID){
  if(missing(parents_ID)) {parents_ID <- extract_parents_pedigree(pedigree = pedigree)[[1]]}
  inds_pos <- unique(c(which(pedigree[[3]] %in% parents_ID)
                       , which(pedigree[[4]] %in% parents_ID)))
  return(pedigree[inds_pos,])
}
# extract_offspring_pedigree(pedigree, c("P2", "P21"))
# extract_offspring_pedigree(pedigree)

# this function takes a pedigree as input and returns a data table showing
# the parents of each family
# of 'per_parent' is set to TRUE, a list will be returned that lists in which
# families each parent is involved
extract_family_parents_pedigree <- function(pedigree, per_parent = FALSE){
  
  out <- unique(pedigree[,-1])
  if(per_parent) {
    unq_parents <- unique(unlist(out[,2:3]))
    # this long expression has to be to cut off all the NAs
    unq_parents <- 
      unq_parents[match(pedigree[[1]]
                        , unq_parents)[!is.na(match(pedigree[[1]], unq_parents))]]
    l <- list()
    for(p in unq_parents) {
      l <- c(l,list(c(out[[1]][out[[2]] %in% p], out[[1]][out[[3]] %in% p])))
    }
    names(l) <- unq_parents
    return(l)
  }
  return(out)
}
#extract_family_parents_pedigree(pedigree, per_parent = FALSE)
#extract_family_parents_pedigree(pedigree, per_parent = TRUE)

# this function takes a pedigree and family names and returns a list
# of pedigrees each according to one family + their respective parents
# 
# the family is extrascted based on th family name, not based on teh parents
# this function works very similar to extract_family_marker_table
# but also returns lists if family names are specified
# if it is wished that only one pedobject should be given as output
# you can run rlist::list.rbind(a) on the output 
extract_family_pedigree <- function(pedigree, family_name, negate = FALSE){
  if(missing(family_name)) {
    out <- lapply(unique(pedigree[[2]]), function(i) {
      extract_family_pedigree(pedigree = pedigree, family_name = i)
    })
    
    # this is here to remove one list layer
    for(i in 1:length(out)){out[i] <- out[[i]]}
    
    which_families <- unlist(lapply(out, function(i){nrow(i) > 0}))
    names(out) <- unique(pedigree[[2]])
    
    out <- out[which_families]
  } else {
    if(negate){
      family_name <- unique(unlist(pedigree[!(pedigree[[2]] %in% family_name),2]))
    }
    ls <- list()
    for(i in 1:length(family_name)){
      namse <- c()
      inds_ID <- pedigree[pedigree[[2]] == family_name[i],][[1]]
      
      father <- extract_individuals_pedigree(pedigree, inds_ID[1])[[3]]
      mother <- extract_individuals_pedigree(pedigree, inds_ID[1])[[4]]
      namse <- c(mother, father, inds_ID)
      ls[[i]] <- extract_individuals_pedigree(pedigree = pedigree, ID_names = namse)
    }
    out <- ls
    names(out) <- family_name
  }
  return(out)
}
# a <- extract_family_pedigree(pedigree = ped, family_name = c('Fam01', 'Fam05'))
# extract_family_pedigree(pedigree = ped, family_name = c('Fam01', 'Fam05'), negate = TRUE)
# a <- extract_family_pedigree(pedigree = ped)

# this function samples offspring randomly per family according to the fraction 
# provided
# Not that the fraction is multiplied with the number of offspring per family
# to get the number of individuals to be sampled. If that doesn't result in an 
# integer, the number if rounded to the next higher integer.
#
# The function has no knowledge about what individuals are HD so it actually
# just samples random offspring
extract_random_HD_offspring_per_fam <- function(pedigree, fraction_HD_offspring_per_family){
  ls_fam_ped <- extract_family_pedigree(pedigree = pedigree)
  parent_IDs <- unlist(extract_family_parents_pedigree(pedigree = pedigree, per_parent = FALSE)[,2:3])
  
  HD_inds <- c()
  for(i in 1:length(ls_fam_ped)){
    IDs <- ls_fam_ped[[i]][[1]]
    offspring_ids <- IDs[!(IDs %in% parent_IDs)]
    # ceiling is rounding up to the next highest integer
    # say, there are only 5 offspring in the family and 10% should be HD, 
    # do I pick 0 or 1?
    HD_inds <- c(HD_inds, sample(offspring_ids, replace = FALSE,
                                 size = ceiling(
                                   fraction_HD_offspring_per_family*length(offspring_ids)
                                 ))
    )
  }
  return(HD_inds)
}
#extract_random_HD_offspring_per_fam(pedigree = pedigree, fraction_HD_offspring_per_family = 0.3)


# extracts all the rows in a mapo object that are on the chromosome specified
# with chrom_number
extract_marker_mapo <- function(mapo, marker_IDs, negate = FALSE) {
  if(negate) {mapo <- mapo[mapo[[1]] %notin% marker_IDs,]}
  else {mapo <- mapo[mapo[[1]] %in% marker_IDs,]}
  return(mapo)
}
#extract_marker_mapo(map, "C1M3")
# extract_marker_mapo(map, map[map[[2]] == 1, 1]) # to get all the names of the 
# markers of the first chromosome

# this function extracts the informative markers of the genotype information from
# a marker table
# strictly speaking, this function does not only extract informative markers but 
# also partially informative marker according to the definition of Gonen et al 2018
# so all markers for which the two parents are not identical (it deletes missing markers)
#
# NOTE: you can supply a marker table with more than two individuals. The function 
# will behave the same but it is more difficult to interprete the 'informativeness' then
extract_informative_marker_mapo <- function(mapo, marker_table){
  marker_table <- order_marker_table_by_mapo(mapo = mapo, marker_table = marker_table)
  mapo <- extract_marker_mapo(mapo = mapo, marker_IDs = marker_table[[1]])
  # this filters out all markers that are 9
  miss_m <- apply(marker_table[,-1], MARGIN = 1, function(i) {any(i == 9)})
  marker_table <- marker_table[!miss_m,]
  mapo <- mapo[!miss_m,]
  
  # all missing data is filtered out.
  # this here checks if at a locus are ONLY heterozygous or homozygous individuals
  # if length is 1, it means that all inds are either 0, or 1, or 2
  e <- unlist(lapply(apply(marker_table[,-1], MARGIN = 1, unique), length)) > 1
  
  # uninformative is if both parent are homozygous for the same allele or are both heterozygous
  return(mapo[e,])
}
# infor <- extract_informative_marker_mapo(mapo, marker_table = mt[,c("marker", "P1", "P2")])


# extracts all rows that are mentioned in marker_IDS
# if negate is TRUE, all rows are extracted, that are NOT mentioned
extract_marker_marker_table <- function(marker_table, marker_IDs, negate = FALSE) {
  if(!negate) {return(marker_table[marker_table[[1]] %in% marker_IDs,])}
  else {return(marker_table[marker_table[[1]] %notin% marker_IDs,])}
}
#extract_marker_marker_table(marker_table, "C1M1")
#extract_marker_marker_table(im_fg, marker_IDs = map_LD3[[1]], negate = TRUE)

# extracts the columns of the individuals specified in inds out of a marker_table
extract_individuals_marker_table <- function(marker_table, ID_inds, negate = FALSE) {
  if(!negate) {
    mt <- subset(marker_table, select = colnames(marker_table) %in% ID_inds)
    mt <- data.table(marker = marker_table[[1]], mt)
  } else {
    pos <- which(colnames(marker_table) %notin% ID_inds)
    mt <- marker_table[, ..pos]
  }
  return(mt)
}
# extract_individuals_marker_table(marker_table, pedigree[[1]])


# provided a family name, this function will look up what individuals belong to this
# family and extract  their columns out of a marker table
# if no family name is provided, the function will return a list of marker tables
# each for one family
extract_family_marker_table <- function(pedigree, marker_table, family_name, negate = FALSE) {
  if(missing(family_name)) {
    out <- lapply(unique(pedigree[[2]]), function(i) {
      extract_family_marker_table(pedigree = pedigree, 
                                  marker_table = marker_table, family_name = i) # recursion
    })
    which_families <- unlist(lapply(out, function(i){ncol(i) > 1}))
    names(out) <- unique(pedigree[[2]])
    
    out <- out[which_families]
  } else {
    if(negate){
      family_name <- unique(unlist(pedigree[!(pedigree[[2]] %in% family_name),2]))
      }
    namse <- c()
    for(i in 1:length(family_name)){
      inds_ID <- pedigree[pedigree[[2]] == family_name[i],][[1]]
      
      father <- extract_individuals_pedigree(pedigree, inds_ID[1])[[3]]
      mother <- extract_individuals_pedigree(pedigree, inds_ID[1])[[4]]
      # if names from another family were already extracted, then their names are alrady in namse
      namse <- c(namse, mother, father, inds_ID)
    }
    # unique is needed here becuase it can be that an individual is parent to more than
    # one family
    out <- extract_individuals_marker_table(marker_table = marker_table
                                            , ID_inds = unique(namse))
  }
  return(out)
}
# extract_family_marker_table(pedigree, marker_table, "Fam01")
# extract_family_marker_table(pedigree, marker_table, c("Fam01", "Fam15"))
#o <- extract_family_marker_table(pedigree, marker_table)
#str(o, max.level = 1)


# this function will check which markers are variable in a marker_table
# and return a marker table that only contains the markers
# for which variation was found
# here, 9 (=missing) counts as variation
extract_variable_loci_marker_table <- function(marker_table){
  # this loops over the rows and checks how many states are present
  # there can be at most 3 states: 0,1,3 (actually 4 if you consider 9 as a state)
  # anyway, even if values are missing, so only 0s and 9s occur, the loci
  # will be left in the marker table, BUT if only 0s or only 2s
  # are at the loci, it will be removed
  # PROBLEM: if you only have 1s, it is not fixed but just one
  # state. For that, the line below is used
  n_genos <- apply(marker_table[,-1], MARGIN = 1, unique)
  if('matrix' %in% class(n_genos)){
    n_genos <- rlist::list.flatten(apply(n_genos, MARGIN = 2, list))
  }
  mult_state <- as.logical(unlist(lapply(n_genos, length))-1)
  # checks if at least one heterozygot is present
  het <- unlist(lapply(n_genos, function(i) any(i == 1)))
  # try this yourself to understand this. Basically, T+T = 2 and as.logical(2) = T
  keep_loci <- as.logical(mult_state + het)
  # here, all loci that are fixed are removed
  marker_table <- marker_table[keep_loci,]
  
  return(marker_table)
}
#extract_variable_loci_marker_table(marker_table = marker_table[1:10, 1:10])

# extracts specified markers from a phase object. Obviously requires that the phase 
# object has marker IDs as column names
extract_marker_phase <- function(phase, marker_IDs, negate = FALSE) {
  if(!negate) {cols <- colnames(phase) %in% marker_IDs}
  else {cols <- colnames(phase) %notin% marker_IDs}
  
  return(cbind(phase[,1], phase[,..cols]))
}

# extracts individuals from a phase object
extract_individuals_phase <- function(phase, ID_inds, negate = FALSE) {
  if(!negate) {return(phase[phase[[1]] %in% ID_inds, ])}
  else {return(phase[phase[[1]] %notin% ID_inds, ])}
}


##### make functions

# produces an object in marker_table format in which all the individuals genotyped at 
# high density are still in high density and all others are in low density.
# the decision which marker should be kept as genotyped is made based on the marker
# names in LD_map. All markers not present in that map are set to 9 (can be adjusted
# with set_NA_to) for individuals not mentioned in the vector HD_indiv.
# The function returns the masked table in the same order of markers and individuals
# as they were in the input object
make_masked_marker_table <- function(HD_marker_table, LD_map, HD_indiv, set_NA_to = 9) {
  m <- HD_marker_table
  mt_LD <- m[m[[1]] %in% LD_map[[1]],]
  mt_LD <- subset(mt_LD, select = (colnames(mt_LD) %notin% HD_indiv))
  mt_HD <- subset(m, select = c(TRUE, (colnames(m)[2:ncol(m)] %in% HD_indiv))) 
  # looks funny because I always want to include the first column (marker column)
  masked_mt <- merge(mt_HD, mt_LD, by = "marker", all.x = TRUE, all.y = TRUE, suffixes = "")
  
  if(!is.na(set_NA_to)) {
    masked_mt[is.na(masked_mt)] <- set_NA_to
  }
  
  masked_mt <- order_marker_table_by_mapo(mapo = HD_marker_table[[1]], 
                                          marker_table = masked_mt)
  masked_mt <- order_marker_table_by_pedigree(pedigree = colnames(HD_marker_table)[-1], 
                                              marker_table = masked_mt)
  return(masked_mt)
}

# hd <- extract_parents_pedigree(pedigree)[[1]]
# make_masked_marker_table(marker_table, map_LD3, hd)

# supply a geno object (output of mk_genotype_object()) and a pedigree object
# with genotype ID (col1), pop ID (col2), father ID(col3) and mother ID (col4)
# a ped object in PLINK format is returned
make_ped_object_PLINK <- function(geno, pedigree_object) {
  ped <- pedigree_object
  ped[is.na(ped)] <- 0
  
  ped_file <- rbind(ped[,c(2,1,3,4)])
  ped_file[[5]] <- rep(0, nrow(ped_file))
  ped_file[[6]] <- rep(0, nrow(ped_file))
  colnames(ped_file)[c(5,6)] <- c("SEX", "PHENOTYPE")
  
  lookup <- c("1" = "C",
              "0" = "T",
              "9" = "0")
  
  namse <- colnames(geno)
  dimens <- dim(geno)
  geno <- apply(geno, 2, as.character)
  
  geno <- unname(lookup[geno])
  dim(geno) <- dimens
  colnames(geno) <- namse
  
  return(as.data.table(cbind(ped_file, geno)))
}

# this function takes a vector of marker info from an individual and the fraction
# of markers that should be incorrect
# it then randomly samples which marker should be incorrect.
# Then, it inserts wrong information atthese positions and return back the vector
# Say if the truth is 2, then a 0 or 1 is inserted
# This function has no problem if there are 9s in the vector. These may be used
# to insert incorrect information too
# this function also works, if the fraction is 0
#
# note that the fraction if incorrect markers is not an average but the precise 
# number of incorrect markers per genotype
make_ind_genotype_error <- function(ind, fraction_incorrect){
  ind <- unname(unlist(ind))
  
  if(fraction_incorrect > 0){
    index_marker <- sample(seq(1:length(ind)), size = length(ind)*fraction_incorrect)
    markers_ind_for_error <- ind[index_marker]
    
    markers_ind_error <- c()
    for(i in 1:length(markers_ind_for_error)){
      markers_ind_error[i] <- sample(c(0,1,2)[c(0,1,2) != markers_ind_for_error[i]],
                                     size = 1, 
                                     replace = TRUE, prob = c(0.5,0.5,0.5)[c(0,1,2) != markers_ind_for_error[i]])
    }
    ind[index_marker] <- markers_ind_error
  }
  
  return(ind)
}
#make_ind_genotype_error(ind = gt, fraction_incorrect = 0.2)

# this is a wrapper function of make_ind_genotype_error that is looping over
# marker table
# note that the below function produces errors independent of the state of the 
# truth. I.e. a '0' is as likely as a '1' if the truth is a '2'. In reality,
# it would probably be more likely to get a '1' if the truth is a '2'
make_genotype_marker_errors_marker_table <- function(marker_table, fraction_incorrect){
  marker_table <- as.data.frame(marker_table)
  for(i in 2:ncol(marker_table)){
    marker_table[,i] <- make_ind_genotype_error(ind = marker_table[,i],
                                                fraction_incorrect = fraction_incorrect)
  }
  return(as.data.table(marker_table))
}
# mt <- make_genotype_errors_marker_table(marker_table = mt, fraction_incorrect = 0.0)

# this function induces missing values per genotype
# it induces missing values per genotype (individual)
# the fraction of missing markers is exact and not an average
# However, if the fraction should be 0.35 and there are only 10 markers,
# 3 markers will actually be set to missing instead of 4. So always
# rounded down to an integer
make_ind_genotype_missing <- function(ind, fraction_missing, set_missing_to = 9){
  ind <- unname(unlist(ind))
  
  # giving the positions at which NA should be inserted
  if(fraction_missing > 0){
    index_marker <- sample(seq(1:length(ind)), size = length(ind)*fraction_missing)
    ind[index_marker] <- set_missing_to
  }
  
  return(ind)
}

# this function is a wrapper function of make_ind_genotype_missing that is looping
# over a marker table
make_genotype_marker_missing_marker_table <- function(marker_table, fraction_missing){
  marker_table <- as.data.frame(marker_table)
  for(i in 2:ncol(marker_table)){
    marker_table[,i] <- make_ind_genotype_missing(ind = marker_table[,i],
                                                  fraction_missing = fraction_missing)
  }
  return(as.data.table(marker_table))
}

#mt <- make_genotype_marker_missing_marker_table(marker_table = mt, fraction_missing = 0.3)


################ write functions
# supply a map in PLINK format
# you can convert mapo to PLINK format with:
# pm <- convert_mapo_to_PLINK_map(mapo = map)
write_PLINK_map <- function(PLINK_map_object, file_prefix) {
  fwrite(PLINK_map_object, quote = FALSE, sep = " ", 
         file = paste(file_prefix, ".map", sep = ""), 
         na = NA, col.names = FALSE, row.names = FALSE)
}
#write_PLINK_map(PLINK_map_object = pm, file_prefix = "test_map")

# takes a ped object (as produced by make_ped_object_PLINK and writes that object
# to a file with the ending '.ped'
write_PLINK_ped <- function(ped_obj, file_prefix) {
  fwrite(ped_obj, file = paste(file_prefix, ".ped", sep = ""), 
         col.names = FALSE, sep = " ")
}


######################## 
# zip functions
# fread and the read functions created here are able to directly read in gz zipped files

# this function unzipps a gz zipped file and does NOT remove the original file
gunzip_gzfile <- function(filename, rm = FALSE) {
  new_filename <- str_remove(filename, pattern = ".gz$")
  if(file.exists(new_filename)) {
    print("Overwrote existing file.")
    file.remove(new_filename)
  }
  return(gunzip(filename, remove=rm))
}
# gunzip_gzfile("out.gt.vcf.gz")


# this function gz zips a file and does not remove the original file
gz_zip_gzfile <- function(filename, rm = FALSE) {
  if(file.exists(paste0(filename, ".gz"))) {
    print("Overwrote existing file.")
    file.remove(paste0(filename, ".gz"))
  }
  return(R.utils::compressFile(filename, ext = "gz", remove = rm, FUN = gzfile))
}
# gz_zip_gzfile("out.gt.vcf")

############################################################
# order functions

# this is the first order function for BEAGLE 
# it was developed to deal with marker tables that are present 
# in a list like they would be after imputing 
# family wise
# 
# this function takes the output form familywise imputation with BEAGLE
# Some individuals (e.g. parents) will be imputed several times according to the number
# of families they are involved in
# this function will return a list of marker table and every marker table contains 
# all predictions for one individual
order_ind_imps_lsmt <- function(list_marker_tables) {
  id_nam <- lapply(list_marker_tables, names)
  id_nam <- unique(unlist(id_nam))
  id_nam <- id_nam[!(id_nam %in% 'marker')]
  
  ls_ind_mt <- replicate(length(id_nam), list())
  names(ls_ind_mt) <- id_nam
  
  for(i in 1:length(list_marker_tables)) {
    col_names <- colnames(list_marker_tables[[i]])
    for(j in 2:ncol(list_marker_tables[[i]])) {
      if(length(ls_ind_mt[[match(col_names[j], names(ls_ind_mt))]]) == 0) {
        nc <- c(1,j)
        ls_ind_mt[[match(col_names[j], names(ls_ind_mt))]] <- list_marker_tables[[i]][,..nc]
      } else {
        ls_ind_mt[[match(col_names[j], names(ls_ind_mt))]] <- cbind(ls_ind_mt[[match(col_names[j], names(ls_ind_mt))]], 
                                                                    list_marker_tables[[i]][,..j])
      }
    }
  }
  return(ls_ind_mt)
}
# order_ind_imps_lsmt(IN_mt)

################################################################################
# extract functions

# this function is a bit different then other extract functions
# when you do familywise imputation, e.g. with 'impute_familywise_BEAGLE()'
# you get a list of phase objects as output
# some individuals are a member in more than one family and depending on the data
# quality and how good the imputation is, may be predicted differently depending on the family
# you can use this function then to check which positions have always been imputed
# the same for every family. This function will return a marker table that only shows
# markers at a genotype that have a unique prediction (so were imputed as '0' in
# every family for example). Others are missing. See a workflow example at the
# end of this function definition
extract_unique_imp_lsmt <- function(list_marker_tables){
  ls_ind_mt <- order_ind_imps_lsmt(list_marker_tables = list_marker_tables)
  
  # only extract those marker tables that only have two columns (marker name and 
  # prediction of that individual)
  # they could also be left in but then the for loop later takes way longer
  one_pred <- ls_ind_mt[(lapply(ls_ind_mt, ncol) == 2)]
  w <- rlist::list.cbind(one_pred)
  nc <- -seq(1,ncol(w),2)
  w <- w[,..nc]
  colnames(w) <- names(one_pred)
  
  # now only extract the ones with at least two predictions
  # two predictions are made if that individual is part 
  # of at least two different families
  more_one_pred <- ls_ind_mt[(lapply(ls_ind_mt, ncol) != 2)]
  
  out <- matrix(data = NA, nrow = nrow(IN_mt[[1]]), ncol = length(more_one_pred))
  
  # here, I go in marker table in the list one by one. Within every marker table
  # I check for every marker if all imputations predicted the same
  # if so, I accept this predictin and use it for the out marker table
  # if not, I set it to missing 
  for(i in 1:length(more_one_pred)) {
    marktab <- more_one_pred[[i]]
    unique_prediction <- apply(marktab[,-1], MARGIN = 1, FUN = unique)
    uniq_pred_marker <- which(unlist(lapply(unique_prediction, length)) == 1)
    out[uniq_pred_marker,i] <- unlist(marktab[uniq_pred_marker,2])
  }
  out <- as.data.table(out)
  colnames(out) <- names(more_one_pred)
  out <- cbind(out, w)
  out[is.na(out)] <- 9
  # this is done to ad the marker names to the out marker table
  out <- data.table(list_marker_tables[[1]][,1], out)
  
  return(out)
}
# workflow example
#
# IN_phase <- impute_familywise_BEAGLE(arguments = arguments, marker_table = mt, pedigree = pedigree
#                               , memory_limit_GB = 2, mapo = mapo, verbose = TRUE)
# IN_mt <- lapply(IN_phase, convert_phase_to_marker_table)
#
# unq_mt <- extract_unique_imp_lsmt(IN_mt)
# unq_mt <- order_marker_table_by_mapo(marker_table = unq_mt, mapo = colnames(mt)[-1])
# unq_mt will be in the same order as mt (the marker table used as input for imputation)
# for every individual at every marker, it will only conatin information if all
# predictions (for every family) predicted the same. If not, it will be missing (9)


################################################################################
# write functions

# provide a map in mapo format. The genetic distance is not important. Numbers from 
# 1 to n are given for bp positions for n markers in increments of 1
# provide a marker_table. A pedigree can be provided but is not required
# make sure that the map and the marker_table object contain the same markers
# e.g. like this:
# marker_table <- extract_marker_marker_table(marker_table = marker_table, marker_IDs = map[[1]])
# map <- extract_marker_mapo(mapo = map, marker_IDs = marker_table[[1]])
# marker_table <- order_marker_table_by_pedigree(pedigree, marker_table)
# if ID_unknown_chrom is NA, it is looked up in the map where markers are on 
# chromosome NA. These are set to '.' in the vcf. You can also specify 99,
# or c(14,99)
# BEAGLE apparently does not accept '.' for
write_marker_table_to_vcf <- function(mapo, marker_table, pedigree = NULL, 
                                      outfile_name = "outfile.vcf", verbose = TRUE,
                                      ID_unknown_chrom = NA, gzipped = TRUE) {
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start), "\n")}
  
  map <- mapo
  
  if(any(marker_table[[1]] %notin% map[[1]]) & verbose) {
    stop("The marker_table object contains markers that are not present in mapo.\n")
  }
  
  if(any(map[[1]] %notin% marker_table[[1]])) {
    warning("The mapo object contains markers not present in the marker_table.\n")
  }
  
  map <- extract_marker_mapo(map, marker_IDs = marker_table[[1]])
  map <- order_mapo_by_position(map)
  marker_table <- order_marker_table_by_mapo(mapo = map, marker_table = marker_table)
  
  dt <- as.data.table(matrix(NA, ncol = 8, nrow = nrow(map)))
  colnames(dt) <- c("#CHROM","POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  dt$POS <- unlist(lapply(table(map[[2]]), function(i) {seq(1, i, 1)}))
  
  #POS <- unlist(unname(lapply(table(map[[2]]), function(i) {seq(1, i, 1)})))
  
  map[map[[2]] %in% ID_unknown_chrom,2] <- "."
  dt$`#CHROM` <- map[[2]]
  dt$ID <- map[[1]]
  
  dt$REF <- rep("T", nrow(marker_table))
  dt$ALT <- rep("C", nrow(marker_table))
  dt$QUAL <- rep(".", nrow(marker_table))
  dt$FILTER <- rep(".", nrow(marker_table))
  dt$FORMAT <- rep("GT", nrow(marker_table))
  
  #colnames(dt)[-c(1:8)] <- colnames(marker_table)[-1]
  
  mt <- marker_table
  mt[mt == 9] <- NA
  AF <- rowMeans(mt[,-1],na.rm = TRUE)/2
  NSM <- unname(apply(mt, MARGIN = 1, function(i) {sum(is.na(i[-1]))}))
  
  mt[mt == 0] <- "0/0"
  mt[mt == 1] <- "0/1"
  mt[mt == 2] <- "1/1"
  
  NSP <- apply(mt[,-1], MARGIN = 1, function(i) {sum(str_detect(i, pattern = "\\|"), na.rm = TRUE)})
  mt[is.na(mt)] <- "./."
  
  dt$INFO <- paste("AF=", AF, ";NSP=", NSP, ";NSM=", NSM, sep = "")
  
  vcf_table <- cbind(dt, mt[,-1])
  
  if(verbose) {cat("Finished creating the vcf table. Procede writing vcf file.\n")}
  ########## Now the part where writing the vcf begins starts
  
  fwrite(list("##fileformat=VCFv4.2"), file = outfile_name)
  if(verbose){cat("Started writing the vcf file.\n")}
  
  fwrite(list(paste0("##fileDate=", str_replace(str_remove_all(start, pattern = "[:punct:]*"),
                                                pattern = " ", replacement = "_"))),
         append = TRUE, file = outfile_name)
  
  metainfo_info <- c(
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Reference Allele Frequency">',
    '##INFO=<ID=NSP,Number=1,Type=Integer,Description="Number of Samples Phased">',
    '##INFO=<ID=NSM,Number=1,Type=Integer,Description="Number of Samples Masked">'
  )
  lapply(metainfo_info, function(i) {fwrite(list(i), append = TRUE, file = outfile_name,
                                            quote = FALSE)})
  
  fwrite(list('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'),
         append = TRUE, file = outfile_name, quote = FALSE)
  
  if(!is.null(pedigree)){
    ped <- apply(pedigree, MARGIN = 1, function(row) {
      paste0("##PEDIGREE=<Descendant=", row[1], ",Parent_1=", row[3], ",Parent_2=", row[4],
             ",Family=", row[2], ">")
    })
    fwrite(list(ped), append = TRUE, file = outfile_name, quote = FALSE)
  }
  
  if(verbose) {cat("Finished writing meta-information. Procede writing vcf table.\n")}
  fwrite(vcf_table, append = TRUE, file = outfile_name, quote = FALSE, sep = "\t",
         col.names = TRUE)
  if(gzipped) {
    gz_zip_gzfile(outfile_name, rm = TRUE)
    outfile_name <- paste0(outfile_name, ".gz")
  }
  
  if(verbose){
    cat('The process needed ',format(Sys.time() - start),
        'and the file is',as.numeric(file.size(outfile_name))/1000000,'Mb',
        'large.\n', 'The file is called:', outfile_name,".\n")
  }
}

# write_marker_table_to_vcf(mapo = map, marker_table)


#### check if it is really necessary to use the map so often

# writes a phase object to vcf file in the working directory
# does the same as the write_marker_table_to_vcf but for phased data so the data
# in the vcf will also be phased
write_phase_to_vcf <- function(mapo, phase, pedigree = NULL, 
                               outfile_name = "outfile.vcf", gzipped = TRUE, verbose = TRUE) {
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start),'\n')}
  
  map <- mapo
  
  if(any(colnames(phase)[-1] %notin% map[[1]]) & verbose) {
    stop("The phase object contains more markers than mapo.\n")
  }
  
  if(any(map[[1]] %notin% colnames(phase)[-1])) {
    warning("The mapo object contains markers not present in phase.\n")
  }
  
  map <- extract_marker_mapo(map, marker_IDs = colnames(phase)[-1])
  map <- order_mapo_by_position(map)
  phase <- order_phase_by_mapo(mapo = map, phase = phase)
  
  dt <- as.data.table(matrix(NA, ncol = 8, nrow = nrow(map)))
  colnames(dt) <- c("#CHROM","POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  dt$POS <- unlist(lapply(table(map[[2]]), function(i) {seq(1, i, 1)}))
  
  #POS <- unlist(unname(lapply(table(map[[2]]), function(i) {seq(1, i, 1)})))
  
  map[map[[2]] == 99, 2] <- "."
  dt$`#CHROM` <- map[[2]]
  dt$ID <- map[[1]]
  
  dt$REF <- rep("T", ncol(phase)-1)
  dt$ALT <- rep("C", ncol(phase)-1)
  dt$QUAL <- rep(".", ncol(phase)-1)
  dt$FILTER <- rep(".", ncol(phase)-1)
  dt$FORMAT <- rep("GT", ncol(phase)-1)
  
  inds <- seq(1, nrow(phase), 2)
  gt_names <- phase[[1]][inds]
  phase <- phase[,-1]
  
  p <- phase[seq(1, nrow(phase),2)] + phase[seq(2, nrow(phase),2)] # this is needed to 
  # determine NSM. 9 + 9 is 18 so all that are 18 are masked
  
  phase[phase == 9] <- NA
  AF <- (colSums(phase,na.rm = TRUE)/(nrow(phase)/2))/2
  NSM <- apply(p, MARGIN = 2, function(i) {sum(i == 18)})
  rm(p)
  
  phase[is.na(phase)] <- "."
  
  NSP <- rep(nrow(phase)/2, nrow(dt))
  
  dt$INFO <- paste("AF=", AF, ";NSP=", NSP, ";NSM=", NSM, sep = "")
  
  ########## Now the part where writing the metainformation begins starts
  fwrite(list("##fileformat=VCFv4.2"), file = outfile_name)
  if(verbose){cat("Started writing the vcf file.\n")}
  
  fwrite(list(paste0("##fileDate=", stringr::str_replace(
    stringr::str_remove_all(start, pattern = "[:punct:]*"),
    pattern = " ", replacement = "_"))),
    append = TRUE, file = outfile_name)
  
  metainfo_info <- c(
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Reference Allele Frequency">',
    '##INFO=<ID=NSP,Number=1,Type=Integer,Description="Number of Samples Phased">',
    '##INFO=<ID=NSM,Number=1,Type=Integer,Description="Number of Samples Masked">'
  )
  lapply(metainfo_info, function(i) {
    fwrite(list(i), append = TRUE, file = outfile_name, quote = FALSE)
  })
  
  fwrite(list('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'),
         append = TRUE, file = outfile_name, quote = FALSE)
  
  if(!is.null(pedigree)){
    ped <- apply(pedigree, MARGIN = 1, function(row) {
      paste0("##PEDIGREE=<Descendant=", row[1], ",Parent_1=", row[3], ",Parent_2=", row[4],
             ",Family=", row[2], ">")
    })
    fwrite(list(ped), append = TRUE, file = outfile_name, quote = FALSE)
  }
  
  if(verbose) {cat("Finished writing metainformation. Procede writing the vcf table.\n")}
  
  fwrite(as.list(c(colnames(dt), gt_names)), quote = FALSE, file = outfile_name, 
         append = TRUE, sep = "\t")
  rm(gt_names)
  
  maternal <- seq(1, nrow(phase),2)
  paternal <- seq(2, nrow(phase),2)
  # now I am putting the phase object into a shape needed for the vcf table
  p_new <- matrix(NA, ncol = ncol(phase), nrow = (nrow(phase)/2))
  for(i in 1:nrow(p_new)){
    p_new[i,] <- c(paste0(phase[maternal[i],],'|', 
                          phase[paternal[i],]))
  }
  rm(phase)
  p_new <- t(p_new)
  vcf_table <- cbind(dt, p_new)
  rm(dt)
  rm(p_new)
  
  fwrite(vcf_table, append = TRUE, file = outfile_name, quote = FALSE, sep = "\t",
         col.names = FALSE)
  if(verbose) {cat('Finished writing the vcf table.\n')}
  
  if(gzipped) {
    gz_zip_gzfile(outfile_name, rm = TRUE)
    outfile_name <- paste0(outfile_name, ".gz")
  }
  
  if(verbose){
    cat('The process needed ',format(Sys.time() - start),
        'and the file is',as.numeric(file.size(outfile_name))/1000000,'Mb',
        'large.\n')
  }
}

# test_map <- extract_marker_mapo(map, marker_IDs = colnames(phase))
# write_phase_to_vcf(mapo = test_map, phase = phase)


# a function that takes a marker_table and a phase object and produces a vcf
# object with both combined

# marker table can contain individuals that are also present in phase
# this function should look which are in phase and take those from phase

# could implement that it is possible to give a marker_table and a phase object
# with different markers and still produces a vcf file and orders markers based 
# on supplied mapo object

# you should supply a marker_table object that contains at least one individual 
# that is not present in the phase object
# the final vcf file will have only those markers that are present in the phase object
write_phase_and_mt_to_vcf <- function(mapo, phase, marker_table, pedigree = NULL, 
                                      outfile_name = "outfile.vcf", verbose = TRUE,
                                      gzipped = TRUE) {
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start))}
  
  map <- mapo
  
  if(any(colnames(phase)[-1] %notin% map[[1]]) & verbose) {
    stop("The phase object contains more markers than mapo.\n")
  }
  
  if(any(map[[1]] %notin% colnames(phase)[-1])) {
    warning("The mapo object contains markers not present in phase.\n")
  }
  
  map <- extract_marker_mapo(map, marker_IDs = colnames(phase)[-1])
  map <- order_mapo_by_position(map)
  phase <- order_phase_by_mapo(mapo = map, phase = phase)
  marker_table <- extract_marker_marker_table(marker_table, marker_IDs = colnames(phase)[-1]) # so only
  # markers that are present in phase will be taken from the marker_table
  
  marker_table <- order_marker_table_by_mapo(marker_table = marker_table, mapo = map)
  
  if(any(colnames(phase)[-1] != marker_table[[1]]) & verbose) {
    stop("The phase and marker_table object do not have the same marker order.\n")
  } # should be false # this checks is the marker 
  # name sin phase are identical to the marker names in marker_table
  
  # trimming the marker_table
  gt_from_mt <- !(colnames(marker_table) %in% phase[[1]])
  mt <- subset(marker_table, select = gt_from_mt)
  mt <- mt[,-1]
  
  dt <- as.data.table(matrix(NA, ncol = 8, nrow = nrow(map)))
  colnames(dt) <- c("#CHROM","POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  dt$POS <- unlist(lapply(table(map[[2]]), function(i) {seq(1, i, 1)}))
  
  #POS <- unlist(unname(lapply(table(map[[2]]), function(i) {seq(1, i, 1)})))
  
  map[map[[2]] == 99, 2] <- "."
  dt$`#CHROM` <- map[[2]]
  dt$ID <- map[[1]]
  
  dt$REF <- rep("T", nrow(map))
  dt$ALT <- rep("C", nrow(map))
  dt$QUAL <- rep(".", nrow(map))
  dt$FILTER <- rep(".", nrow(map))
  dt$FORMAT <- rep("GT", nrow(map))
  
  inds <- seq(1, nrow(phase), 2)
  gt_names <- c(phase[[1]][inds], colnames(mt))
  phase <- phase[,-1]
  
  p <- phase[seq(1, nrow(phase),2)] + phase[seq(2, nrow(phase),2)] # this is needed to 
  # determine NSM. 9 + 9 is 18 so all that are 18 are masked
  
  NSM_ph <- apply(p, MARGIN = 2, function(i) {sum(i == 18)})
  NSM_mt <- apply(mt, MARGIN = 1, function(i) {sum(i == 9)})
  if(length(NSM_ph) == 0) {NSM_ph <- rep(0, nrow(map))} 
  if(length(NSM_mt) == 0) {NSM_mt <- rep(0, nrow(map))}
  # if nothing is in phase or marker_table, it will be integer
  # probably worrying too much
  
  NSM <- NSM_ph + NSM_mt
  rm(p)
  
  phase[phase == 9] <- NA
  mt[mt == 9] <- NA
  
  tot <- ncol(mt) + (nrow(phase)/2)
  AF_ph <- (colSums(phase,na.rm = TRUE)/(nrow(phase)/2))/2
  AF_mt <- (rowSums(mt, na.rm = TRUE)/ncol(mt))/2
  if(length(AF_ph) == 0) {AF_ph <- rep(0, nrow(map))} 
  if(length(AF_mt) == 0) {AF_mt <- rep(0, nrow(map))}
  
  AF <- AF_ph * ((nrow(phase)/2)/tot) + AF_mt * (ncol(mt)/tot)
  
  phase[is.na(phase)] <- "."
  NSP <- rep(nrow(phase)/2, nrow(dt))
  
  dt$INFO <- paste("AF=", AF, ";NSP=", NSP, ";NSM=", NSM, sep = "")
  
  ########## Now the part where writing the metainformation begins starts
  fwrite(list("##fileformat=VCFv4.2"), file = outfile_name)
  if(verbose){cat("Started writing the vcf file.\n")}
  
  fwrite(list(paste0("##fileDate=", str_replace(str_remove_all(start, pattern = "[:punct:]*"),
                                                pattern = " ", replacement = "_"))),
         append = TRUE, file = outfile_name)
  
  metainfo_info <- c(
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Reference Allele Frequency">',
    '##INFO=<ID=NSP,Number=1,Type=Integer,Description="Number of Samples Phased">',
    '##INFO=<ID=NSM,Number=1,Type=Integer,Description="Number of Samples Masked">'
  )
  lapply(metainfo_info, function(i) {
    fwrite(list(i), append = TRUE, file = outfile_name, quote = FALSE)
  })
  
  fwrite(list('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'),
         append = TRUE, file = outfile_name, quote = FALSE)
  
  if(!is.null(pedigree)){
    ped <- apply(pedigree, MARGIN = 1, function(row) {
      paste0("##PEDIGREE=<Descendant=", row[1], ",Parent_1=", row[3], ",Parent_2=", row[4],
             ",Family=", row[2], ">")
    })
    fwrite(list(ped), append = TRUE, file = outfile_name, quote = FALSE)
  }
  
  if(verbose) {cat("Finished writing metainformation. Procede writing the vcf table.\n")}
  
  fwrite(as.list(c(colnames(dt), gt_names)), quote = FALSE, file = outfile_name, 
         append = TRUE, sep = "\t")
  rm(gt_names)
  
  mt[mt == 0] <- "0/0"
  mt[mt == 1] <- "0/1"
  mt[mt == 2] <- "1/1"
  mt[mt == NA] <- "./."
  
  maternal <- seq(1, nrow(phase),2)
  paternal <- seq(2, nrow(phase),2)
  
  # this is were the magic happens. Since making a big table in R first and then 
  # writing it to the outfile was not possible on my pc (too much RAM needed),
  # this is the solution: this loop loops over the markers in the phase object,
  # produces the line of the vcf file for that marker and then writes it. Then this
  # is repeated for the next marker. So actually it is like writing line by line to 
  # the vcf file instead of one large table at once
  for(i in 1:ncol(phase)){
    marker <- (unname(c(dt[i,], 
                        paste0(phase[[i]][maternal], "|", phase[[i]][paternal]),
                        mt[i,])))
    
    fwrite(marker, append = TRUE, file = outfile_name, quote = FALSE, sep = "\t")
    
    if(verbose) {cat("Finished", (i/ncol(phase))*100, "%.\n")}
  }
  
  rm(phase)
  rm(dt)
  
  if(verbose) {cat("Finished writing vcf table.\n")}
  if(gzipped) {gz_zip_gzfile(outfile_name, rm = TRUE)}
  
  if(verbose){
    cat('The process needed ',format(Sys.time() - start),
        'and the file is',as.numeric(file.size(outfile_name))/1000000,'Mb',
        'large.\n')
  }
}
# assuming maker_table and phase have the same individuals
# phase <- phase[-c(1:8),] # to delete some individuals in phase
# write_phase_and_mt_to_vcf(mapo = map, phase = phase, marker_table = marker_table)

################################################################################
# read functions

# or make_arguments_createlib_API2 and make_arguments_impute_API2
# works for both
read_arguments <- function(file_name) {
  IN <- readLines(file_name)[-1]
  arg <- paste(IN, collapse = " ")
  return(arg)
}
#read_arguments("arguments_API2.txt")
#read_arguments(file_name = "arguments_BEAGLE.txt")

# reads vcf to marker_table
# also able to read .vcf.gz (zip compressed files) directly to marker_table
read_vcf_to_marker_table <- function(file, verbose = TRUE, set_missing_gt_to = 9) {
  # this function is mainly copied from Johannes class (class 8) and modified in the end
  
  # record start time
  start <- Sys.time()
  
  # read in file
  IN <- data.table::fread(file,
                          skip = '#CHROM',
                          sep = '\t', sep2=NULL,
                          na.strings = '.')
  
  nheader_cols <- sum(!str_detect(IN[1,], "[/\\|]"))
  
  # generate initial structure with map
  VCF <- list(MAP = IN[,1:nheader_cols])
  
  # define function to structure your info field
  structure_info <- function(x){
    x <- stringr::str_split(x,';') # split info filed by elements
    
    # define a function which seperates name and content and produces named list
    name_list <- function(x){
      if(is.na(x[[1]])) {
        out <- '.'
      } else {
        out <- as.list(as.numeric(str_split(x,'=',simplify = TRUE)[,2])) # extract content of info entry
        names(out) <- str_split(x,'=',simplify = TRUE)[,1] # extract name of info entry
        names(out) <- str_extract(x,'^[^=]+') # alternative using str_extract
      }
      return(out)
    }
    
    # apply function to list elements - output is a list of length(x)
    lapply(x,name_list)
  }
  
  VCF$MAP[,INFO := structure_info(INFO)]
  
  # remove unnecessary columns from IN and transform to matrix
  IN <- as.matrix(IN[,-1:-nheader_cols])
  
  # extract samples
  VCF$SAMPLES <- colnames(IN)
  
  # extract genotypes
  VCF$GT <- str_sub(IN,1,3)
  lookup <- c('0/0' = 0,
              '0|0' = 0,
              '0/1' = 1,
              '1/0' = 1,
              '0|1' = 1,
              '1|0' = 1,
              '1/1' = 2,
              '1|1' = 2,
              './.' = set_missing_gt_to,
              '.|.' = set_missing_gt_to,
              '.|0' = set_missing_gt_to,
              './0' = set_missing_gt_to,
              '.|1' = set_missing_gt_to,
              './1' = set_missing_gt_to,
              '0|.' = set_missing_gt_to,
              '0/.' = set_missing_gt_to,
              '1|.' = set_missing_gt_to,
              '1/.' = set_missing_gt_to)
  VCF$GT <- unname(lookup[VCF$GT])
  dim(VCF$GT) <- dim(IN)
  
  ###### The following is for if you have more information per individual. Maybe 
  # build a function that checks what formats are present based on the FORMAT column
  # and then executes the respective commands
  
  # extract AD of alternative allele
  # VCF$ALT_AD <- str_split(IN,':',simplify = TRUE)[,2]
  # VCF$ALT_AD <- str_split(VCF$ALT_AD,',',simplify = TRUE)[,2]
  # dim(VCF$ALT_AD) <- dim(IN)
  # storage.mode(VCF$ALT_AD) <- 'integer'
  
  # extract depth
  # VCF$DP <- str_split(IN,':',simplify = TRUE)[,3]
  # dim(VCF$DP) <- dim(IN)
  # storage.mode(VCF$DP) <- 'integer' # can produce NA when '.'
  
  # assign class VCF to final object
  class(VCF) <- 'VCF'
  
  marker_table <- VCF$GT
  
  marker_table <- apply(marker_table, MARGIN = 2, as.integer)
  marker_table <- data.table(VCF$MAP$ID, marker_table)
  colnames(marker_table) <- c("marker", VCF$SAMPLES)
  
  # add runtime message
  if(verbose){
    cat('The process needed ',format(Sys.time() - start),
        'and the object is',format(object.size(marker_table),'Mb'),
        'large.\n')
  }
  
  # return final object
  return(marker_table)
}
# mt <- read_vcf_to_marker_table("outfile.vcf)


# reads in a vcf file and creates a phase object based on it
# the vcf file must be in the working directory. 
# not sure if function still works (on personal laptop) if "set_missing_gt_to"
# is set to any other than an integer
# set_missing_gt_to a parameter with which it can be told to what missing genotypes 
# should be set to. Something is missing, if e.g. "0/1" (unphased) or for first haplotype
# if (.|1) (phased but first marker info missing). In the latter case, the second
# haplotype would still be 1
read_vcf_to_phase <- function(file, verbose = TRUE, set_missing_gt_to = 9) {
  # this function is mainly copied from Johannes class (class 8) and modified in the end
  
  # record start time
  start <- Sys.time()
  
  # read in file
  IN <- fread(file,
              skip = '#CHROM',
              sep = '\t', sep2=NULL,
              na.strings = '.')
  
  nheader_cols <- sum(!str_detect(IN[1,], "[/\\|]"))
  
  # generate initial structure with map
  VCF <- list(MAP = IN[,1:nheader_cols])
  
  # define function to structure your info field
  structure_info <- function(x){
    x <- str_split(x,';') # split info filed by elements
    
    # define a function which seperates name and content and produces named list
    name_list <- function(x){
      if(is.na(x[[1]])) {
        out <- '.'
      } else {
        out <- as.list(as.numeric(str_split(x,'=',simplify = TRUE)[,2])) # extract content of info entry
        names(out) <- str_split(x,'=',simplify = TRUE)[,1] # extract name of info entry
        names(out) <- str_extract(x,'^[^=]+') # alternative using str_extract
      }
      return(out)
    }
    
    # apply function to list elements - output is a list of length(x)
    lapply(x,name_list)
  }
  
  VCF$MAP[,INFO := structure_info(INFO)]
  
  nheader_cols <- sum(!str_detect(IN[1,], "[/\\|]"))
  
  # remove unnecessary columns from IN and transform to matrix
  IN <- as.matrix(IN[,-1:-nheader_cols])
  
  # extract samples
  VCF$SAMPLES <- colnames(IN)
  
  # extract genotypes
  VCF$GT <- str_sub(IN,1,3)
  
  # this lookup table is for the maternal haplotype
  lookup_mat <- c('0/0' = 0,
                  '0|0' = 0,
                  '0/1' = set_missing_gt_to,
                  '0|1' = 0,
                  '1/0' = set_missing_gt_to,
                  '1|0' = 1,
                  '1/1' = 1,
                  '1|1' = 1,
                  './.' = set_missing_gt_to,
                  '.|.' = set_missing_gt_to, 
                  '0/.' = set_missing_gt_to,
                  '0|.' = 0,
                  '1/.' = set_missing_gt_to,
                  '1|.' = 1,
                  './0' = set_missing_gt_to,
                  '.|0' = set_missing_gt_to,
                  './1' = set_missing_gt_to,
                  '.|1' = set_missing_gt_to)
  
  # this lookup table is for the paternal haplotype
  lookup_pat <- c('0/0' = 0,
                  '0|0' = 0,
                  '0/1' = set_missing_gt_to,
                  '0|1' = 1,
                  '1/0' = set_missing_gt_to,
                  '1|0' = 0,
                  '1/1' = 1,
                  '1|1' = 1,
                  './.' = set_missing_gt_to,
                  '.|.' = set_missing_gt_to, 
                  '0/.' = set_missing_gt_to,
                  '0|.' = set_missing_gt_to,
                  '1/.' = set_missing_gt_to,
                  '1|.' = set_missing_gt_to,
                  './0' = set_missing_gt_to,
                  '.|0' = 0,
                  './1' = set_missing_gt_to,
                  '.|1' = 1)
  
  phase_mat <- unname(lookup_mat[VCF$GT])
  phase_mat <- matrix(data = as.integer(phase_mat), nrow = ncol(IN), byrow = TRUE)
  
  phase_pat <- unname(lookup_pat[VCF$GT])
  phase_pat <- matrix(data = as.integer(phase_pat), nrow = ncol(IN), byrow = TRUE)
  rm(IN)
  
  phase <- matrix(data = as.integer(NA), nrow = nrow(phase_mat)*2, ncol = ncol(phase_mat))
  
  maternal <- seq(1, nrow(phase), 2)
  paternal <- seq(2, nrow(phase), 2)
  
  phase[maternal,] <- phase_mat
  rm(phase_mat)
  phase[paternal,] <- phase_pat
  rm(phase_pat)
  
  phase <- as.data.table(phase)
  
  ###### The following is for if you have more information per individual. Maybe 
  # build a function that checks what formats are present based on the FORMAT column
  # and then executes the respective commands
  
  # extract AD of alternative allele
  # VCF$ALT_AD <- str_split(IN,':',simplify = TRUE)[,2]
  # VCF$ALT_AD <- str_split(VCF$ALT_AD,',',simplify = TRUE)[,2]
  # dim(VCF$ALT_AD) <- dim(IN)
  # storage.mode(VCF$ALT_AD) <- 'integer'
  
  # extract depth
  # VCF$DP <- str_split(IN,':',simplify = TRUE)[,3]
  # dim(VCF$DP) <- dim(IN)
  # storage.mode(VCF$DP) <- 'integer' # can produce NA when '.'
  
  gt_names <- VCF$SAMPLES[sort(c(1:length(VCF$SAMPLES), 1:length(VCF$SAMPLES)))]
  # because I need every genotye name twice, I extract it twice.
  # to extract it twice, I just numerate all names twice and then sort the indices
  
  phase <- data.table(gt_names, phase)
  colnames(phase) <- c("Genotype", VCF$MAP$ID)
  
  # add runtime message
  if(verbose){
    cat('The process needed ',format(Sys.time() - start),
        'and the object is',format(object.size(phase),'Mb'),
        'large.\n')
  }
  
  # return final object
  return(phase)
}
# read_vcf_to_phase("outfile.vcf")


# reads in a vcf file. Function from class8_solution, Applied effective programming 
# (Johannes class)
read_vcf <- function(file, verbose = TRUE){ # add verbose to allow for messages
  
  # record start time
  start <- Sys.time()
  
  # read in file
  IN <- fread(file,
              skip = '#CHROM',
              sep = '\t', sep2=NULL,
              na.strings = '.')
  
  # generate initial structure with map
  VCF <- list(MAP = IN[,1:8])
  
  # define function to structure your info field
  structure_info <- function(x){
    x <- str_split(x,';') # split info filed by elements
    
    # define a function which seperates name and content and produces named list
    name_list <- function(x){
      out <- as.list(as.numeric(str_split(x,'=',simplify = TRUE)[,2])) # extract content of info entry
      names(out) <- str_split(x,'=',simplify = TRUE)[,1] # extract name of info entry
      names(out) <- str_extract(x,'^[^=]+') # alternative using str_extract
      return(out)
    }
    
    # apply function to list elements - output is a list of length(x)
    lapply(x,name_list)
  }
  
  VCF$MAP[,INFO := structure_info(INFO)]
  VCF$MAP$INFO[[1]]$DP
  
  nheader_cols <- sum(!str_detect(IN[1,], "[/\\|]"))
  
  # remove unnecessary columns from IN and transform to matrix
  IN <- as.matrix(IN[,-1:-nheader_cols])
  
  # extract samples
  VCF$SAMPLES <- colnames(IN)
  
  # extract genotypes
  VCF$GT <- str_sub(IN,1,3)
  lookup <- c('0/0' = 0,
              '0|0' = 0,
              '0/1' = 1,
              '0|1' = 1,
              '1|0' = 1,
              '1/1' = 2,
              '1|1' = 2,
              './.' = NA,
              '.|.' = NA)
  VCF$GT <- unname(lookup[VCF$GT])
  dim(VCF$GT) <- dim(IN)
  
  ###### The following is for if you have more information per individual. Maybe 
  # build a function that checks what formats are present based on the FORMAT column
  # and then executes the respective commands
  
  # extract AD of alternative allele
  # VCF$ALT_AD <- str_split(IN,':',simplify = TRUE)[,2]
  # VCF$ALT_AD <- str_split(VCF$ALT_AD,',',simplify = TRUE)[,2]
  # dim(VCF$ALT_AD) <- dim(IN)
  # storage.mode(VCF$ALT_AD) <- 'integer'
  
  # extract depth
  # VCF$DP <- str_split(IN,':',simplify = TRUE)[,3]
  # dim(VCF$DP) <- dim(IN)
  # storage.mode(VCF$DP) <- 'integer' # can produce NA when '.'
  
  # assign class VCF to final object
  class(VCF) <- 'VCF'
  
  # add runtime message
  if(verbose){
    cat('The process needed ',format(Sys.time() - start),
        'and the object is',format(object.size(VCF),'Mb'),
        'large.\n')
  }
  
  # return final object
  return(VCF)
}

# res <- read_vcf('subsetchickenAppliedR.vcf.gz')
# str(res,max.level = 2)
