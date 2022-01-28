library("multcomp")

################ 
# calc functions

# calculates the accuracy per marker
# if either in the true or imputed marker table (or both) information is missing, 
# this genotype_marker combination will not be considered for calculation
# only pairwise complete observations will be used
calc_accuracy_pmarker_imputation_mt <- function(true_marker_table, imputed_marker_table, 
                                                verbose = TRUE) {
  acc <- NULL
  err <- 1/1000000
  true_marker_table[true_marker_table == 9] <- NA
  imputed_marker_table[imputed_marker_table == 9] <- NA
  
  for(i in 1:nrow(imputed_marker_table)) {
    acc[i] <- cor(c(unlist(imputed_marker_table[i, -1]), err), 
                  c(unlist(true_marker_table[i, -1]), err), method = "pearson",
                  use = "pairwise.complete.obs")
    if(verbose){cat('finished', i/nrow(imputed_marker_table), '\n')}
  }
  return(data.table(marker = true_marker_table[[1]], acc))
}
#calc_accuracy_pmarker_imputation_mt(true, imputed)

# this function calculates the correlation between two genotypes as in Gonen et al.
# and in AlphaPlantImpute and AlphaPlantImpute2
# be aware that if information is missing for the true marker_table no comparison
# can be made
calc_accuracy_ind_imputation_mt <- function(true_marker_table, imputed_marker_table) {
  acc <- NULL
  err <- 1/1000000
  true_marker_table[true_marker_table == 9] <- NA
  imputed_marker_table[imputed_marker_table == 9] <- NA
  for(i in 2:ncol(imputed_marker_table)) {
    acc[i-1] <- cor(c(unlist(imputed_marker_table[, ..i]), err), 
                    c(unlist(true_marker_table[, ..i]), err), 
                    method = "pearson", use = "pairwise.complete.obs")
  }
  return(data.table(Genotype = colnames(true_marker_table)[-1], acc))
}
# calc_accuracy_ind_imputation_mt(true, imputed)

# this function compares how many markers in the imputed table have the same information
# as in the true table and returns it as a fraction of all markers with information
# this is the same as in AlphaPlantImpute 'GenotypeCorrectFilled'
calc_correct_imp_ind_imputation_mt <- function(true_marker_table, imputed_marker_table) {
  perc_correct <- NULL
  true_marker_table[true_marker_table == 9] <- NA
  imputed_marker_table[imputed_marker_table == 9] <- NA
  na_tab <- colSums(0 < (is.na(true_marker_table[,-1]) + is.na(imputed_marker_table[,-1])))
  
  l <- nrow(imputed_marker_table)
  for(i in 2:ncol(imputed_marker_table)) {
    nNA <- na_tab[i-1]
    perc_correct[i-1] <- 
      sum(true_marker_table[, ..i] == imputed_marker_table[, ..i], na.rm = TRUE)/(l-nNA)
  }
  return(data.table(Genotype = colnames(true_marker_table)[-1], perc_correct))
}
#calc_correct_imp_ind_imputation_mt(true, imputed)

# this function does the same as the above just on a per marker level
# the output is a marker table reporting the fraction of correctly imputed
# individuals per marker
# it compares if the value of the true marker table and the imputed marker table 
# are identical. Before it does that, it sets all missing values (i.e. 9) to NA
# because they cannot be compared meaningfully and if either the true or the 
# imputed marker table have a missing value, this one is not considered in the 
# calculation of number of markers
calc_correct_imp_marker_imputation_mt <- function(true_marker_table, imputed_marker_table){
  true_marker_table[true_marker_table == 9] <- NA
  imputed_marker_table[imputed_marker_table == 9] <- NA
  # this gives the number of missing values per marker
  na_tab <- rowSums(0 < (is.na(true_marker_table[,-1]) + is.na(imputed_marker_table[,-1])))
  
  # minus 1 because the marker name is positioned on the first column
  l <- ncol(imputed_marker_table)-1
  anal <- true_marker_table == imputed_marker_table
  perc_correct <- 
    unlist(apply(anal[,-1], MARGIN = 1
                 , function(i){sum(i, na.rm = TRUE)}))/(l - na_tab)
  
  return(data.table(Genotype = true_marker_table[[1]], perc_correct))
}
#calc_correct_imp_marker_imputation_mt(true_marker_table, imputed_marker_table)

# this function calculates the fraction of markers that don't have missing information
calc_nonmissing_marker_ind_marker_table <- function(marker_table) {
  l <- nrow(marker_table)
  geno_yield <- unlist(apply(marker_table[,-1], MARGIN = 2, function(ind) {
    sum(ind != 9)
  }))
  geno_yield <- geno_yield/l
  return(data.table(Genotype = colnames(marker_table)[-1], geno_yield))
}
# calc_nonmissing_marker_ind_marker_table(marker_table)

# calculates the distance between the true and imputed sequence
# The distance is e.g. true(2,1) imputed(1,0) -> difference ((2-1), (1-0)) = (1,1) = 2
# the maximally possible distance is 2+1 = 3 because the imputation could be
# (0,0) or (0,2). This measure is invented by me and just for fun
# it is better to have imputed 1 than 0 if the truth is 2 but the calculation
# of correct prediction doesn't account for that
calc_difference_ind_imputation_mt <- function(true_marker_table, imputed_marker_table) {
  true_marker_table[true_marker_table == 9] <- NA
  imputed_marker_table[imputed_marker_table == 9] <- NA
  
  het <- unlist(apply(true_marker_table[,-1], MARGIN = 2, 
                      function(ind) {sum(ind == 1, na.rm = TRUE)}))
  hom <- unlist(apply(true_marker_table[,-1], MARGIN = 2, 
                      function(ind) {sum(ind == 0, na.rm = TRUE) + 
                          sum(ind == 2, na.rm = TRUE)}))
  max_poss_diff <- het + hom*2
  
  diff_abs <- NULL
  for(i in 2:ncol(imputed_marker_table)) {
    diff_abs[i-1] <- sum(abs(true_marker_table[[i]] - imputed_marker_table[[i]]), 
                         na.rm = TRUE)
  }
  difference <- (diff_abs/max_poss_diff)
  return(data.table(Genotype = colnames(imputed_marker_table)[-1], difference))
}
# calc_difference_ind_imputation_mt(true, imputed)

# this function calculates the distance, here called difference, between the true 
# and the imputed marker table per marker
calc_difference_marker_imputation_mt <- function(true_marker_table, imputed_marker_table){
  true_marker_table[true_marker_table == 9] <- NA
  imputed_marker_table[imputed_marker_table == 9] <- NA
  
  het <- 
    apply(true_marker_table[,-1], MARGIN = 1, function(i) {sum(i == 1, na.rm = TRUE)})
  hom0 <- 
    apply(true_marker_table[,-1], MARGIN = 1, function(i) {sum(i == 0, na.rm = TRUE)})
  hom2 <- 
    apply(true_marker_table[,-1], MARGIN = 1, function(i) {sum(i == 2, na.rm = TRUE)})
  #because you can only be one wrong if the truth is 1 but two steps wrong
  # if the truth is 0 or 2
  max_poss_diff <- het + hom0*2 + hom2*2
  diff <- abs(true_marker_table[,-1] - imputed_marker_table[,-1])
  diff_abs <- apply(diff, MARGIN = 1, function(i){sum(i, na.rm = TRUE)})
  difference <- diff_abs/max_poss_diff
  return(data.table(Genotype = imputed_marker_table[[1]], difference))
}
#calc_difference_marker_imputation_mt(true_marker_table, imputed_marker_table)

#################################################
################### analyze functions


# analyses marker in a marker table (minor allele count, minor allele frequency,
# freq of ungenotyped (NA) individuals, number of individuals that are not NA)
analyze_marker_marker_table <- function(marker_table, inds_exclude = NULL) {

  marker_table <- extract_individuals_marker_table(marker_table = marker_table
                                                   ,ID_inds = inds_exclude
                                                   , negate = TRUE)
  marker_table[marker_table == 9] <- NA
  occNA_in_marker <- unlist(
    apply(marker_table[,-1], MARGIN = 1, function(marker){sum(is.na(marker))}))

  # checks how often NA occurs per marker
  max_count <- ((ncol(marker_table)-1)-occNA_in_marker)*2 # calcs how often an allele could theoretically be present
  count <- rowSums(marker_table[,-1], na.rm = TRUE)
  count[count > (max_count/2)] <- max_count[count > (max_count/2)] - count[count > (max_count/2)] 

  # checks if an allele is present more often than half of possible and if so, 
  # returns the other side (so that you get minor allele count)
  MAC <- count
  MAF <- count/max_count
  nind_nonNA <- max_count/2
  freq_NA <- occNA_in_marker/(ncol(marker_table)-1)
  
  return(data.table(marker = marker_table[[1]], MAC, MAF, freq_NA, nind_nonNA))
}
#analyze_marker_marker_table(marker_table = mt, inds_exlcude = C('P2', 'P34'))



analyze_marker_imputation <- function(true_marker_table, imputed_marker_table
                                      , inds_exclude = NULL){
  IN <- analyze_marker_marker_table(true_marker_table, inds_exclude = inds_exclude)
  true_MAF <- IN[[3]]
  true_freq_NA <- IN[[4]]
  
  IN <- analyze_marker_marker_table(imputed_marker_table)
  imp_MAF <- IN[[3]]
  imp_freq_NA <- IN[[4]]
  
  acc <- calc_accuracy_pmarker_imputation_mt(true_marker_table = true_marker_table, 
                                             imputed_marker_table = imputed_marker_table,
                                             verbose = FALSE)[[2]]
  crrct <- calc_correct_imp_marker_imputation_mt(true_marker_table = true_marker_table, 
                                                 imputed_marker_table = imputed_marker_table)[[2]]
  diff <- calc_difference_marker_imputation_mt(true_marker_table = true_marker_table, 
                                               imputed_marker_table = imputed_marker_table)[[2]]
  
  out <- data.table(MarkerID = IN[[1]], MarkerAccuracy = acc,
                    MarkerCorrect = crrct, MarkerDifference = diff, 
                    MarkerMissingGenotypeTrue = true_freq_NA, 
                    MarkerMissingGenotypeImputed = imp_freq_NA,
                    MarkerMAFTrue = true_MAF, MarkerMAFImputed = imp_MAF)
  return(out)
}
# analyze_marker_marker_table(marker_table, inds_exlcude = 'P1')



# analyzes the individuals in a marker table
# gives frequency of NA markers (9 in the marker table indicates NA) and proportion
# of the genome being homozygous based on marker information
# if run with default settings, this function is just as fast as it's 
# predecessor (the one called 'analyze' within this function definition)
#
# This function offers more options. The parents can be excluded from the analysis.
# The analysis can be done for only the variable loci.
#
# choosing 'per_family' can result in confusing results. If set to TRUE, the 
# variable loci are found per family and the number of variabl eloci might differ
# between families
#
# could be extended to exclude specified individuals, not just parents
# difficulty: how to deal with multiple predictions?
analyze_inds_marker_table <- function(marker_table, pedigree = NULL, parents_exclude = FALSE, 
                                       only_variable_loci = FALSE, per_family = FALSE){
  
  if(parents_exclude){
    if(is.null(pedigree)){stop("Can't exclude parents if no pedigree is specified.")}
  }
  
  ##############################################
  # This is the actual function that does the analysis and everything around it 
  # is just to feed in right data
  # originally, only this part was the function for analyzing a marker table 
  # I did not make to separate functions in which one would call the other because
  # I want to keep the number of functions down especially when they are 
  # doing very similar things
  analyze <- function(marker_table) {
    nmarker <- nrow(marker_table)
    nNA_inds <- apply(marker_table[,-1], MARGIN = 2, 
                      function(ind) {sum(ind == 9, na.rm = TRUE)})
    het_inds <- apply(marker_table[,-1], MARGIN = 2, 
                      function(ind) {sum(ind == 1, na.rm = TRUE)})
    
    nNA_inds[is.na(nNA_inds)] <- 0
    freq_NA <- nNA_inds/nmarker
    
    prop_homozygous <- 1-(het_inds/(nmarker - nNA_inds))
    
    return(data.table(Genotype = colnames(marker_table)[-1], freq_NA, prop_homozygous))
  }
  ###############################################
  
  
  if(per_family & !is.null(pedigree)) {
    ls_mt <- extract_family_marker_table(pedigree = pedigree, marker_table = marker_table)
    
    out <- list()
    # I could also remove the parents before BUT
    # then no population called 'parents' can be analyzed
    for(i in 1:length(ls_mt)){
      # if TRUE, removes all loci that have no variation left (9 (=missing) is
      # is considered as variation)
      if(only_variable_loci) {
        ls_mt[[i]] <- 
          extract_variable_loci_marker_table(marker_table = ls_mt[[i]])
      }
      
      # this bit extracts the names of parents and removes them from the
      # families they are involved in. Otherwise it will be confusing later
      # Parents will still be analyzed if they are part of the 'parents' population
      # 
      # I could also loop over ls_mt and remove all the parents first
      # and then in a second loop do the actual analysis
      parents <- 
        extract_parents_pedigree(pedigree = pedigree, 
                                 ID_names = colnames(ls_mt[[i]]))[[1]]
      ls_mt[[i]] <- 
        extract_individuals_marker_table(marker_table = ls_mt[[i]], 
                                         ID_inds = parents, negate = TRUE)
      
      out <- c(
        out, list(analyze(marker_table = ls_mt[[i]]))
      )
    }
    
    out <- data.table::rbindlist(out) 
  } 
  # else, if no pedigree has been specified
  else{
    if(only_variable_loci) {
      marker_table <- 
        extract_variable_loci_marker_table(marker_table = marker_table)
    }
    
    out <- analyze(marker_table = marker_table)
  }
  
  if(parents_exclude){
    parents <- 
      extract_parents_pedigree(pedigree = pedigree)[[1]]
    out <- out[!(out[[1]] %in% parents),]
  }
  
  return(out)
}
# m <- matrix(c(1,2,2,9,2,0,0,1,1,0,0,0), ncol = 3, byrow = TRUE)
# m <- data.frame(marker_ID  = c("M1", "M2", "M3", "M4"), m)
# analyze_inds_marker_table(m)
# analyze_inds_marker_table(marker_table = m, only_variable_loci = TRUE)
# analyze_inds_marker_table(marker_table)

# a function that provides a very simple analysis of the pedigree
# it basically just counts the number of offspring every individual has
# # and reports this, the genotype name and the family name in a data table
analyze_inds_pedigree <- function(pedigree){
  dt <- pedigree[, 1:3]
  colnames(dt) <- c('Genotype', 'Pop', 'n_offspring')
  for(i in 1:nrow(pedigree)){
    ind <- pedigree[[1]][i]
    dt[[3]][i] <- nrow(extract_offspring_pedigree(pedigree = pedigree, parents_ID = ind))
  }
  return(dt)
}
#o <- analyze_inds_pedigree(pedigree)


# this function is to analyze the imputation
# description same as 'analyze_ind_marker_table'
# suggested to leave 'per_family' as FALSE
analyze_inds_imputation <- function(true_marker_table, imputed_marker_table,
                                     pedigree = NULL, parents_exclude = FALSE, 
                                     per_family = FALSE,
                                     only_variable_loci = FALSE){
  
  if(parents_exclude){
    if(is.null(pedigree)){stop("Can't exclude parents if no pedigree is specified.")}
  }
  
  ##############################################
  # This is the actual function that does the analysis and everything around it 
  # is just to feed in right data
  # originally, only this part was the function for analyzing a marker table 
  # I did not make to separate functions in which one would call the other because
  # I want to keep the number of functions down especially when they are 
  # doing very similar things
  analyze <- function(true_marker_table, imputed_marker_table) {
    fnmiss_marker_true <-   1-calc_nonmissing_marker_ind_marker_table(true_marker_table)[[2]]
    fnmiss_marker_imputed <-  1-calc_nonmissing_marker_ind_marker_table(imputed_marker_table)[[2]]
    
    acc <- calc_accuracy_ind_imputation_mt(true_marker_table, imputed_marker_table)
    crrct <- calc_correct_imp_ind_imputation_mt(true_marker_table, imputed_marker_table)[[2]]
    diff <- calc_difference_ind_imputation_mt(true_marker_table, imputed_marker_table)[[2]]
    
    out <- data.table(GenotypeID = acc[[1]], GenotypeAccuracy = acc[[2]], 
                      GenotypeCorrect = crrct, GenotypeDifference = diff, 
                      GenotypeMissingMarkerTrue = fnmiss_marker_true,
                      GenotypeMissingMarkerImputed = fnmiss_marker_imputed)
    return(out)
  }
  ###############################################
  
  
  if(per_family & !is.null(pedigree)) {
    ls_mt_true <- 
      extract_family_marker_table(pedigree = pedigree, marker_table = true_marker_table)
    ls_mt_imp <- 
      extract_family_marker_table(pedigree = pedigree, marker_table = imputed_marker_table)
    
    
    out <- list()
    # I could also remove the parents before BUT
    # then no population called 'parents' can be analyzed
    for(i in 1:length(ls_mt_true)){
      # if TRUE, removes all loci that have no variation left (9 (=missing) is
      # is considered as variation)
      if(only_variable_loci) {
        ls_mt_true[[i]] <- 
          extract_variable_loci_marker_table(marker_table = ls_mt_true[[i]])
        # it is possible that loci are fixed in the true mt that are not fixed 
        # in the imputed mt. Therefore, I filter only for the variable loci
        # in the true_mt and then extract these markers from the imputed
        # marker_table. Otherwise a comparison is not possible
        ls_mt_imp[[i]] <- 
          extract_marker_marker_table(marker_table = ls_mt_imp[[i]], 
                                      marker_IDs = ls_mt_true[[i]][[1]])
        
      }
      
      # this bit extracts the names of parents and removes them from the
      # families they are involved in. Otherwise it will be confusing later
      # Parents will still be analyzed if they are part of the 'parents' population
      # 
      # I could also loop over ls_mt and remove all the parents first
      # and then in a second loop do the actual analysis
      
      parents <- 
        extract_parents_pedigree(pedigree = pedigree, 
                                 ID_names = colnames(ls_mt_true[[i]]))[[1]]
      ls_mt_true[[i]] <- 
        extract_individuals_marker_table(marker_table = ls_mt_true[[i]], 
                                         ID_inds = parents, negate = TRUE)
      
      ls_mt_imp[[i]] <- 
        extract_individuals_marker_table(marker_table = ls_mt_imp[[i]], 
                                         ID_inds = parents, negate = TRUE)
      
      out <- c(
        out, list(analyze(true_marker_table = ls_mt_true[[i]],
                          imputed_marker_table = ls_mt_imp[[i]]))
      )
    }
    
    out <- data.table::rbindlist(out) 
  } 
  # else, if the analysis should not be done per family
  else{
    if(only_variable_loci) {
      true_marker_table <- 
        extract_variable_loci_marker_table(marker_table = true_marker_table)
      # extract the loci that are variable in the true marker table
      imputed_marker_table <- 
        extract_marker_marker_table(marker_table = imputed_marker_table, 
                                    marker_IDs = true_marker_table[[1]])
    }
    
    out <- analyze(true_marker_table = true_marker_table, 
                   imputed_marker_table = imputed_marker_table)
  }
  
  if(parents_exclude){
    parents <- 
      extract_parents_pedigree(pedigree = pedigree)[[1]]
    out <- out[!(out[[1]] %in% parents),]
  }
  
  return(out)
}
# analyze_inds_imputation(true, imputed)

# analyzes a whole population at once and returns averages, precision, and 
# number of individuals
# highly suggested to leave 'per_family' as 'FALSE' unless you really know what's
# going on
analyze_pop_imputation <- function(true_marker_table, imputed_marker_table,
                                    pedigree = NULL, parents_exclude = FALSE, 
                                    per_family = FALSE,
                                    only_variable_loci = FALSE) {
  IN <- 
    analyze_inds_imputation(true_marker_table = true_marker_table
                            , imputed_marker_table = imputed_marker_table
                            , pedigree = pedigree
                            , parents_exclude = parents_exclude
                            , per_family = per_family
                            , only_variable_loci = only_variable_loci)
  av <- colMeans(IN[,2:ncol(IN)], na.rm = TRUE)
  precision <- log(1/var(IN[[2]])) # calculated after Gonen et al. 2018 (Analysis)
  nGenotypes <- nrow(IN)
  out <- as.list(c(nGenotypes, precision, unlist(av)))
  names(out) <- c('nGenotypes', 'Precision', paste0("Mean", names(av)))
  
  return(out)
}
# analzye_pop_imputation(true, imputed)

# thsi function does the same as analyze_inds_imputation, just per marker
# additionally it also returns the allele frequency of the true and the imputed
# marker table. This might be useful if one wants to see whether an imputation method
# favor a particular allele
analyze_marker_imputation <- function(true_marker_table, imputed_marker_table
                                      , inds_exclude = NULL){
  IN <- analyze_marker_marker_table(true_marker_table, inds_exclude = inds_exclude)
  true_MAF <- IN[[3]]
  true_freq_NA <- IN[[4]]
  
  IN <- analyze_marker_marker_table(imputed_marker_table, inds_exclude = inds_exclude)
  imp_MAF <- IN[[3]]
  imp_freq_NA <- IN[[4]]
  
  acc <- calc_accuracy_pmarker_imputation_mt(true_marker_table = true_marker_table, 
                                             imputed_marker_table = imputed_marker_table,
                                             verbose = FALSE)[[2]]
  crrct <- calc_correct_imp_marker_imputation_mt(true_marker_table = true_marker_table, 
                                                 imputed_marker_table = imputed_marker_table)[[2]]
  diff <- calc_difference_marker_imputation_mt(true_marker_table = true_marker_table, 
                                               imputed_marker_table = imputed_marker_table)[[2]]
  
  out <- data.table(MarkerID = IN[[1]], MarkerAccuracy = acc,
                    MarkerCorrect = crrct, MarkerDifference = diff, 
                    MarkerMissingGenotypeTrue = true_freq_NA, 
                    MarkerMissingGenotypeImputed = imp_freq_NA,
                    MarkerMAFTrue = true_MAF, MarkerMAFImputed = imp_MAF)
  return(out)
}
#analyze_marker_imputation(true_mt, mt)

# this function does the same as analyze_pop_imputation just on a per marker
# basis. It runs analyze_marker_imputation and takes the means of its output and 
# calculates teh precision
analyze_marker_panel_imputation <- function(true_marker_table, imputed_marker_table
                                            , inds_exclude = NULL) {
  IN <- analyze_marker_imputation(true_marker_table = true_marker_table,
                                  imputed_marker_table = imputed_marker_table
                                  , inds_exclude = inds_exclude)
  av <- colMeans(IN[,2:ncol(IN)], na.rm = TRUE)
  precision <- log(1/var(IN[[2]])) # calculated after Gonen et al. 2018 (Analysis)
  # BUT note that Gonen et al calculate this per individual, not per marker
  nMarkers <- nrow(IN)
  out <- as.list(c(nMarkers, precision, unlist(av)))
  names(out) <- c('nMarkers', 'Precision', paste0("Mean", names(av)))
  return(out)
}
#analyze_marker_panel_imputation(true_marker_table = true_mt, imputed_marker_table = mt)


# loops through pedigree and check for every individual
# family identifier is not needed
# so far, this function only works well if the offspring are no direct offspring
# if they are, more information can be said
#
# the output will be a marker table like object that shows for every individual 
# at every marker whether it is consistent with the parent markers (TRUE) or not (FALSE)
# if it is not possible to check, say because one of the parents or offspring is 9
# or because the parents are 0 and 2 or 1 and 1, NA is given
# this function is used by others
analyze_mendelian_consistency <- function(marker_table, pedigree) {
  out_mt <- data.table(matrix(NA, nrow = nrow(marker_table), ncol = ncol(marker_table)))
  out_mt[[1]] <- marker_table[[1]]
  colnames(out_mt) <- colnames(marker_table)
  for(i in 1:nrow(pedigree)) {
    parents <- unlist(pedigree[i, 3:4])
    
    # execute the following only if both parents are known
    if(!any(parents == 0)){
      #reduce marker_table to only the parents
      a <- marker_table[,..parents]
      # find column of individual in question
      ind <- match(pedigree[i,][[1]], colnames(marker_table))
      # find at which positions the parents are homozygous for the same allele and offspring is not 9
      pos <- ((a[[1]] == a[[2]]) & (a[[1]] %in% c(0,2))) & marker_table[[ind]] != 9
      # check if ind is equal with parents at loci for which they are homozygous for the same allele
      out_mt[[ind]][pos] <- a[[1]][pos] == marker_table[,..ind][pos]## mistake
    }
  }
  return(out_mt)
}
#o <- analyze_mendelian_consistency(marker_table = marker_table, pedigree = pedigree)
#o[19500:19519, c('P17','P35','G48', 'G49')] # it is known that G48 is incorrent at marker C9M1953 in the true data set

# this functions uses the above one
# this function will tell you how many mistakes there are per individual
# and automatically prints a histogram
#
# Update 20.02.2021: I now give additional output stating the fraction of 
# inconsistencies per all positions at which it was possible to make a comparison.
# Say both parents are 1, it is impossible to check for mendelian inconsistency here.
analyze_mendel_nerr_ind <- function(marker_table, pedigree) {
  o <- analyze_mendelian_consistency(marker_table = marker_table, pedigree = pedigree)
  nerrors <- apply(o[,-1], MARGIN = 2, function(i){sum(!i, na.rm = TRUE)})
  
  # now checking at how many positions I could actually make a check
  w <- lapply(o[,-1], is.na)
  npos_check <- nrow(marker_table) - unlist(lapply(w, sum))
  ratio_nerror_npos.checkable <- nerrors/npos_check
  
  nerrors <- data.table(Genotype = names(nerrors), nerrors, ratio_nerror_npos.checkable)
  hist(nerrors[[2]], ylim = c(0,40), breaks = 40, main = "Mendelian errors per ind", xlab = 'errors per ind')
  return(nerrors)
}
# w <- analyze_mendel_nerr_ind(marker_table, pedigree)
# sum(w[[2]] > 0)/length(w[[2]]) # reports what fraction of individuals have at least one error

# this function does the same as the above just per marker
# it will return a data table that shows how many mistakes there are per marker
#
# Update 20.02.2021: now the ratio of errors per individual at which a comparison
# could be made is included
analyze_mendel_nerr_marker <- function(marker_table, pedigree) {
  o <- analyze_mendelian_consistency(marker_table = marker_table, pedigree = pedigree)
  nerrors <- apply(o[,-1], MARGIN = 1, function(i){sum(!i, na.rm = TRUE)})
  
  # now checking at how many individuals I could actually make a check
  # I subtract the number of parents because the way I set it, a parent
  # cannot be declared having a mendelian inconsistency
  w <- lapply(as.data.frame(t(o[,-1])), is.na)
  nparents <- nrow(extract_parents_pedigree(pedigree))
  ninds_check <- (ncol(marker_table)-nparents) - (unlist(lapply(w, sum))-nparents)
  ratio_nerror_ninds.checkable <- nerrors/ninds_check
  
  f <- data.table(Marker = o[[1]], nerrors, ratio_nerror_ninds.checkable)
  hist(f[[2]], ylim = c(0,20), xlim = c(0,max(f[[2]])+10), main = 'Mendelian errors per marker'
       , xlab = "n errors per markers")
  return(f)
}
# e <- analyze_mendel_nerr_marker(marker_table, pedigree)
# sum(e[[2]] > 0)/length(e[[2]]) # this shows the fraction of markers that have at least one error
