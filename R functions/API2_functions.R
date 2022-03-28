source("general_functions.R")
#######################
# make functions

# this function makes a founders object as required b y API2 and to be used by
# 'write_founders_file'
make_founders_obj <- function(pedigree) {
  out <- sapply(pedigree[[1]], function(i) {extract_ancestors_pedigree(pedigree, ID_names = i)})
  out <- data.frame(lapply(out, "length<-", max(lengths(out))))
  out <- as.data.table(cbind(colnames(out), t(out)))
  
  return(out)
}
# make_founders_obj(pedigree)


# this function creates the arguments that can then be directly passed to the 
# command line
make_arguments_impute_API2 <- function(out, libphase, genotypes, founders, seed, 
                                       n_haplotypes = 100, n_windows = 5, 
                                       overwrite = FALSE, write_args_to_txt = FALSE,
                                       arg_filename = "arguments_impute_API2.txt",
                                       ...) {
  if(missing(out)) stop("You need to provide an output file prefix for the 'out'
                         argument.")
  if(missing(libphase)) stop("You need to provide the name of the haplotype 
                              library for the 'libphase' option 
                              (.phase or .phase).")
  if(missing(genotypes)) stop("You need to provide the name of the genotypes file
                               for which you want to perform imputation 
                               (.genotypes or .txt).")
  if(missing(seed)) {
    seed <- round(runif(1, 1, 1000000))
    cat(paste0("A seed was randomly set as ", seed,".\n"))
  }
  
  arg <- c(as.list(environment()), list(...))
  arg <- arg[!names(arg) %in% c("write_args_to_txt", "arg_filename")]
  arg <- arg[as.logical(str_length(as.character(arg)))]
  # I filter out the arguments that do not have a value. A positive number is TRUE
  
  ifelse(overwrite, (arg$overwrite <- ""), arg <- arg[!(names(arg) == "overwrite")])
  
  arg <- as.list(paste0("-", names(arg), " ", arg))
  if(write_args_to_txt){
    fwrite(as.list(format(Sys.time())), arg_filename, sep = "\n")
    fwrite(as.list("-impute"), arg_filename, sep = "\n", append = TRUE)
    fwrite(arg, arg_filename, sep = "\n", append = TRUE)
    }
  arg <- paste0("-impute ", paste(arg, collapse = " "))
  return(arg)
}
#make_arguments_impute_API2(out = "imputed", libphase = "lib500.phase", genotypes = "masked.genotypes")


# makes the arguments that can be passed to the command line to create
# a haplotype library with API2
# 'founders' can only be specified for the '-impute' option but not for
# the '-createlib' option
make_arguments_createlib_API2 <- function(out, genotypes, hd_threshold = 0,
                                          n_sample_rounds = 5, seed,
                                          n_haplotypes = 100, n_windows = 5,
                                          overwrite = FALSE, write_args_to_txt = FALSE,
                                          arg_filename = "arguments_library_creation_API2.txt",
                                          ...) {
  if(missing(out)) stop("You need to provide an output file prefix for the 'out'
                         argument.")
  if(missing(genotypes)) stop("You need to provide the name of the genotypes file
                               for which you want to perform imputation 
                               (.genotypes or .txt). Use 
                              convert_marker_table_to_GenotypeFile().")
  if(missing(seed)) {
    seed <- round(runif(1, 1, 1000000))
    cat(paste0("A seed was randomly choosen as ", seed,".\n"))
  }
  
  arg <- c(as.list(environment()), list(...))
  arg <- arg[!names(arg) %in% c("write_args_to_txt", "arg_filename")]
  arg <- arg[as.logical(stringr::str_length(as.character(arg)))]
  # I filter out the arguments that do not have a value. A positive number is TRUE
  
  ifelse(overwrite, (arg$overwrite <- ""), arg <- arg[!(names(arg) == "overwrite")])
  
  arg <- as.list(paste0("-", names(arg), " ", arg))
  if(write_args_to_txt){
    data.table::fwrite(as.list(format(Sys.time())), arg_filename, sep = "\n")
    data.table::fwrite(as.list("-createlib"), arg_filename, sep = "\n", append = TRUE)
    data.table::fwrite(arg, arg_filename, sep = "\n", append = TRUE)
    }
  arg <- paste0("-createlib ", paste(arg, collapse = " "))
  return(arg)
}
#make_arguments_createlib_API2(out = "lib500", genotypes = "masked.genotypes", 
#hd_threshold = 0.9, n_sample_rounds = 5)

####################################
# write functions


# this function takes a founders object and writes it to a founder file so that 
# it can be used by API2
write_founders_file <- function(founders_obj, file_name = "founders.txt") {
  data.table::fwrite(founders_obj, file_name, sep = " ", quote = FALSE, col.names = FALSE)
}

# writes a GenotypeFile
write_GenotypesFile <- function(GenotypesFile, filename) {
  data.table::fwrite(GenotypesFile, quote = FALSE, sep = " ", file = filename, col.names = FALSE,
         row.names = FALSE)
}
# gtf <- convert_marker_table_to_GenotypeFile(mt)
# write_GenotypesFile(gtf, "maskedGenotypes.txt")

############
# check functions
#checks if the python version is at least a value specified by user
# default is "3.7"
check_version_python <- function(at_least_vers = "3.7") {
  # check if python version is at least 3.7
  cmd_out <- system2(command = "python3", args = "--version", stdout = TRUE)
  cmd_out <- stringr::str_extract(cmd_out[1], pattern = "[[:digit:]].*[[:digit:]]")
  
  is_sufficient <- cmd_out >= at_least_vers
  return(is_sufficient)
}
#check_version_python()


####################
# impute functions
# just supply the arguments and this function passed them to the cmd
impute_single_chrom_API2 <- function(arguments_impute) {
  if(!check_version_python(at_least_vers = "3.7")) stop("Python version is not 
                                                       sufficient. Please use 
                                                       python 3.7")
  begin <- Sys.time()
  system2(command = "AlphaPlantImpute2", 
          args = arguments_impute,
          stdout = TRUE)
  time_needed <- difftime(Sys.time(), begin, units = "mins")
  return(time_needed)
}

#impute_API2(arguments_impute)

# kinda an impute function

# the execution of createlib and impute are exactly the same. only other arguments
# are passed to the cmd
# to not confuse the user, I created this function but they can be used interchangeably 
# (sind austauschbar)
createlib_single_chrom_API2 <- function(arguments_createlib){
  run_time <- impute_single_chrom_API2(arguments_impute = arguments_createlib)
  return(run_time)
}
#createlib_single_chrom_API2(arguments_createlib = arg_createlib)


# this function is very similar in it's behavior to impute_API
# it is a wrapper for impute_single_chrom_API2 that loops over the chromosomes
impute_API2 <- function(arguments_impute, genotypes, founders
                            , GeneticMapFile, libphase, verbose = TRUE) {
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start), "\n")}
  
  original_dir <- getwd()
  
  if(missing(genotypes)) stop("You need to provide the name of the
                                      GenotypeFile object.")
  if(missing(libphase)) stop("You need to provide the name of the
                                      library object.")
  if(ncol(libphase) != ncol(genotypes)) stop("The genotypes object and the library 
                                               object have an unequal number 
                                               of columns,")
  
  dir.create("temporary")
  setwd("./temporary")
  
  # extracts only the markers of the mapfile that are present in the TrueGenotypesFile
  GeneticMapFile <- GeneticMapFile[GeneticMapFile[[2]] %in% colnames(genotypes),]
  genotypes <- extract_marker_phase(genotypes, marker_IDs = GeneticMapFile[[2]])
  genotypes <- order_phase_by_mapo(phase = genotypes, mapo = GeneticMapFile[[2]])
  
  param <- stringr::str_split(arguments_impute, pattern = " ")[[1]]
  # this file does not need to be rewritten every time
  if(any(param == "-founders")) {
    name_founders <- param[which(param == "-founders")+1]
    write_founders_file(founders, name_founders)
  }
  
  tab <- table(GeneticMapFile[[1]])
  s <- 1
  # looping over all chromosomes and imputing them
  for(i in 1:length(names(tab))) {
    # check which columns need to be extracted
    ncol <- seq(s, sum(tab[1:i]))
    s <- sum(tab[1:i]) + 1
    
    GeneticMapFile_chrom <- GeneticMapFile[ncol,] 
    
    ncol <- c(1, ncol + 1)
    
    # extract the markers belonging to the respective chromosome
    GenotypeFile_chrom <- genotypes[,..ncol]
    libphase_chrom <- libphase[,..ncol]
    
    name_GenotypesFile_chrom <- param[which(param == "-genotypes")+1]
    write_GenotypesFile(GenotypesFile = GenotypeFile_chrom
                        , filename = name_GenotypesFile_chrom)
    
    name_libphase <- param[which(param == "-libphase")+1]
    
    write_GenotypesFile(GenotypesFile = libphase_chrom
                        , filename = name_libphase)
    if(i == 1){
      run_time <- impute_single_chrom_API2(arguments_impute = arguments_impute)
    } else {
      run_time <- run_time + impute_single_chrom_API2(arguments_impute = arguments_impute)
    }

    
    name_out <- paste0(param[which(param == "-out")+1], ".genotypes")
    imputed_genotypes_chrom <- data.table::fread(name_out)
    
    if(i == 1) {imputed_genotypes <- imputed_genotypes_chrom} 
    else {imputed_genotypes <- cbind(imputed_genotypes, imputed_genotypes_chrom[,-1])}
    
    if(verbose){print(paste0("Finished imputing chromosome ", names(tab)[i],
                             " of ", length(names(tab)), " chromosomes."))}
  }
  colnames(imputed_genotypes) <- colnames(genotypes)
  
  # this is done to include the parents in the output so that the input and 
  # output datasets are identical in dimensions
  missing_g <- unlist(genotypes[!(genotypes[[1]] %in% imputed_genotypes[[1]]),1])
  missing_g <- extract_individuals_phase(phase = libphase, ID_inds = missing_g)
  
  # sometimes the parents are already included in 'imputed_genotypes' so 
  # 'missing_g' would be empty then which causes problems
  if(nrow(missing_g) > 0){
    missing_g <- convert_marker_table_to_GenotypeFile(convert_phase_to_marker_table(missing_g))
  }  
  
  colnames(missing_g) <- colnames(genotypes)
  out <- convert_GenotypeFile_to_marker_table(rbind(missing_g, imputed_genotypes))
  out <- order_marker_table_by_pedigree(marker_table = out, pedigree = genotypes[[1]])

  #### after reading in the results, the subdirectory and its files can be deleted
  setwd("..")
  unlink("./temporary", recursive = TRUE)
  
  if(verbose){cat('The process needed ', format(Sys.time() - start))}
  
  return(list(imputed_genotypes = out, run_time = run_time))
}
#impute_API2(arguments_impute = arguments_impute, genotypes = genotypes,founders = founders_obj, GeneticMapFile = gmap, libphase = phase)

# this function was developed after I found out that API2 should only be run on 
# a per chromosome basis
# this function is very similar to the 'impute_API' function of API.
# the createlib option does no accept a founders file
#
# Update 13.02.2021: in advanced test it was noticed that there are issues with 
# the way this function is putting the haplotype libraries of every chromosome
# together. This used to be done with a cbind() command.
# However, since the arguments for API2 also include 'hd_threshold'
# it is possible that some individuals are below this threshold for one chromosome
# but above it for another chromosome. This results in libraries that contain 
# different individuals. And when fused together with cbind(), they don't fit 
# because the libraries have different sizes.
# The new way of merging the libraries together is accounting for that.
# I could not use the merge function on the two libraries directly because
# every individual has two rows which means that when merged together, I will 
# get four rows for an individual that is present in the old library and the 
# chromosome library. But I only want to have two rows per individual. 
# Technically, thw two rows are not the ground truth in a phase object
# because who is to say that the first chromsome of the first pair has to be 
# in front of the first chromosome of the secodn chromsome pair?
# They are not linked. But since this is taken into account whenever I deal
# with a phase object, I does not matter.
# Also, I have tested if API2 can handle libaries that have missing values and 
# it is fine. 
# This new way coded here is better than any other option because only the 
# chromosomes of a population that are sufficiently HD are used for referencing.
createlib_API2 <- function(arguments_createlib, genotypes
                           , GeneticMapFile, verbose = TRUE) {
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start), "\n")}
  
  original_dir <- getwd()
  
  if(missing(genotypes)) stop("You need to provide the name of the
                                      GenotypeFile object.")
  dir.create("temporary")
  setwd("./temporary")
  
  # extracts only the markers of the mapfile that are present in the TrueGenotypesFile
  GeneticMapFile <- GeneticMapFile[GeneticMapFile[[2]] %in% colnames(genotypes),]
  genotypes <- extract_marker_phase(genotypes, marker_IDs = GeneticMapFile[[2]])
  
  genotypes <- order_phase_by_mapo(phase = genotypes, mapo = GeneticMapFile[[2]])
  
  tab <- table(GeneticMapFile[[1]])
  s <- 1
  # looping over all chromosomes and imputing them
  for(i in 1:length(names(tab))) {
    # check which columns need to be extracted
    ncol <- seq(s, sum(tab[1:i]))
    s <- sum(tab[1:i]) + 1
    
    GeneticMapFile_chrom <- GeneticMapFile[ncol,] # comment 13.02.2021: this is 
    # unnecessarily complicated with
    # s and ncol
    
    ncol <- c(1, ncol + 1)
    
    # extract the markers belonging to the respective chromosome
    GenotypeFile_chrom <- genotypes[,..ncol]
    
    param <- stringr::str_split(arguments_createlib, pattern = " ")[[1]]
    name_GenotypesFile_chrom <- param[which(param == "-genotypes")+1]
    write_GenotypesFile(GenotypesFile = GenotypeFile_chrom, filename = name_GenotypesFile_chrom)
    
    if(i == 1){
      run_time <- createlib_single_chrom_API2(arguments_createlib = arguments_createlib)
    } else {
      run_time <- run_time + createlib_single_chrom_API2(arguments_createlib = arguments_createlib)
    }
    
    name_out <- paste0(param[which(param == "-out")+1], ".phase")
    library_phase_chrom <- data.table::fread(name_out)
    colnames(library_phase_chrom) <- c('Genotype', GeneticMapFile_chrom[[2]])
    
    if(i == 1) {library_phase <- library_phase_chrom} 
    else {
      # this here was changed on 13.02.2021
      lib_old_1 <- library_phase[!duplicated(library_phase[[1]]),]
      lib_chrom_1 <- library_phase_chrom[!duplicated(library_phase_chrom[[1]]),]
      lib_both_1 <- data.table::merge.data.table(lib_old_1, lib_chrom_1, all = TRUE, sort = FALSE)
      
      lib_old_2 <- library_phase[duplicated(library_phase[[1]]),]
      lib_chrom_2 <- library_phase_chrom[duplicated(library_phase_chrom[[1]]),]
      lib_both_2 <- data.table::merge.data.table(lib_old_2, lib_chrom_2, all = TRUE, sort = FALSE)
      
      lib_merged <- data.table::merge.data.table(lib_both_1, lib_both_2, all = TRUE, sort = FALSE)
      
      lib_merged <- order_phase_by_pedigree(phase = lib_merged, pedigree = genotypes[[1]])
      lib_merged[,2:ncol(lib_merged)][is.na(lib_merged[,2:ncol(lib_merged)])] <- 9
      library_phase <- lib_merged
      rm(lib_merged)
    }
    
    if(verbose) {print(paste0("Finished creating library for chromosome ", names(tab)[i],
                              " of ", length(names(tab)), " chromosomes."))}
  }
  colnames(library_phase) <- colnames(genotypes)
  
  #### after reading in the results, the subdirectory and its files can be deleted
  setwd("..")
  unlink("./temporary", recursive = TRUE)
  
  if(verbose){cat('The process needed ', format(Sys.time() - start),'./n')}
  
  return(list(library_phase = library_phase, run_time = run_time))
}
#createlib_API2(arguments_createlib = arguments_library, genotypes = genotypes, GeneticMapFile = gmap, verbose = TRUE)