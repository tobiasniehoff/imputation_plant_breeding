source("general_functions.R")

######################################
#################### make functions

# this function takes BEAGLE arguments as input parameters and converts them
# to one large space separated string that can be passed to BEAGLE on the cmd
# if write_args_to_txt is set to "TRUE", a txt file is written that contains 
# the arguents for BEAGLE line separated (each argument on a new line) for 
# better readablility
# the arguments stored like this can be read in an converted to BEAGLE accepted style
# with read_arguments
make_arguments_BEAGLE <- function(gt, out, write_args_to_txt = FALSE,
                                  arg_filename = "arguments_BEAGLE.txt",
                                  ref, map, seed = -99999,
                                  burnin = 6, iterations = 12,
                                  phase_states = 280, impute = "true",
                                  imp_states = 1600, imp_segment = 6.0,
                                  imp_step = 0.1, imp_nsteps = 7,
                                  cluster = 0.005, ap = "false",
                                  gp = "false", ne = 1000000,
                                  window = 40.0, overlap = 4.0, ...){
  arg <- c(as.list(environment()), list(...))
  arg <- arg[!names(arg) %in% c("write_args_to_txt", "arg_filename")]
  arg <- arg[as.logical(str_length(as.character(arg)))]
  names(arg) <- str_replace_all(names(arg), "\\_", "-")
  
  arg <- as.list(paste0(names(arg), "=", arg))
  
  if(write_args_to_txt){
    fwrite(as.list(format(Sys.time())), arg_filename, sep = "\n")
    fwrite(arg, arg_filename, sep = "\n", append = TRUE)
    # this writes the arguments_BEAGLE.txt file with the arguments in separate rows 
    # for better visualization
  }
  
  for_BEAGLE <- paste(arg, collapse = " ")
  return(for_BEAGLE)
}
# make_arguments_BEAGLE(gt = 3, out = "test.gt", write_args_to_txt = TRUE)

################################
# check functions
# checks if the java version is at least a value specified by user
# default is "1.8" required by BEAGLE 5.1
check_version_java <- function(at_least_vers = "1.8") {
  cmd_out <- system2(command = "java", args = "-version", stdout = TRUE
                     , stderr = TRUE)
  cmd_out <- stringr::str_extract(cmd_out[1], pattern = "[[:digit:]].*[[:digit:]]")
  
  is_sufficient <- cmd_out >= at_least_vers
  return(is_sufficient)
}
# check_version_java()


###############################
# Impute functions
# the files that should be used as input must be present in the same directory as 
# this script
# arguments <- read_arguments("arguments_BEAGLE.txt")
# this function is now a bit different from it's input requirements
# compared to API and API2
impute_BEAGLE <- function(arguments, memory_limit_GB, verbose = TRUE) {
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start), "\n")}
  
  original_dir <- getwd()
  os <- Sys.info()["sysname"]
  if(!str_detect(os, pattern = "[[inux][indows]]")) { 
    #checks if the operating system is not windows or Linux
    stop("Operating system is neither Windows nor Linux. Function was written for
         Windows and Linux only.")
  }
  
  # check if java version is at least 1.8
  if(!check_version_java(at_least_vers = "1.8")) {
    stop("Java version is lower than 1.8. Please install at least java version
         1.8 to proceed.")
  }
  
  # checks if the BEAGLE program is present in the current directory or any of 
  # its subdirectories.If not, the loop goes to the parent directory and checks 
  # again if beagle is present in it and any of its subdirectories. If not,
  # it steps to the grandparent directory and searches again everywhere.
  # If BEAGLE still cannot be found, an error message is issued
  # If BEAGLE is found, the program is copied to the working directory of this 
  # script.
  # Note: if the directory tree doesn't change and many analyses are conducted it 
  # is recommended to put BEAGLE in the directory from which this script is executed
  # so that it does not need to search for beagle every time
  for(i in 1:3) {
    lf <- list.files(recursive = TRUE, full.names = TRUE)
    beagle_file <- lf[str_detect(lf, pattern = "/beagle.18May20.d20.jar$")]
    if(length(beagle_file) == 0) {setwd("..")} 
    else {break}
  }
  rm(lf)
  if(length(beagle_file) == 0) {
    stop("The BEAGLE software should be present in the directory of this script or 
        in a subdirectory of two parent directories backwards.
        The name should be 'beagle.18May20.d20.jar'.")
  }
  file.copy(beagle_file, original_dir) # copies beagle to the original directory
  setwd(original_dir)
  
  # for limiting the amount of memory used
  limit <- NULL
  if(!is.na(memory_limit_GB)) {
    limit <- paste0("-Xmx", memory_limit_GB, "g") 
  }
  
  begin <- Sys.time()
  system2(command = "java", 
          args = paste(limit, "-jar", "beagle.18May20.d20.jar", arguments), 
          stdout = TRUE)
  time_needed <- difftime(Sys.time(), begin, units = "mins")
  
  file.remove("beagle.18May20.d20.jar")
  
  if(verbose){cat('The process needed', format(Sys.time() - start), "\n")}
  return(time_needed)
}
#arguments <- make_arguments_BEAGLE(gt = "gt_test_family05.vcf", out = "out.gt", write_args_to_txt = F)
# impute_BEAGLE(arguments = arguments, memory_limit_GB = 2)

# This function does familywise imputation. Therefore, it needs a pedigree
# to determine the families. It uses the name of the family in 'pedigree'
# to determine which individuals belong to what family
# the parents of a family are also included in the imputation
# a reference in phase format can also be provided. This function will then only
# take those individuals out of the reference phase that are in that particular 
# family that is to be imputed (including parents of course but also other phased
# offspring).
# This function returns a list of lists of phase objects, each belonging to one
# family. 
# Theoretically, if an individual is parent in two crosses and no phase 
# information of it is provided and imputation and phasing is bad, that can 
# result in two different outcomes for that parent
#
# handlingwise is this function more like 'impute_API' or 'impute_API2' by that
# you don't have to produce the input vcf files manually but just provide the 
# marker table and the function does everything by itself
#
# if you would like to impute a whole population at once with this function,
# make sure that you modify the pedigree object in such a way that all individuals 
# come from the same family
#
# NOTE: it seems as if BEAGLE wants to have an integer number for the memory limitation
#
# this function returns a list of phase objects for each family
#
#
# Update 19.02.2021: the argument 'restrict_ref_to_fam' was added with default 
# 'FALSE'. Previously, this function only used the individuals  of the 
# provided reference that are also present in the provided marker table.
# Now, months after the development of this function, this behaviour is 
# not wanted anymore as I may want to provide references of individuals that 
# are not present in the marker table. Code in the scripts does not have to 
# be changed
# However, filtering out individuals that are not present in the marker table
# is possible by setting 'restrict_ref_to_fam' to 'TRUE'
impute_familywise_BEAGLE <- function(arguments, marker_table, pedigree,
                                     memory_limit_GB = 2, reference_phase, mapo, 
                                     verbose = TRUE, restrict_ref_to_fam = FALSE) {######### update here the parameters
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start), "\n")}
  
  marker_table <- order_marker_table_by_pedigree(pedigree = pedigree, marker_table = marker_table)
  marker_table <- order_marker_table_by_mapo(mapo = mapo, marker_table = marker_table)
  ls_fam <- extract_family_marker_table(pedigree = pedigree, marker_table = marker_table)
  
  ##### create a subdirectory for this analysis
  dir.create("temporary")
  setwd("./temporary")
  ######################
  # the following is done to copy beagle in the current directory
  # for an explanation, look at 'impute_BEGLE'
  # this function would also work without these lines but it seems that it is slower 
  # then
  # I create two "temporary" files, one within the other. I do this so that I can 
  # save BEAGLE in one of them and run in the subdirectory. The problem is,
  # that 'impute_Beagle' will delete the beagle jar file if it is in its directory
  # after execution
  original_dir <- getwd()
  
  for(i in 1:3) {
    lf <- list.files(recursive = TRUE, full.names = TRUE)
    beagle_file <- lf[str_detect(lf, pattern = "/beagle.18May20.d20.jar$")]
    if(length(beagle_file) == 0) {setwd("..")} 
    else {break}
  }
  rm(lf)
  if(length(beagle_file) == 0) {
    stop("The BEAGLE software should be present in the directory of this script or 
        in a subdirectory of two parent directories backwards.
        The name should be 'beagle.18May20.d20.jar'.")
  }
  file.copy(beagle_file, original_dir) # copies beagle to the original directory
  
  setwd(original_dir)
  dir.create("temporary2")
  setwd("./temporary2")
  ########################################
  
  param <- stringr::str_split(arguments, pattern = " ")[[1]]
  name_map <- 
    stringr::str_remove(param[stringr::str_detect(param, pattern = "^map=")], "^map=")
  name_map <- stringr::str_remove(name_map, pattern = ".map$")
  map_PLINK <- convert_mapo_to_PLINK_map(mapo = mapo)
  write_PLINK_map(PLINK_map_object = map_PLINK, file_prefix = name_map) 
  
  use_ref <- any(stringr::str_detect(param, pattern = "^ref="))
  for(i in 1:length(ls_fam)) {
    if(verbose){
      start2 <- Sys.time()
      cat('Started imputing family', names(ls_fam)[[i]], 'of', length(ls_fam), 'families.\n')
      cat("The subprocess started at", format(start2), "\n")
    }
    mt_fam <- ls_fam[[i]]
    
    if(use_ref) {
      if(missing(reference_phase)) stop("You specified a reference in the arguments 
                                      but did not provide the reference phase object.")
      else {
        ref_phase_fam <- order_phase_by_mapo(mapo = mapo, phase = reference_phase)
        # here only the individuals of the family that are in the reference phase 
        # will be given to BEAGLE as a reference
        if(restrict_ref_to_fam){ # Updated 19.02. to allow to give references that are not in the family
          # this allows the code in the script for automatic imputation
          # to stay the same while still keeping the functionality
          # in case it it wanted
        ref_phase_fam <- 
          extract_individuals_phase(phase = ref_phase_fam
                                    , ID_inds = ref_phase_fam[ref_phase_fam[[1]] %in% colnames(mt_fam),1][[1]])
        }
        
        ref_name <- param[stringr::str_detect(param, pattern = "^ref=")]
        ref_name = stringr::str_remove(ref_name, "^ref=")
        as_gz <- stringr::str_detect(ref_name, pattern = '.gz$') 
        if(as_gz){ref_name <- stringr::str_remove(ref_name, pattern = ".gz$")}
        write_phase_to_vcf(mapo = mapo, phase = ref_phase_fam,
                           outfile_name = ref_name, gzipped = as_gz, verbose = FALSE)
      }
    }
    
    gt_name <- param[stringr::str_detect(param, "^gt=")]
    gt_name <- stringr::str_remove(gt_name, pattern = "^gt=")
    as_gz <- stringr::str_detect(gt_name, pattern = '.gz$')
    if(as_gz){gt_name <- stringr::str_remove(gt_name, pattern = ".gz$")}
    
    write_marker_table_to_vcf(marker_table = mt_fam, mapo = mapo, 
                              outfile_name = gt_name, gzipped = as_gz, 
                              verbose = FALSE)
    if(i == 1){
      run_time <- impute_BEAGLE(arguments = arguments
                                , memory_limit_GB = memory_limit_GB
                                , verbose = FALSE)
    }
    run_time <- run_time + impute_BEAGLE(arguments = arguments
                                         , memory_limit_GB = memory_limit_GB
                                         , verbose = FALSE)
    
    out_name <- param[stringr::str_detect(param, "^out=")]
    out_name <- paste0(stringr::str_remove(out_name, pattern = "^out="), ".vcf.gz")
    
    imputed_phase <- read_vcf_to_phase(file = out_name, set_missing_gt_to = 9, verbose = FALSE)
    
    if(i == 1) {
      out_phase <- list(imputed_phase)
    } else{
      out_phase <- c(out_phase, list(imputed_phase))
    }
    if(verbose){cat('The subprocess needed', format(Sys.time() - start2), "\n")}
  }
  
  #### after reading in the results, the subdirectory and its files can be deleted
  setwd("..")
  unlink("./temporary2", recursive = TRUE)
  setwd("..")
  unlink("./temporary", recursive = TRUE)
  
  names(out_phase) <- names(ls_fam)
  if(verbose){cat('The process needed', format(Sys.time() - start), "\n")}
  return(list(out_phase = out_phase, run_time = run_time))
}

#p <- impute_familywise_BEAGLE(arguments, marker_table = mt, pedigree = ped, mapo = mapo, memory_limit_GB = 2, verbose = TRUE, reference_phase = phase)

################################################################################
# Update 20.01.2022
# below, functions are defined for imputation with Beagle 5.2
# they are basically identical with their counterparts above

make_arguments_BEAGLE5.2 <- function(gt, out, write_args_to_txt = FALSE,
                                     arg_filename = "arguments_BEAGLE.txt",
                                     ref, map, seed = -99999,
                                     burnin = 3, iterations = 12,
                                     phase_states = 280, impute = "true",
                                     imp_states = 1600, imp_segment = 6.0,
                                     imp_step = 0.1, imp_nsteps = 7,
                                     cluster = 0.005, ap = "false",
                                     gp = "false", ne = 1000000,
                                     window = 40.0, overlap = 2.0, ...){
  arg <- c(as.list(environment()), list(...))
  arg <- arg[!names(arg) %in% c("write_args_to_txt", "arg_filename")]
  arg <- arg[as.logical(str_length(as.character(arg)))]
  names(arg) <- str_replace_all(names(arg), "\\_", "-")
  
  arg <- as.list(paste0(names(arg), "=", arg))
  
  if(write_args_to_txt){
    fwrite(as.list(format(Sys.time())), arg_filename, sep = "\n")
    fwrite(arg, arg_filename, sep = "\n", append = TRUE)
    # this writes the arguments_BEAGLE.txt file with the arguments in separate rows 
    # for better visualization
  }
  
  for_BEAGLE <- paste(arg, collapse = " ")
  return(for_BEAGLE)
}


# be aware that Beagle 5.2 returns an error when individuals that are in the
# reference are also included in the gt set
# they are detected based on their ID
impute_BEAGLE5.2 <- function(arguments, memory_limit_GB, verbose = TRUE) {
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start), "\n")}
  
  original_dir <- getwd()
  os <- Sys.info()["sysname"]
  if(!str_detect(os, pattern = "[[inux][indows]]")) { 
    #checks if the operating system is not windows or Linux
    stop("Operating system is neither Windows nor Linux. Function was written for
         Windows and Linux only.")
  }
  
  # check if java version is at least 1.8
  if(!check_version_java(at_least_vers = "1.8")) {
    stop("Java version is lower than 1.8. Please install at least java version
         1.8 to proceed.")
  }
  
  # checks if the BEAGLE program is present in the current directory or any of 
  # its subdirectories.If not, the loop goes to the parent directory and checks 
  # again if beagle is present in it and any of its subdirectories. If not,
  # it steps to the grandparent directory and searches again everywhere.
  # If BEAGLE still cannot be found, an error message is issued
  # If BEAGLE is found, the program is copied to the working directory of this 
  # script.
  # Note: if the directory tree doesn't change and many analyses are conducted it 
  # is recommended to put BEAGLE in the directory from which this script is executed
  # so that it does not need to search for beagle every time
  for(i in 1:3) {
    lf <- list.files(recursive = TRUE, full.names = TRUE)
    beagle_file <- lf[str_detect(lf, pattern = "/beagle.28Jun21.220.jar$")]
    if(length(beagle_file) == 0) {setwd("..")} 
    else {break}
  }
  rm(lf)
  if(length(beagle_file) == 0) {
    stop("The BEAGLE software should be present in the directory of this script or 
        in a subdirectory of two parent directories backwards.
        The name should be 'beagle.28Jun21.220.jar'.")
  }
  file.copy(beagle_file, original_dir) # copies beagle to the original directory
  setwd(original_dir)
  
  # for limiting the amount of memory used
  limit <- NULL
  if(!is.na(memory_limit_GB)) {
    limit <- paste0("-Xmx", memory_limit_GB, "g") 
  }
  
  begin <- Sys.time()
  system2(command = "java", 
          args = paste(limit, "-jar", "beagle.28Jun21.220.jar", arguments), 
          stdout = TRUE)
  time_needed <- difftime(Sys.time(), begin, units = "mins")
  
  file.remove("beagle.28Jun21.220.jar")
  
  if(verbose){cat('The process needed', format(Sys.time() - start), "\n")}
  return(time_needed)
}

# the only difference to impute_familywise_Beagle5.2() is that this function
# sets the IDs of the reference set to fake IDs so that they do't match
# ID in the gt set. In that case, Beagle 5.2 would complain
impute_familywise_BEAGLE5.2 <- function(arguments, marker_table, pedigree,
                                        memory_limit_GB = 2, reference_phase, mapo, 
                                        verbose = TRUE, restrict_ref_to_fam = FALSE) {######### update here the parameters
  start <- Sys.time()
  if(verbose) {cat("The process started at", format(start), "\n")}
  
  marker_table <- order_marker_table_by_pedigree(pedigree = pedigree, marker_table = marker_table)
  marker_table <- order_marker_table_by_mapo(mapo = mapo, marker_table = marker_table)
  ls_fam <- extract_family_marker_table(pedigree = pedigree, marker_table = marker_table)
  
  ##### create a subdirectory for this analysis
  dir.create("temporary")
  setwd("./temporary")
  ######################
  # the following is done to copy beagle in the current directory
  # for an explanation, look at 'impute_BEGLE'
  # this function would also work without these lines but it seems that it is slower 
  # then
  # I create two "temporary" files, one within the other. I do this so that I can 
  # save BEAGLE in one of them and run in the subdirectory. The problem is,
  # that 'impute_Beagle' will delete the beagle jar file if it is in its directory
  # after execution
  original_dir <- getwd()
  
  for(i in 1:3) {
    lf <- list.files(recursive = TRUE, full.names = TRUE)
    beagle_file <- lf[str_detect(lf, pattern = "/beagle.28Jun21.220.jar$")]
    if(length(beagle_file) == 0) {setwd("..")} 
    else {break}
  }
  rm(lf)
  if(length(beagle_file) == 0) {
    stop("The BEAGLE software should be present in the directory of this script or 
        in a subdirectory of two parent directories backwards.
        The name should be 'beagle.28Jun21.220.jar'.")
  }
  file.copy(beagle_file, original_dir) # copies beagle to the original directory
  
  setwd(original_dir)
  dir.create("temporary2")
  setwd("./temporary2")
  ########################################
  
  param <- stringr::str_split(arguments, pattern = " ")[[1]]
  name_map <- 
    stringr::str_remove(param[stringr::str_detect(param, pattern = "^map=")], "^map=")
  name_map <- stringr::str_remove(name_map, pattern = ".map$")
  map_PLINK <- convert_mapo_to_PLINK_map(mapo = mapo)
  write_PLINK_map(PLINK_map_object = map_PLINK, file_prefix = name_map) 
  
  use_ref <- any(stringr::str_detect(param, pattern = "^ref="))
  for(i in 1:length(ls_fam)) {
    if(verbose){
      start2 <- Sys.time()
      cat('Started imputing family', names(ls_fam)[[i]], 'of', length(ls_fam), 'families.\n')
      cat("The subprocess started at", format(start2), "\n")
    }
    mt_fam <- ls_fam[[i]]
    
    if(use_ref) {
      if(missing(reference_phase)) stop("You specified a reference in the arguments 
                                      but did not provide the reference phase object.")
      else {
        ref_phase_fam <- order_phase_by_mapo(mapo = mapo, phase = reference_phase)
        # here only the individuals of the family that are in the reference phase 
        # will be given to BEAGLE as a reference
        if(restrict_ref_to_fam){ # Updated 19.02. to allow to give references that are not in the family
          # this allows the code in the script for automatic imputation
          # to stay the same while still keeping the functionality
          # in case it it wanted
          ref_phase_fam <- 
            extract_individuals_phase(phase = ref_phase_fam
                                      , ID_inds = ref_phase_fam[ref_phase_fam[[1]] %in% colnames(mt_fam),1][[1]])
        }
        
        ref_name <- param[stringr::str_detect(param, pattern = "^ref=")]
        ref_name = stringr::str_remove(ref_name, "^ref=")
        as_gz <- stringr::str_detect(ref_name, pattern = '.gz$') 
        if(as_gz){ref_name <- stringr::str_remove(ref_name, pattern = ".gz$")}
        
        fakenames <- paste0('fakeID_',sort(rep(1:(nrow(ref_phase_fam)/2),2)))
        ref_phase_fam[[1]] <- fakenames
        write_phase_to_vcf(mapo = mapo, phase = ref_phase_fam,
                           outfile_name = ref_name, gzipped = as_gz, verbose = FALSE)
      }
    }
    
    gt_name <- param[stringr::str_detect(param, "^gt=")]
    gt_name <- stringr::str_remove(gt_name, pattern = "^gt=")
    as_gz <- stringr::str_detect(gt_name, pattern = '.gz$')
    if(as_gz){gt_name <- stringr::str_remove(gt_name, pattern = ".gz$")}
    
    write_marker_table_to_vcf(marker_table = mt_fam, mapo = mapo, 
                              outfile_name = gt_name, gzipped = as_gz, 
                              verbose = FALSE)
    if(i == 1){
      run_time <- impute_BEAGLE5.2(arguments = arguments
                                   , memory_limit_GB = memory_limit_GB
                                   , verbose = FALSE)
    }
    run_time <- run_time + impute_BEAGLE5.2(arguments = arguments
                                            , memory_limit_GB = memory_limit_GB
                                            , verbose = FALSE)
    
    out_name <- param[stringr::str_detect(param, "^out=")]
    out_name <- paste0(stringr::str_remove(out_name, pattern = "^out="), ".vcf.gz")
    
    imputed_phase <- read_vcf_to_phase(file = out_name, set_missing_gt_to = 9, verbose = FALSE)
    
    if(i == 1) {
      out_phase <- list(imputed_phase)
    } else{
      out_phase <- c(out_phase, list(imputed_phase))
    }
    if(verbose){cat('The subprocess needed', format(Sys.time() - start2), "\n")}
  }
  
  #### after reading in the results, the subdirectory and its files can be deleted
  setwd("..")
  unlink("./temporary2", recursive = TRUE)
  setwd("..")
  unlink("./temporary", recursive = TRUE)
  
  names(out_phase) <- names(ls_fam)
  if(verbose){cat('The process needed', format(Sys.time() - start), "\n")}
  return(list(out_phase = out_phase, run_time = run_time))
}