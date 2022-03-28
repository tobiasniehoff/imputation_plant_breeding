source('general_functions.R')
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

##############################
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


##############################
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

################################
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

# same as make_arguments_BEAGLE5.2() but with the addition of the em parameter
# at its default value. Using make_arguments_BEAGLE5.2() and specifying the 
# em parameter works as well
make_arguments_BEAGLE5.3 <- function(gt, out, write_args_to_txt = FALSE,
                                     arg_filename = "arguments_BEAGLE.txt",
                                     ref, map, seed = -99999,
                                     burnin = 3, iterations = 12,
                                     phase_states = 280, impute = "true",
                                     imp_states = 1600, imp_segment = 6.0,
                                     imp_step = 0.1, imp_nsteps = 7,
                                     cluster = 0.005, ap = "false",
                                     gp = "false", ne = 1000000,
                                     window = 40.0, overlap = 2.0, 
                                     em='true',...){
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

impute_BEAGLE5.3 <- function(arguments, memory_limit_GB, verbose = TRUE) {
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
    beagle_file <- lf[str_detect(lf, pattern = "/beagle.08Feb22.fa4.jar$")]
    if(length(beagle_file) == 0) {setwd("..")} 
    else {break}
  }
  rm(lf)
  if(length(beagle_file) == 0) {
    stop("The BEAGLE software should be present in the directory of this script or 
        in a subdirectory of two parent directories backwards.
        The name should be 'beagle.08Feb22.fa4.jar'.")
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
          args = paste(limit, "-jar", "beagle.08Feb22.fa4.jar", arguments), 
          stdout = TRUE)
  time_needed <- difftime(Sys.time(), begin, units = "mins")
  
  file.remove("beagle.08Feb22.fa4.jar")
  
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

# same as impute_familywise_BEAGLE5.2(). Only difference is that this function is looking for the file 
# beagle.08Feb22.fa4.jar which is BEAGLE 5.3
impute_familywise_BEAGLE5.3 <- function(arguments, marker_table, pedigree,
                                        memory_limit_GB = 2, reference_phase, mapo, 
                                        verbose = TRUE, 
                                        restrict_ref_to_fam = FALSE) {######### update here the parameters
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
    beagle_file <- lf[str_detect(lf, pattern = "/beagle.08Feb22.fa4.jar$")]
    if(length(beagle_file) == 0) {setwd("..")} 
    else {break}
  }
  rm(lf)
  if(length(beagle_file) == 0) {
    stop("The BEAGLE software should be present in the directory of this script or 
        in a subdirectory of two parent directories backwards.
        The name should be 'beagle.08Feb22.fa4.jar'.")
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
      run_time <- impute_BEAGLE5.3(arguments = arguments
                                   , memory_limit_GB = memory_limit_GB
                                   , verbose = FALSE)
    }
    run_time <- run_time + impute_BEAGLE5.3(arguments = arguments
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