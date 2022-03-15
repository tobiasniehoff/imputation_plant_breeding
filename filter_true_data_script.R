rm(list = ls())

# load packages
library(data.table)
library(microbenchmark)
library(stringr)
library(rlist)
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
setwd('./../data')
(odir <- getwd())

# load in functions
setwd("../R functions")
source("analysis.R")
source("general_functions.R")
source("simulation_functions.R")
source("BEAGLE_functions.R")
setwd(odir)

# read in data
marker_table <- fread("Sugar-beet-data1_Marker_012.txt")
pedigree <- fread("Sugar-beet-data1_ped.txt")
pedigree[is.na(pedigree)] <- 0

map <- read.table("Sugar-beet-data1_map.txt", header = TRUE)
map_LD1 <- read.table("Sugar-beet-data1_map-LD1.txt", header = TRUE)
map_LD2 <- read.table("Sugar-beet-data1_map-LD2.txt", header = TRUE)
map_LD3 <- read.table("Sugar-beet-data1_map-LD3.txt", header = TRUE)

#####################################################

mt <- marker_table
mapo <- map[map[[2]] != 99,]
mapo <- order_mapo_by_position(mapo)
ped <- extract_individuals_pedigree(pedigree, ID_names = colnames(mt)[-1])

mt <- extract_marker_marker_table(marker_table = mt, marker_IDs = mapo[[1]])
mt <- order_marker_table_by_mapo(mapo = map, marker_table = mt)
mt <- order_marker_table_by_pedigree(pedigree = ped, marker_table = mt)
mapo <- extract_marker_mapo(mapo = mapo, marker_IDs = mt[[1]])

######################################################
anal_mt <- analyze_inds_marker_table(mt)

##filter out all individuals with too much missing data
# find all individuals with more than 10% data missing
(ind_miss_filter <- anal_mt[[1]][which(anal_mt[[2]] > 0.1)])

# this shows in which families each parent is involved in
invol_parent <- extract_family_parents_pedigree(ped, per_parent = TRUE)

# this returns the family names that have at least one parent that had too much 
# missing data
miss_p_fam <- unique(unname(unlist(invol_parent[names(invol_parent) %in% ind_miss_filter])))

# removes all offspring that belong to a family with a parent of too much missing
# data

# removes all offspring that are in a family that has a parent with highly missing data
new_mt <- extract_family_marker_table(pedigree = ped, marker_table = mt, 
                            family_name = miss_p_fam, negate = TRUE)
# removes all individuals with highly missing data, i.e. the parents
new_mt <- extract_individuals_marker_table(marker_table = new_mt, ID_inds = ind_miss_filter
                                       ,negate = TRUE)

# takes the names of the marker table and only leaves these individuals in the 
# pedigree object, i.e., only genotypes available in the marker table should be 
# kept in the pedigree
new_ped <- extract_individuals_pedigree(pedigree = ped, ID_names = colnames(new_mt))

# check if pedigree is still resolvable, i.e., exclude loops, make sure all 
# are also entered on a row
check_resolvable_pedigree(new_ped)[[1]]


########################
# now all individuals with too much missing data and all their offspring have been 
# filtered out
a <- analyze_inds_marker_table(new_mt)
sum(a[[2]] > 0.1) # 0 individuals have 10% or more data missing

# mendelian error check
# (the variable names have no deeper meaning)
s <- analyze_mendelian_consistency(new_mt, new_ped)
e <- analyze_mendel_nerr_ind(new_mt, new_ped)
m <- analyze_mendel_nerr_marker(new_mt, new_ped)

# check how many individuals have more than 100 errors
sum(e[[2]] >= 100)
(x <- e[e[[2]] >= 100,])

data.table(extract_individuals_pedigree(pedigree = new_ped, ID_names = x[[1]]),Mend_err = x[[2]])
# one can see that the individuals with high numbers of incorrect data are coming 
# from only three families: Fam03, Fam11, Fam37

# this produces a list of marker tables. For each family one
ls_mt <- extract_family_marker_table(pedigree = new_ped, marker_table = new_mt)
names(ls_mt)

# check for family Fam03
(y <- analyze_mendel_nerr_ind(marker_table = ls_mt$Fam03, 
    pedigree = extract_individuals_pedigree(pedigree = new_ped,
                                            ID_names = colnames(ls_mt$Fam03))))
# look at the output and you can see that only G48 has such a high number of errors
# -> remove G48
new_ped <- extract_individuals_pedigree(new_ped, ID_names = y[y$nerrors > 100, 'Genotype'], negate = TRUE)


# check for family Fam11
(y <- analyze_mendel_nerr_ind(marker_table = ls_mt$Fam11, 
                              pedigree = extract_individuals_pedigree(pedigree = new_ped,
                                                                      ID_names = colnames(ls_mt$Fam11))))
# look at the output and you can see that only the individuals G471 to G478
# have such high numbers of errors. The others have a very low number.
new_ped <- extract_individuals_pedigree(new_ped, ID_names = y[y$nerrors > 100, 'Genotype'], negate = TRUE)


# check for family Fam37
(y <- analyze_mendel_nerr_ind(marker_table = ls_mt$Fam37, 
                              pedigree = extract_individuals_pedigree(pedigree = new_ped,
                                                                      ID_names = colnames(ls_mt$Fam37))))
# look at the output and you can see that only G1322 and G1329
# have very many mistakes. The others have none.
# -> remove G1322 and G1329
new_ped <- extract_individuals_pedigree(new_ped, ID_names = y[y$nerrors > 100, 'Genotype'], negate = TRUE)

## adjust marker table
new_mt <- extract_individuals_marker_table(marker_table = new_mt, 
                                           ID_inds = new_ped[[1]])

## check again

d <- analyze_mendel_nerr_ind(new_mt, new_ped)
# you can see that the number of errors went quite down
mean(e[[2]]) # data set including inds with many errors (number is average error per ind)
mean(d[[2]]) # data set excluding inds with many errors (number is average error per ind)

# now write the cleaned data set for later use
setwd('filtered_data/')
fwrite(new_mt, file = 'filtered_Sugar-beet-data1_Marker_012.txt', sep = '\t'
       ,quote = FALSE, col.names = TRUE)

fwrite(new_ped, file = 'filtered_Sugar-beet-data1_ped.txt', sep = '\t'
       ,quote = FALSE, col.names = TRUE)

fwrite(mapo, file = 'filtered_Sugar-beet-data1_map.txt', sep = '\t'
       ,quote = FALSE, col.names = TRUE)
# here I also save the low density maps again. Just to be sure
map_LD1 <- extract_marker_mapo(mapo = map_LD1, marker_IDs = mapo[[1]])
fwrite(map_LD1, file = 'filtered_Sugar-beet-data1_map1.txt', sep = '\t'
       ,quote = FALSE, col.names = TRUE)
map_LD2 <- extract_marker_mapo(mapo = map_LD2, marker_IDs = mapo[[1]])
fwrite(map_LD2, file = 'filtered_Sugar-beet-data1_map2.txt', sep = '\t'
       ,quote = FALSE, col.names = TRUE)
map_LD3 <- extract_marker_mapo(mapo = map_LD3, marker_IDs = mapo[[1]])
fwrite(map_LD3, file = 'filtered_Sugar-beet-data1_map3.txt', sep = '\t'
       ,quote = FALSE, col.names = TRUE)