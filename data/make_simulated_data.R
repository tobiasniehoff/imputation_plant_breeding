# the path should point to the directory in which this script is located
original_directory <- getwd()

# sourcing functions
setwd("../R functions")

source("analysis.R")
source("general_functions.R")
source("BEAGLE_functions.R") # really needed?
source("simulation_functions.R")

setwd(original_directory)

# this script will simulate some population history with MoBPS
# starting with pop size 1000 and linearly decreasing to 100 after 100 generations
# followed by 4 generations of linear population expansion from 100 to 1000 individuals
# followed by 10 generations of phenotypic selection based on F3 plants
# with heritability of 0.5 and selecting 30 best plants out of 2000; random mating
# then, 30 plants are randomly selected and selfed to times to get F3
# these are used as parents in the dummy data set
# the 20k marker for the high-density chip are randomly sampled from markers that
# have a minor allele frequency of higher than 0.01 in the last 2000 plant population
# the low-density panel are nested and markers are selected randomly

# In Windows use
# RandomFieldsUtils version 0.6.6
# miraculix version 1.0.0.1

# In Linux use
# RandomFieldsUtils version 1.0.6
# miraculix version 1.0.5
library(MoBPS)
library(stringr)


set.seed(1)
n.snps <- 63000
n.qtls <- 2700
n.chr <- 9
ninds <- 1000
n.snps <- n.snps+n.qtls # total number of positions that need to be simulated

# effect sizes are sampled from a gamma distribution with shape parameter 0.4
effect0 <- effect2 <- qtl.effects <- rgamma(n = n.qtls, shape = 0.4)

# sample which SNPs should be positions of true QTLs
qtl.snps <- sort(sample(1:n.snps, size = n.qtls))

SNP <- qtl.snps
SNP <- (SNP-1)%%(n.snps/n.chr)+1

pos <- sample(x = n.qtls, size = n.qtls/2,replace = F)
# effects assigned for the 0 allele
effect0[pos] <- qtl.effects[pos]
# effects assigned for the 1 allele
effect2[pos] <- 0
effect0[-pos] <- 0
effect2[-pos] <- qtl.effects[-pos]
effect1 <- rowMeans(cbind(effect0, effect2))
chromosome <- sort(rep(1:n.chr,n.qtls/n.chr))
effect.matrix <- cbind(SNP, chromosome, effect0, effect1, effect2)

# creating the founder population
# with 9 chromosomes
# and each chromosome is 80 centi Morgan long
population <- creating.diploid(nsnp = n.snps, 
                               nindi = ninds, "random",
                               chr.nr = n.chr, chromosome.length = 0.8,
                               snps.equidistant = T, 
                               name.cohort = "Founder", sex.quota = 0)

# adding the first trait 
population <- creating.trait(population, real.bv.add = effect.matrix)


# simulating 100 generations of random mating with decreasing 
# population size over time
n.gen <- 100
d <- data.frame(x=c(1,n.gen),y=c(1000,100))
l <- lm(y~x, d)

for(i in 1:n.gen){
  ninds <- round(predict(l, rbind(d,c(i,NA)))[[3]],0)
  population <- breeding.diploid(population, breeding.size = c(ninds,0))
  population <- new.base.generation(population, 
                                    base.gen = length(population$breeding), 
                                    delete.previous.gen = T, 
                                    delete.bve.data = T, delete.breeding.totals = T)
}
pop.store.100 <- population

# simulating sudden population expansion in 4 generations from 100 to 1000
# population size over time
n.gen <- 4
d <- data.frame(x=c(1,n.gen),y=c(100,1000))
l <- lm(y~x, d)

for(i in 1:n.gen){
  ninds <- round(predict(l, rbind(d,c(i,NA)))[[3]],0)
  population <- breeding.diploid(population, breeding.size = c(ninds,0))
  population <- new.base.generation(population, 
                                    base.gen = length(population$breeding), 
                                    delete.previous.gen = T, 
                                    delete.bve.data = T, delete.breeding.totals = T)
}
pop.store.100.4 <- population

population <- pop.store.100.4
# simulating 10 generations of phenotypic selection 
# based on 1 trait
# with heritability 0.5
# with the best 30 F3 lines selected out of 2000
# with random mating
# with distinct generations

n.gen <- 10
for(i in 1:n.gen){
  
  # performing phenotyping
  population <- breeding.diploid(population, 
                                 heritability = 0.5, 
                                 phenotyping.gen = length(population$breeding),
                                 verbose = F)
  
  population <- breeding.diploid(population,
                                 selection.size = c(30,0),
                                 breeding.size = c(2000,0),
                                 selection.criteria = "pheno")
  
  # creating F3 plants from F1
  for(index in 1:2){ 
    population <- breeding.diploid(population,
                                   selfing.mating = TRUE)
  }
  
  population <- new.base.generation(population, 
                                    base.gen = length(population$breeding), 
                                    delete.previous.gen = T, 
                                    delete.bve.data = T, delete.breeding.totals = T)
}
pop.store.100.4.10 <- population

geno <- get.geno(population, gen = length(population$breeding))
# calculates allele frequency per marker
AF <- rowMeans(geno)/2

SNP.IDs <- rownames(geno)

# removing all SNPs with a minor allele frequency below 0.01
SNP.IDs <- SNP.IDs[AF > 0.01 & AF < (1-0.01)]

# Now creating the genetic maps
map <- get.map(population, use.snp.nr = F)
map <- data.frame(Marker = map[,2], 
                  Chrom = as.integer(map[,1]), 
                  Pos_start = as.numeric(map[,4]))
# position is in bp. In the end, the unit should be cM
# 100,000,000 bp is equal to 1 Morgan

# Calculating in Morgan
map[,3] <- map[,3]/100000000
# Calculating in centi Morgan
map[,3] <- map[,3]*100

# selecting SNPs for the high-density 20k panel
SNP.20k <- SNP.IDs[sort(sample(1:length(SNP.IDs), size = 20000))]

map.20k <- map[map[,1] %in% SNP.20k,]

# Now creation of the 30 parents
# random selection of 30 plants and F1 creation
population <- breeding.diploid(population,
                                  selection.size = c(30,0),
                                  breeding.size = c(30,0))

# generating of F3 (S2) inbreds by 2 selfings
selfing_population <- population 
for(index in 1:2){ 
  selfing_population <- breeding.diploid(selfing_population, 
                                         breeding.size = c(30,0), 
                                         selfing.mating = TRUE) 
}

# writing the data to vcf format
get.vcf(selfing_population, path = "founder_population", gen = length(selfing_population$breeding))

# reading in the data and converting to phase (custom) format
parents_phase <- read_vcf_to_phase("founder_population.vcf", verbose = T)
file.remove("founder_population.vcf")

# restrict the phase object to only those markers that are on the 20k chip
parents_phase <- extract_marker_phase(parents_phase, marker_IDs = map.20k[,1])

# renaming the SNPs
for(i in as.integer(unique(str_extract(SNP.20k, pattern = "[:digit:]")))){
  positions <- str_detect(SNP.20k, pattern = paste0("Chr", i))
  SNP.20k[positions] <- paste0("Chr",i,"SNP",1:sum(positions))
}

# selecting SNPs for the low-density 3k panel
SNP.3k <- SNP.20k[sort(sample(1:length(SNP.20k), size = 3000))]

# selecting SNPs for the low-density 2k panel
SNP.2k <- SNP.3k[sort(sample(1:length(SNP.3k), size = 2000))]

# selecting SNPs for the low-density 1k panel
SNP.1k <- SNP.2k[sort(sample(1:length(SNP.2k), size = 1000))]

map.20k[,1] <- SNP.20k
map.3k <- map.20k[map.20k[,1] %in% SNP.3k,]
map.2k <- map.20k[map.20k[,1] %in% SNP.2k,]
map.1k <- map.20k[map.20k[,1] %in% SNP.1k,]

write.table(map.20k, file = "simulated_map.txt", quote = F, row.names = F)
write.table(map.3k, file = "simulated_map3.txt", quote = F, row.names = F)
write.table(map.2k, file = "simulated_map2.txt", quote = F, row.names = F)
write.table(map.1k, file = "simulated_map1.txt", quote = F, row.names = F)

colnames(parents_phase)[-1] <- map.20k[,1]
parents_marker_table <- convert_phase_to_marker_table(parents_phase)
anal_parents <- analyze_marker_marker_table(parents_marker_table)
# number of markers that are fixed in the parents
sum(anal_parents$MAF == 0)
# if equates to 0, then every parent is genetically unique
sum(duplicated(t(parents_marker_table[,-1])))

# change parent identifiers 
ped <- read.table("simulated_data_ped.txt", header = T)
parents_phase[,1] <- rep(ped[1:30,1],each = 2)

# this function is defined in "simulation_functions.R"
# briefly, this function detects all unique pairings in the provided pedigree
# then it executes every unique cross for as many times as there should be
# offspring in the family (n_off.family)
# the offspring are creating according to the desired type (F3 means 2 selfings
# based on F1)
sim_data <- make_simulated_data(offspring_type = "F3", 
                                parent_type = "real",
                                n_off.family = 50, 
                                phase = parents_phase, 
                                pedigree = ped, 
                                map = map.20k, 
                                seed = 100)

sim_phase <- sim_data$sim_phase
sim_ped <- sim_data$sim_ped

write_phase_to_vcf(phase = sim_phase, mapo = map.20k, 
                   outfile_name = "simulated_data.vcf", gzipped = T)
write.table(sim_ped, file = "simulated_data_ped.txt", quote = F, row.names = F)

# done creating the dummy data set
################################################################################
# investigating the genetic properties of the created population

marker_table <- convert_phase_to_marker_table(sim_phase)

# get overview of MAF distribution
analysis.marker <- analyze_marker_marker_table(marker_table)
hist(analysis.marker$MAF, breaks = 15)

# get overview of SNP homozygosity distribution
analysis.individuals <- analyze_inds_marker_table(marker_table)
hist(analysis.individuals$prop_homozygous)