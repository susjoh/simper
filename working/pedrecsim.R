# Simulate recombination events in R


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(reshape)
library(plyr)
library(kinship2)


load("test.Rdata")

error.rate   <- 1e-5
missing.rate <- 0.001
minimum.r    <- 0.00001

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Prepare pedigree information                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Check the pedigree format and rename variables for downstream analysis

ped.name.rules <- writeLines("Pedigree columns must be named as follows:
  Column 1 should be named ID or ANIMAL
  Column 2 should be MOTHER, MUM, MOM or DAM
  Column 3 should be FATHER, DAD, POP or SIRE")

pedigree.check <- function(ped){
  
  ped <- ped[,1:3]
  names(ped) <- toupper(names(ped))
  if(!names(ped)[1] %in% c("ID", "ANIMAL"))                 stop(ped.name.rules())
  if(!names(ped)[2] %in% c("MUM", "MOM", "MOTHER", "DAM"))  stop(ped.name.rules())
  if(!names(ped)[3] %in% c("DAD", "POP", "FATHER", "SIRE")) stop(ped.name.rules())
  
}

pedigree.format <- function(ped){
  
  pedigree.check(ped)
  
  names(ped) <- c("ANIMAL", "MOTHER", "FATHER")
  
  for(i in 1:3) ped[which(is.na(ped[,i])),i] <- 0
  
}

pedigree.format(ped)

#~~ melt pedfile to get a unique row for each gamete transfer

transped <- melt(ped, id.vars = "ANIMAL")
transped$variable <- as.character(transped$variable)

#~~ assign pedfile IDs to cohort and merge with transped

cohorts <- data.frame(ANIMAL = ped[,1],
                      Cohort = kindepth(ped[,1], ped[,3], ped[,2]))

transped <- join(transped, cohorts)

#~~ Redefine columns

names(transped) <- c("Offspring.ID", "Parent.ID.SEX", "Parent.ID", "Cohort")
transped$Key <- paste(transped$Parent.ID, transped$Offspring.ID, sep = "_")


#~~ Recode founder gametes to 0 e.g. if one parent is unknown

transped$Cohort[which(transped$Parent.ID == 0)] <- 0



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Run simulation                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Create objects in which to put recombination count, errors and missing genotypes

error.tab <- NULL
missing.tab <- NULL

#~~ recode r = 0 if minimum.r is specified

if(!is.null(minimum.r)){
  female.r[which(female.r == 0)] <- minimum.r
  male.r  [which(male.r   == 0)] <- minimum.r
}

#~~ Create a recombination template by sampling the probability of a crossover within an interval
#   NB. THIS DOES NOT MODEL CROSSOVER INTERFERENCE


template.list <- sapply(transped$Parent.ID.SEX, function(x){
  if(x == "MOTHER") remp.temp <- which(((runif(length(female.r)) < female.r) + 0L) == 1)
  if(x == "FATHER") remp.temp <- which(((runif(length(male.r  )) < male.r  ) + 0L) == 1)
  remp.temp
})

#~~ Add key to list

names(template.list) <- transped$Key

#~~ Calculate recombination count

transped$RecombCount <- unlist(lapply(template.list, length))
transped$RecombCount[which(transped$Cohort == 0)] <- NA


#~~ Create founder haplotypes by sampling allele frequencies

head(transped)

#~~ Create a list to sample haplotypes

haplo.list <- list()
haplo.list[1:length(unique(transped$Offspring.ID))] <- list(list(MOTHER = NA, FATHER = NA))
names(haplo.list) <- unique(transped$Offspring.ID)

#~~ Sample the founder haplotypes

system.time(for(i in which(transped$Cohort == 0)){
  if(transped$Parent.ID.SEX[i] == "MOTHER") haplo.list[transped$Offspring[i]][[1]]$MOTHER <- (runif(length(snp.mafs)) < snp.mafs) + 0L
  if(transped$Parent.ID.SEX[i] == "FATHER") haplo.list[transped$Offspring[i]][[1]]$FATHER <- (runif(length(snp.mafs)) < snp.mafs) + 0L
})


#~~ Sample the non-founder haplotypes by cohort. This loops through cohorts sequentially
#   parental haplotypes must exist before sampling.  

for(cohort in 1:max(transped$Cohort)){
  
  print(paste("Simulating haplotypes for cohort", cohort, "of", max(transped$Cohort)))
  
  for(j in which(transped$Cohort == cohort)){
    
    #~~ Pull out recombination information from template.list
    
    rec.pos <- template.list[transped$Key[j]][[1]]
    
    #~~ If there are no recombination events, sample one of the parental haplotypes at random
    
    if(length(rec.pos) == 0){
      haplo.list[transped$Offspring.ID[j]][[1]][transped$Parent.ID.SEX[j]] <- haplo.list[transped$Parent.ID[j]][[1]][sample.int(2, 1)]
    }
    
    #~~ If there are recombination events, sample the order of the parental haplotypes,
    #   exchange haplotypes and then transmit haplotype to offspring
    
    if(length(rec.pos) > 0){
      
      parental.haplotypes <- haplo.list[transped$Parent.ID[j]][[1]][sample.int(2, 2, replace = F)]
      
      start.pos <- c(1, rec.pos[1:(length(rec.pos))] + 1)
      stop.pos <- c(rec.pos, length(female.r) + 1)
      
      fragments <- list()
      
      for(k in 1:length(start.pos)){
        if(k %% 2 != 0) fragments[[k]] <- parental.haplotypes[1][[1]][start.pos[k]:stop.pos[k]]
        if(k %% 2 == 0) fragments[[k]] <- parental.haplotypes[2][[1]][start.pos[k]:stop.pos[k]]  
      }
      
      haplo.list[transped$Offspring.ID[j]][[1]][transped$Parent.ID.SEX[j]] <- list(unlist(fragments))
      
    }
  }
}

#~~ Condense haplotypes into genotypes

str(haplo.list)
geno.format <- "phased" # "ordered"
out.format <- "PLINKlike" # ?

genotype.list <- list()

for(i in 1:length(haplo.list)){
  vec <- rep(NA, (length(female.r) + 1)*2)
  
  if(geno.format == "phased"){
    vec[seq(1, length(vec), 2)] <- haplo.list[[i]][1][[1]]
    vec[seq(2, length(vec), 2)] <- haplo.list[[i]][2][[1]]
  }
  
  if(geno.format == "ordered"){
    haplo.tab <- data.frame(A = haplo.list[[i]][1][[1]],
                            B = haplo.list[[i]][2][[1]])
    
    vec <- as.vector(apply(haplo.tab, 1, sort))
  }
   
  genotype.list[[i]] <- vec
}

names(genotype.list) <- names(haplo.list)


x <- t(data.frame(genotype.list))
row.names(x) <- names(genotype.list)


if(!is.null(snp.names)) colnames(x) <- rep(snp.names, each = 2)
  
