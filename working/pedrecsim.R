# Simulate recombination events in R


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(reshape)
library(plyr)
library(kinship2)

load("working/test.Rdata")
rm(female.r, male.r, snp.mafs, snp.names)

# Possible inputs

error.rate   <- 1e-5
missing.rate <- 0.001

minimum.r        <- 0.00001
minimum.r.male   <- 0.00001
minimum.r.female <- 0.000005

r        <- c(0.10, 0.10, 0.10, 0.2, 0.2, 0.2, 0.2, 0.2, 0)
r.male   <- c(0.20, 0.20, 0.20, 0.4, 0.4, 0.4, 0.4, 0.4, 0)
r.female <- c(0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0)

snp.names <- paste0("SNP", 1:10)
snp.mafs  <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2)

cM        <- c(10, 20, 30,  40,  60,  80, 100, 120, 140, 140)    # Give message saying that this is just an approximation
cM.male   <- c(20, 40, 60, 100, 140, 180, 220, 260, 300, 300)
cM.female <- c( 5, 10, 15,  25,  35,  45,  55,  65,  75,  75)

# Crossover intereference parameter

xover.min.r      <- NULL
xover.min.cM     <- NULL
xover.min.marker <- NULL

xover.min.r.male       <- NULL
xover.min.cM.male      <- NULL
xover.min.marker.male  <- NULL

xover.min.r.female      <- NULL
xover.min.cM.female     <- NULL
xover.min.marker.female <- NULL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Prepare pedigree information                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Check the pedigree format and rename variables for downstream analysis

pedigree.format <- function(ped, pedigree.type = "simple"){   # "plink"
  
  simple.ped.name.rules <- function(){
    writeLines("Pedigree columns must be named as follows:
    ID should be named ID or ANIMAL
    Mother should be MOTHER, MUM, MOM or DAM
    Father should be FATHER, DAD, POP or SIRE")
  }
  
  if(pedigree.type == "simple"){
    
    ped <- ped[,1:3]
    names(ped) <- toupper(names(ped))
    if(!names(ped)[1] %in% c("ID", "ANIMAL"))                 stop(simple.ped.name.rules())
    if(!names(ped)[2] %in% c("MUM", "MOM", "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE"))  stop(simple.ped.name.rules())
    if(!names(ped)[3] %in% c("MUM", "MOM", "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE")) stop(simple.ped.name.rules())
    
    names(ped)[which(names(ped) %in% c("MUM", "MOM", "MOTHER", "DAM"))] <- "MOTHER"
    names(ped)[which(names(ped) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"
    
    ped <- ped[,c("ANIMAL", "MOTHER", "FATHER")]

    for(i in 1:3) ped[which(is.na(ped[,i])),i] <- 0
    
  }
  
  if(pedigree.type == "plink"){
    
    if(ncol(ped) == 5){
      names(ped) <- c("FAMILY", "ANIMAL", "FATHER", "MOTHER", "SEX") #(1=male; 2=female; other=unknown)
      print(paste("Assuming columns ordered as:", paste(names(ped), collapse = " ")))
    }
    
    if(ncol(ped) == 6){
      names(ped) <- c("FAMILY", "ANIMAL", "FATHER", "MOTHER", "PHENOTYPE") #(1=male; 2=female; other=unknown)
      print(paste("Assuming columns ordered as:", paste(names(ped), collapse = " ")))
    }
    
    if(!ncol(ped) %in% c(5, 6)){
      stop("Number of columns does not match those expected of PLINK format")
    }
    
    # re-code missing values
    
    for(i in 2:4) ped[which(is.na(ped[,i])),i] <- 0
    for(i in 2:4) ped[which(ped[,i] == -9),i]  <- 0
    
    
  }
  
}

convert.cM.r <- function(cM){
  writeLines("Assuming 1cM is equivalent to r = 0.01")
  diff(cM)/100
}

# simulate.pedigree <- function(ped, r, r.female, r.male, cM

ped2 <- pedigree.format(ped, pedigree.type = "simple")

#~~ melt pedfile to get a unique row for each gamete transfer

transped <- melt(ped[,c("ANIMAL", "MOTHER", "FATHER")], id.vars = "ANIMAL")
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
  r.female[which(r.female == 0)] <- minimum.r
  r.male  [which(r.male   == 0)] <- minimum.r
}

#~~ Create a recombination template by sampling the probability of a crossover within an interval




template.list <- sapply(transped$Parent.ID.SEX, function(x){
  if(x == "MOTHER") remp.temp <- which(((runif(length(r.female)) < r.female) + 0L) == 1)
  if(x == "FATHER") remp.temp <- which(((runif(length(r.male  )) < r.male  ) + 0L) == 1)
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
      stop.pos <- c(rec.pos, length(r.female) + 1)
      
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
  vec <- rep(NA, (length(r.female) + 1)*2)
  
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
  
