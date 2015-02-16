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

error.rate   <- 1e-4
missing.rate <- 0.001

cM        <- c(0,  5, 10,  15,   25,   35,   45,   55,   65,   65)
cM.male   <- c(0, 10, 20,  30,   50,   70,   90,  110,  130,  130)
cM.female <- c(0, 2.5,  5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 32.5)

snp.names     <- paste0("SNP", 1:10)
founder.mafs  <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2)

# Crossover intereference parameter

xover.min.cM            <- 25
xover.min.cM.male       <- NULL
xover.min.cM.female     <- NULL



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Prepare pedigree information                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#simperSNP <- function(ped, 


pedigree.type <- "simple"
phased.output <-  TRUE
out.format    <- "PLINKlike" # ?

cM.female <- NULL
cM.male <- NULL


#~~ Check and format the recombination units

if( is.null(cM) &  is.null(cM.male) &  is.null(cM.female)) stop   ("No recombination rates specified.")
if(!is.null(cM) &  is.null(cM.male) & !is.null(cM.female)) warning("Recombination rate only defined in one sex. Missing sex defaults to r.")
if(!is.null(cM) & !is.null(cM.male) &  is.null(cM.female)) warning("Recombination rate only defined in one sex. Missing sex defaults to r.")
if( is.null(cM) &  is.null(cM.male) & !is.null(cM.female)) stop   ("Recombination rate only defined in one sex. Use r= ... ")
if( is.null(cM) & !is.null(cM.male) &  is.null(cM.female)) stop   ("Recombination rate only defined in one sex. Use r= ...")

if(!is.null(cM) &  is.null(cM.male) &  is.null(cM.female)) message("Assuming that recombination rates are equal in both sexes.")

#~~ Use cM if no sex specific information given

if(is.null(cM.male  )) cM.male   <- cM
if(is.null(cM.female)) cM.female <- cM


#~~ Convert to R

message("Assuming 1cM is equivalent to r = 0.01")
r        <- diff(cM/100)
r.male   <- diff(cM.male/100)
r.female <- diff(cM.female/100)

#~~ Format the pedigree

ped2 <- pedigree.format(ped, pedigree.type = pedigree.type)

#~~ melt pedfile to get a unique row for each gamete transfer

transped <- melt(ped[,c("ANIMAL", "MOTHER", "FATHER")], id.vars = "ANIMAL")
transped$variable <- as.character(transped$variable)

#~~ assign pedfile IDs to cohort and merge with transped

cohorts <- data.frame(ANIMAL = ped[,1],
                      Cohort = kindepth(ped[,1], ped[,3], ped[,2]))

transped <- join(transped, cohorts)

#~~ Redefine columns

names(transped) <- c("Offspring.ID", "Parent.ID.SEX", "Parent.ID", "Cohort")
transped$Key    <- paste(transped$Parent.ID, transped$Offspring.ID, sep = "_")

#~~ Recode founder gametes to 0 e.g. if one parent is unknown

transped$Cohort[which(transped$Parent.ID == 0)] <- 0


#~~ Create objects in which to put recombination count, errors and missing genotypes

error.tab   <- NULL
missing.tab <- NULL

#~~ Tidy up crossover interference parameters

if( is.null(xover.min.cM) &  is.null(xover.min.cM.male) &  is.null(xover.min.cM.female)) message ("Assuming no crossover interference.")
if(!is.null(xover.min.cM) &  is.null(xover.min.cM.male) & !is.null(xover.min.cM.female)) warning("Crossover interference parameter only defined in one sex. Missing sex defaults to xover.min.cM.")
if(!is.null(xover.min.cM) & !is.null(xover.min.cM.male) &  is.null(xover.min.cM.female)) warning("Crossover interference parameter only defined in one sex. Missing sex defaults to xover.min.cM.")
if( is.null(xover.min.cM) &  is.null(xover.min.cM.male) & !is.null(xover.min.cM.female)) stop   ("Crossover interference parameter only defined in one sex. Use xover.min.cM= ... ")
if( is.null(xover.min.cM) & !is.null(xover.min.cM.male) &  is.null(xover.min.cM.female)) stop   ("Crossover interference parameter only defined in one sex. Use xover.min.cM= ...")

if(!is.null(xover.min.cM) &  is.null(xover.min.cM.male) &  is.null(xover.min.cM.female)) message("Assuming that crossover interference is equal in both sexes.")

#~~ Use r if no sex specific information given

if(is.null(xover.min.cM.male  )) xover.min.cM.male   <- xover.min.cM
if(is.null(xover.min.cM.female)) xover.min.cM.female <- xover.min.cM

#~~ Convert to r

xover.min.r        <- xover.min.cM/100
xover.min.r.male   <- xover.min.cM.male/100
xover.min.r.female <- xover.min.cM.female/100


#~~ Create a recombination template by sampling the probability of a crossover within an interval

r.cumu.female <- cumsum(r.female)
r.cumu.male   <- cumsum(r.male)

suppressWarnings({
  
  template.list <- sapply(transped$Parent.ID.SEX, function(x){
    
    if(x == "MOTHER"){
      
      remp.temp <- which(((runif(length(r.female)) < r.female) + 0L) == 1)
      
      if(length(remp.temp) > 1 & !is.null(xover.min.r.female)){
        while(min(diff(r.cumu.female[remp.temp])) < xover.min.r.female) {
          remp.temp <- which(((runif(length(r.female)) < r.female) + 0L) == 1)
        }
      }
    }
    
    
    if(x == "FATHER"){
      remp.temp <- which(((runif(length(r.male  )) < r.male  ) + 0L) == 1)
      
      if(length(remp.temp) > 1 & !is.null(xover.min.r.male)){
        while(min(diff(r.cumu.female[remp.temp])) < xover.min.r.male) {
          remp.temp <- which(((runif(length(r.male)) < r.male) + 0L) == 1)
        }
      }
    }
    
    remp.temp
  })

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
  if(transped$Parent.ID.SEX[i] == "MOTHER") haplo.list[transped$Offspring[i]][[1]]$MOTHER <- (runif(length(founder.mafs)) < founder.mafs) + 0L
  if(transped$Parent.ID.SEX[i] == "FATHER") haplo.list[transped$Offspring[i]][[1]]$FATHER <- (runif(length(founder.mafs)) < founder.mafs) + 0L
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

genotype.list <- list()

if(phased.output == TRUE) message("Genotypes are phased - Maternal, Paternal")

for(i in 1:length(haplo.list)){
  
  vec <- rep(NA, (length(r.female) + 1)*2)
  
  if(phased.output == TRUE){
    
    vec[seq(1, length(vec), 2)] <- haplo.list[[i]][1][[1]]
    vec[seq(2, length(vec), 2)] <- haplo.list[[i]][2][[1]]
    
  } else {
    
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


#~~ Create errors in the dataset

geno.count <- nrow(x) * (ncol(x)/2)

if(as.integer(geno.count * error.rate) > 0){

error.pos <- sample(1:geno.count, size = as.integer(geno.count * error.rate))

error.row <- error.pos/ncol(x)


