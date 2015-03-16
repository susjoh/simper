# Simulate gene dropping events in R



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(reshape)
library(plyr)
library(kinship2)

pedigree <- read.table("data/test_pedigree.txt", header = T)
pedigree <- subset(pedigree, ID != 8160)

cohort   <- pedigree$BirthYear
minimum.cohort <- 1989
minimum.na.pause <- 3
minimum.cohort.gap <- 2
initial.sample.lag <- 4
nitt <- 10

#~~ Get founder allele frequency
x <- table(subset(pedigree, BirthYear == minimum.cohort)$Genotype)
founder.allele.freq <- (x[1] + 0.5*x[2])/sum(x)
attributes(founder.allele.freq) <- NULL



#~~ Organise the pedigree

ped <- pedigree.format.genedrop(pedigree)

#~~ Check that no parents were born in the same cohort as their offspring.Remove offspring if case

tempCohort <- rbind(data.frame(MinCohort = tapply(ped$Cohort, ped$MOTHER, min, na.rm = T)),
                    data.frame(MinCohort = tapply(ped$Cohort, ped$FATHER, min, na.rm = T)))
tempCohort$ANIMAL <- row.names(tempCohort)
tempCohort <- join(tempCohort, ped)

if(length(which(tempCohort$MinCohort == tempCohort$Cohort) > 0)){
  stop(paste(tempCohort[which(tempCohort$MinCohort == tempCohort$Cohort),"ANIMAL"],
             "has offspring born in the same cohort"))
}


tempCohort$MinCohort <- tempCohort$MinCohort - minimum.cohort.gap


tempCohort <- subset(tempCohort, ANIMAL %in% ped$ANIMAL[which(is.na(ped$Cohort))])

if(any(is.na(ped$Cohort))){
  for(i in 1:nrow(tempCohort)) ped$Cohort[which(ped$ANIMAL == tempCohort$ANIMAL[i])] <- tempCohort$MinCohort[i]
}


#~~ Subset the pedigree based on the minimum.cohort. NB. NA's need to be dealt with! NA individuals are
#   included if offspring are born minimum.na.pause years after minimum.cohort. They are placed in the
#   cohort two years previous to their first offspring


ped <- subset(ped, Cohort >= minimum.cohort)

if(any(!ped$MOTHER %in% ped$ANIMAL)) ped$MOTHER[which(!ped$MOTHER %in% ped$ANIMAL)] <- 0
if(any(!ped$FATHER %in% ped$ANIMAL)) ped$FATHER[which(!ped$FATHER %in% ped$ANIMAL)] <- 0

#~~ melt pedfile to get a unique row for each gamete transfer

transped <- melt(ped, id.vars = c("ANIMAL", "Cohort"))
transped$variable <- as.character(transped$variable)

#~~ Redefine columns

names(transped) <- c("Offspring.ID",  "Cohort", "Parent.ID.SEX", "Parent.ID")
transped$Offspring.ID <- as.character(transped$Offspring.ID)
transped$Parent.ID    <- as.character(transped$Parent.ID)


#~~ Create a list to store sampled genotypes

results.list <- list()

for(iteration in 1:nitt){
  
  print(paste("Running iteration", iteration, "of", nitt))
  
  
  haplo.list <- list()
  haplo.list[1:length(unique(transped$Offspring.ID))] <- list(list(MOTHER = NA, FATHER = NA))
  names(haplo.list) <- unique(transped$Offspring.ID)
  
  
  #~~ Sample genotypes over cohorts
  
  for(cohort in sort(unique(transped$Cohort))){
    
    if(nitt == 1) print(paste("Simulating haplotypes for cohort", cohort))
    
    #~~ sample founder cohort
    
    if(cohort == minimum.cohort){
      for(i in which(transped$Cohort == minimum.cohort)){
        
        if(transped$Parent.ID.SEX[i] == "MOTHER") haplo.list[transped$Offspring[i]][[1]]$MOTHER <- (runif(length(founder.allele.freq)) < founder.allele.freq) + 0L
        if(transped$Parent.ID.SEX[i] == "FATHER") haplo.list[transped$Offspring[i]][[1]]$FATHER <- (runif(length(founder.allele.freq)) < founder.allele.freq) + 0L
        
      }
    } else {
      
      #~~ sample IDs with parents from cohort
      
      samp.vec <- 0
      
      for(j in which(transped$Cohort == cohort & transped$Parent.ID != 0)){
        samp <- haplo.list[transped$Parent.ID[j]][[1]][sample.int(2, 1)][[1]]
        haplo.list[transped$Offspring.ID[j]][[1]][transped$Parent.ID.SEX[j]] <- samp
        samp.vec <- samp.vec + samp
      }    
      
      if(cohort - minimum.cohort < initial.sample.lag) sample.allele.freq <- founder.allele.freq
      if(cohort - minimum.cohort >= initial.sample.lag) sample.allele.freq <- 1 - samp.vec/length(which(transped$Cohort == cohort & transped$Parent.ID != 0))
      
      for(k in which(transped$Cohort == cohort & transped$Parent.ID == 0)){
        
        if(transped$Parent.ID.SEX[k] == "MOTHER") haplo.list[transped$Offspring[k]][[1]]$MOTHER <- (runif(length(sample.allele.freq)) < sample.allele.freq) + 0L
        if(transped$Parent.ID.SEX[k] == "FATHER") haplo.list[transped$Offspring[k]][[1]]$FATHER <- (runif(length(sample.allele.freq)) < sample.allele.freq) + 0L
        
      }     
    }
  }
  
  
  #~~ Condense haplotypes into genotypes
  
  genotype.list <- list()
  
  for(i in 1:length(haplo.list)){
    
    vec <- rep(NA, 2)
    
    vec[seq(1, length(vec), 2)] <- haplo.list[[i]][1][[1]]
    vec[seq(2, length(vec), 2)] <- haplo.list[[i]][2][[1]]
    
    genotype.list[[i]] <- vec
  }
  
  names(genotype.list) <- names(haplo.list)
  
  
  x <- t(data.frame(genotype.list))
  row.names(x) <- names(genotype.list)
  
  results.list[[iteration]] <- x
  
  
}

beepr::beep()

# all(row.names(x) == ped$ANIMAL)
# 
# ped$Genotype <- rowSums(x)
# 
# head(ped)
# 
# newped <- pedigree[,c(1, 3, 2)]














system.time({
  
  test2 <- simgenotypes(list(loc1 = c(founder.allele.freq, 1 - founder.allele.freq)), ped = newped)
  test3 <- data.frame(cbind(MBgeno = test2$G$loc1, ANIMAL = test2$Gid))
  test3$MBgeno <- test3$MBgeno - 1
  
})
ped <- join(ped, test3)
1 - (tapply(ped$Genotype, ped$Cohort, sum)/2)/tapply(ped$Genotype, ped$Cohort, length)
1 - (tapply(ped$MBgeno, ped$Cohort, sum)/2)/tapply(ped$MBgeno, ped$Cohort, length)

