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

# minimum.r        <- 0.00001
# minimum.r.male   <- 0.00001
# minimum.r.female <- 0.000005

r        <- c(0.10, 0.10, 0.10, 0.2, 0.2, 0.2, 0.2, 0.2, 0)
r.male   <- c(0.20, 0.20, 0.20, 0.4, 0.4, 0.4, 0.4, 0.4, 0)
r.female <- c(0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0)

snp.names <- paste0("SNP", 1:10)
snp.mafs  <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2)

# Crossover intereference parameter

xover.min.r            <- 0.5
xover.min.r.male       <- NULL
xover.min.r.female     <- NULL



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Prepare pedigree information                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

r             = NULL
r.male        = NULL
r.female      = NULL
pedigree.type = "simple"
minimum.r = NULL
phased.output <- TRUE
out.format <- "PLINKlike" # ?

#~~ Check and format the recombination units

if(units == "cM"){
  message("Assuming 1cM is equivalent to r = 0.01")
  if(!is.null(r       )) r        <- r/100
  if(!is.null(r.male  )) r.male   <- r.male/100
  if(!is.null(r.female)) r.female <- r.female/100
}


if(units = "r" & max(na.omit(r)) > 0.5)){
  stop("Maximum r is greater than 0.5. Check that values are not in cM difference, or use units = \"cM\"")
}

if( is.null(r) &  is.null(r.male) &  is.null(r.female)) stop   ("No recombination rates specified.")
if(!is.null(r) &  is.null(r.male) & !is.null(r.female)) warning("Recombination rate only defined in one sex. Missing sex defaults to r.")
if(!is.null(r) & !is.null(r.male) &  is.null(r.female)) warning("Recombination rate only defined in one sex. Missing sex defaults to r.")
if( is.null(r) &  is.null(r.male) & !is.null(r.female)) stop   ("Recombination rate only defined in one sex. Use r= ... ")
if( is.null(r) & !is.null(r.male) &  is.null(r.female)) stop   ("Recombination rate only defined in one sex. Use r= ...")

if(!is.null(r) &  is.null(r.male) &  is.null(r.female)) message("Assuming that recombination rates are equal in both sexes.")

#~~ Use r if no sex specific information given

if(is.null(r.male  )) r.male   <- r
if(is.null(r.female)) r.female <- r

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
  r.female[which(r.female < minimum.r)] <- minimum.r
  r.male  [which(r.male  < minimum.r)] <- minimum.r
}

#~~ Tidy up crossover interference parameters

if( is.null(xover.min.r) &  is.null(xover.min.r.male) &  is.null(xover.min.r.female)) message ("Assuming no crossover interference.")
if(!is.null(xover.min.r) &  is.null(xover.min.r.male) & !is.null(xover.min.r.female)) warning("Crossover interference parameter only defined in one sex. Missing sex defaults to xover.min.r.")
if(!is.null(xover.min.r) & !is.null(xover.min.r.male) &  is.null(xover.min.r.female)) warning("Crossover interference parameter only defined in one sex. Missing sex defaults to xover.min.r.")
if( is.null(xover.min.r) &  is.null(xover.min.r.male) & !is.null(xover.min.r.female)) stop   ("Crossover interference parameter only defined in one sex. Use xover.min.r= ... ")
if( is.null(xover.min.r) & !is.null(xover.min.r.male) &  is.null(xover.min.r.female)) stop   ("Crossover interference parameter only defined in one sex. Use xover.min.r= ...")

if(!is.null(xover.min.r) &  is.null(xover.min.r.male) &  is.null(xover.min.r.female)) message("Assuming that crossover interference is equal in both sexes.")

#~~ Use r if no sex specific information given

if(is.null(xover.min.r.male  )) xover.min.r.male   <- xover.min.r
if(is.null(xover.min.r.female)) xover.min.r.female <- xover.min.r


#~~ Create a recombination template by sampling the probability of a crossover within an interval

r.cumu.female <- cumsum(r.female)
r.cumu.male   <- cumsum(r.male)




template.list <- sapply(transped$Parent.ID.SEX, function(x){
  
  if(is.null(xover.min.r.female) & x == "MOTHER") remp.temp <- which(((runif(length(r.female)) < r.female) + 0L) == 1)
  if(is.null(xover.min.r.male)   & x == "FATHER") remp.temp <- which(((runif(length(r.male  )) < r.male  ) + 0L) == 1)
    
  if(!is.null(xover.min.r.female))
    
    if(length(remp.temp) > 1 & !is.null(xover.min.r.female) & min(diff(r.cumu.female[remp.temp])) < xover.min.r.female){
      
      while(!is.null(xover.min.r.female) & min(diff(r.cumu.female[remp.temp])) < xover.min.r.female) {
        remp.temp <- which(((runif(length(r.female)) < r.female) + 0L) == 1)
        print(remp.temp)
      }
      
    }
      
      
  

  
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
  
head(x)
