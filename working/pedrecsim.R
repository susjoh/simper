# Simulate recombination events in R


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = F)

setwd("C:/Users/Susan Johnston/Desktop/PedRecSim/")

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

#~~ Recode founder gametes to 0 e.g. if one parent is unknown

transped$Cohort[which(transped$Parent.ID == 0)] <- 0

#~~ Redefine columns

names(transped) <- c("Offspring.ID", "Parent.ID.SEX", "Parent.ID", "Cohort")
transped$Key <- paste(transped$Parent.ID, transped$Offspring.ID, sep = "_")


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

transped.nonfounders <- transped[which(transped$Cohort != 0),]

template.list <- sapply(transped$Parent.ID.SEX, function(x){
  if(x == "MOTHER") remp.temp <- which(((runif(length(female.r)) < female.r) + 0L) == 1)
  if(x == "FATHER") remp.temp <- which(((runif(length(male.r  )) < male.r  ) + 0L) == 1)
  remp.temp
})

#~~ Add key to list

names(template.list) <- transped.nonfounders$Key

#~~ Calculate recombination count

transped.nonfounders$RecombCount <- unlist(lapply(template.list, length))

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




transped[j,]
haplo.list[transped$Parent.ID[j]]

#~~ Pull out recombination information

rec.pos <- template.list[transped$Key[j]][[1]]

#~~ If there are no recombination events, sample one of the parental haplotypes at random

if(length(rec.pos) == 0){
  haplo.list[transped$Offspring.ID[j]][[1]][transped$Parent.ID.SEX[j]] <- haplo.list[transped$Parent.ID[j]][[1]][sample.int(2, 1)]
}

if


  
  
  sample.int(2, 1)
  
  library(rbenchmark)
  benchmark((runif(1) < 0.5) + 0L, 
            sample.int(2, 1),
            sample(c(0, 1), 1), replications = 100000)


haplo.list[transped$Parent.ID[j]][[1]][sample.int(2, 1)]


transped$Key[which(transped$Offspring.ID == transped$Parent.ID[j])]


which(transped$Cohort == 1)



for(i in 1:max(simped2.chr$Cohort)){
  
  print(paste("Analysing cohort", i, "of", max(simped2.chr$Cohort)))
  
  cohortped <- subset(simped2.chr, Cohort == i & RRID != 0)
  
  if(length(which(cohortped$HaploKey %in% names(haplo.list))) > 0) break("Already has Haplo Information")
  
  temp.haplo.list <- mapply(CreateHaplotypes, cohortped$Key, cohortped$ParentKey1, cohortped$ParentKey2, SIMPLIFY = F)
  
  names(temp.haplo.list) <- cohortped$HaploKey
  
  haplo.list <- c(haplo.list, temp.haplo.list)
  
  rm(cohortped, temp.haplo.list)
}














transped[which(transped$Cohort == 0),]




simped2.chr <- subset(simped2, Chr == chrsim)

#~~ Run function on simped2 with cohort == 0 to create initial haplo.list # 6 seconds

haplo.list <- lapply(simped2.chr$Chr[which(simped2.chr$RRID == 0)], FounderHaploSamp, lmap = lmap)
names(haplo.list) <- simped2.chr$HaploKey[which(simped2.chr$RRID == 0)]

#~~ generate a temporary key

simped2.chr$ParentKey1 <- paste(simped2.chr$RRID, simped2.chr$Chr, "MOTHER", sep = "_")
simped2.chr$ParentKey2 <- paste(simped2.chr$RRID, simped2.chr$Chr, "FATHER", sep = "_")

#~~ Now loop the function over each cohort

for(i in 1:max(simped2.chr$Cohort)){
  
  print(paste("Analysing cohort", i, "of", max(simped2.chr$Cohort)))
  
  cohortped <- subset(simped2.chr, Cohort == i & RRID != 0)
  
  if(length(which(cohortped$HaploKey %in% names(haplo.list))) > 0) break("Already has Haplo Information")
  
  temp.haplo.list <- mapply(CreateHaplotypes, cohortped$Key, cohortped$ParentKey1, cohortped$ParentKey2, SIMPLIFY = F)
  
  names(temp.haplo.list) <- cohortped$HaploKey
  
  haplo.list <- c(haplo.list, temp.haplo.list)
  
  rm(cohortped, temp.haplo.list)
}

save(recped.chr, simped2.chr, haplo.list, template.list, file = paste0("simulation_chr_", chrsim, "_", analysisid, ".Rdata"))

table(unlist(lapply(haplo.list, function(x) length(x))))
table(unlist(lapply(template.list, function(x) length(x))))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Convert to precrimap format to create crimap files          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ introduce error to the haplo.list

dataset.sim <- length(haplo.list) * length(haplo.list[[1]])

no.errors <- round(dataset.sim*error.rate)

error.samp.key <- sample(1:length(haplo.list),      replace = F,   size = no.errors)
error.samp.loc <- sample(1:length(haplo.list[[1]]), replace = T,   size = no.errors)

error.tab <- rbind(error.tab, data.frame(Key = names(haplo.list)[error.samp.key],
                                         LocusOrder = error.samp.loc,
                                         LocusID = lmap$SNP.Name[which(lmap$Chr == chrsim)][error.samp.loc],
                                         Chr = chrsim))

for(i in 1:no.errors){
  x <- haplo.list[[error.samp.key[i]]][error.samp.loc[i]]
  haplo.list[[error.samp.key[i]]][error.samp.loc[i]] <- ifelse(x == 1, 2, 1)
}

table(unlist(lapply(haplo.list, function(x) length(x))))

#~~ create missing data information for next stage

haplo.list2 <- haplo.list

haplo.list <- haplo.list2

no.missing <- round(dataset.sim*missing.rate/2)

missing.samp.key <- sample(1:length(haplo.list),      replace = T, size = no.missing)
missing.samp.loc <- sample(1:length(haplo.list[[1]]), replace = T, size = no.missing)

missing.tab <- rbind(missing.tab, data.frame(Key = names(haplo.list)[missing.samp.key],
                                             LocusOrder = missing.samp.loc,
                                             LocusID = lmap$SNP.Name[which(lmap$Chr == chrsim)][missing.samp.loc],
                                             Chr = chrsim))

head(missing.tab)

for(i in 1:no.missing) haplo.list[[missing.samp.key[i]]][missing.samp.loc[i]] <- 0

table(unlist(lapply(haplo.list, function(x) length(x))))

#~~ need to create a genotype for each individual in the famped used in CRIMAP

crimprep1 <- simped2.chr[,c("Offspring.ID", "HaploKey")]

#~~ Shuffle the rows

crimprep1 <- crimprep1[sample(nrow(crimprep1)),]

#~~ Cast the table to create a single line for each Offspring.ID

crimprep2 <- cast(crimprep1, Offspring.ID ~ ., function(x) x)

#~~ Merge with family information

names(crimprep2)[1] <- "ANIMAL" 

crimprep2 <- join(famped, crimprep2)

crimprep2 <- crimprep2[,c("Family", "ANIMAL", "MOTHER", "FATHER",  "SEX", "X1", "X2")]

#~~ Create genotype files based on the haplotypes

CreateGenoLineFunc <- function(x, y){
  
  miss.loc <- unique(c(which(haplo.list[[x]] == 0), which(haplo.list[[y]] == 0)))
  
  if(length(miss.loc) > 0)  {
    haplo.list[[x]][miss.loc] <- 0
    haplo.list[[y]][miss.loc] <- 0
  }
  
  paste(paste(haplo.list[[x]], haplo.list[[y]]), collapse = " ")
}

crimprep2$Geno <- as.vector(mapply(CreateGenoLineFunc, crimprep2$X1, crimprep2$X2))

table(sapply(crimprep2$Geno, nchar))
#~~ remove the X1 and X2 columns

crimprep2 <- crimprep2[,c("Family", "ANIMAL", "MOTHER", "FATHER",  "SEX", "Geno")]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Create crimap files                                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

crimprep2 <- as.matrix(crimprep2)

print(paste("Creating Crimap File for Chromosome", chrsim, "of 26"))

system("cmd", input = paste0("del chr", chrsim, "_", analysisid, "*"))

nfamilies <- length(unique(crimprep2[,"Family"]))  
nloci <- length(which(lmap$Chr == chrsim))
locus.names <- lmap$SNP.Name[which(lmap$Chr == chrsim)]

outfile <- paste("chr", chrsim, "_", analysisid, ".gen", sep = "")

write.table(nfamilies,   outfile, row.names = F, quote = F, col.names = F)
write.table(nloci,       outfile, row.names = F, quote = F, col.names = F, append=T)
write.table("",          outfile, row.names = F, quote = F, col.names = F, append=T)
write.table(locus.names, outfile, row.names = F, quote = F, col.names = F, append=T)
write.table("",          outfile, row.names = F, quote = F, col.names = F, append=T)

counter <- 0

for(j in unique(crimprep2[,"Family"])){
  
  counter <- counter + 1
  if(i %in% counter %in% seq(1, length(unique(crimprep2[,"Family"])), 100)){ 
    print(paste("Analysing Family", counter, "of", length(unique(crimprep2[,"Family"]))))
  }
  
  famtab <- crimprep2[which(crimprep2[,"Family"] == j),-1]
  famsize <- nrow(famtab)
  
  write.table(j,       outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table(famsize, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("",      outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table(famtab,  outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("",      outfile, row.names = F, quote = F, col.names = F, append=T)
  
}

rm(famtab, famsize, nfamilies, nloci, locus.names, outfile)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Run Crimap and Chrompic                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)

system("cmd", input = paste0("\"../crimap2504.exe\" ", chrsim, "_", analysisid, " prepare < crimapinput1"))

system("cmd", input = paste0("\"../crimap2504.exe\" ", chrsim, "_", analysisid, " chrompic > chr", chrsim, "_", analysisid, ".cmp"))

system("cmd", input = "del *.cg")

}

write.table(missing.tab, paste0("simulation_", analysisid, "_missingsim.txt"))
write.table(error.tab,   paste0("simulation_", analysisid, "_errorsim.txt"))
}


#~~ Add Recomb Count information if crimap was run previously

if(length(which(is.na(recped$RecombCount))) == nrow(recped)){
  
  rec.info <- NULL
  
  for(chrsim in 1:26){
    
    load(paste0("simulation_chr_", chrsim, "_", analysisid, ".Rdata"))
    
    rec.info <- rbind(rec.info, recped.chr[,c("Key", "RecombCount")])
    
    rm(recped.chr, simped2.chr, haplo.list, template.list)
  }
  
  recped <- subset(recped, select = -RecombCount)
  
  recped <- join(recped, rec.info)
  
  rm(rec.info)  
}


save(recped, lmap, famped, AnalysisSuffix, analysisid, file = "../../results/5.1_Post_Simulation_Input.Rdata")


















# lmap <- read.table("C:/Users/Susan Johnston/Desktop/Recombination Rate Study/results/2_Merged_map_f.txt", header = T)
# head(lmap)
# lmap <- subset(lmap, Chr == 24)
# 
# mapfile <- read.table("map.txt", header = T)
# head(mapfile)
# mapfile <- mapfile[, c(1, 4)]
# 
# mapfile <- join(mapfile, lmap)
# 
# snp.names <- mapfile$SNP.Name
# snp.mafs  <- mapfile$MAF
# female.r  <- mapfile$Female.r
# male.r    <- mapfile$Male.r
# female.r  <- female.r[1:(length(female.r) - 1)]
# male.r    <- male.r[1:(length(male.r) - 1)]
# 
# 
# save(snp.names, snp.mafs, female.r, male.r, ped, file = "test.Rdata")
