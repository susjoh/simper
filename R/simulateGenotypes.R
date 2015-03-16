
#' Simulate recombination events in R
#' 
#' This function simulates genotypes within a pedigree, based on centiMorgan map
#' distances. Outputs genotypes as a GenABEL gwaa object.
#' @param ped data.frame of pedigree in "simple" format (Three columns for 
#'   ANIMAL, MOTHER and FATHER) or in "plink" format (Five to Six columns for 
#'   FAMILY, ANIMAL, FATHER, MOTHER, SEX and Phenotype, where the phenotype 
#'   column is optional). Is formatted within this function using 
#'   pedigree.format()
#' @param pedigree.type Character. "simple" or "plink"
#' @param cM Numeric vector of sex-averaged map distances in centiMorgans. Not 
#'   required if sex specific map distances are defined. First position should 
#'   be given as 0cM
#' @param cM.male Numeric vector of Male map distances in centiMorgans
#' @param cM.female Numeric vector of Female map distances in centiMorgans
#' @param founder.mafs Numeric vector of Minor allele frequencies for each 
#'   marker
#' @param chromosome.ids Vector of chromosome or linkage group IDs for each 
#'   marker
#' @param snp.names Optional character vector of SNP names
#' @param map.distances Optional vector of map distances in base pairs for 
#'   GenABEL object. Not used in the simulation.
#' @param error.rate numeric. Default is 1e-4. Future versions will allow vector.
#' @param missing.rate numeric. Default is 0.001. Future versions will allow vector.
#' @param xover.min.cM numeric. Minimum cM distance between two crossovers
#' @param xover.min.cM.male numeric. Minimum cM distance between two crossovers in males
#' @param xover.min.cM.female numeric. Minimum cM distance between two crossovers in females



simulateGenos <- function(ped,
                          pedigree.type = "simple",
                          cM = NULL,
                          cM.male = NULL,
                          cM.female = NULL,
                          founder.mafs = NULL,
                          chromosome.ids = NULL,
                          snp.names = NULL,
                          map.distances = NULL,
                          error.rate = 1e-4,
                          missing.rate = 0.001,
                          xover.min.cM = NULL,
                          xover.min.cM.male = NULL,
                          xover.min.cM.female = NULL){
  #   
  #   cM        = c(0,  5, 10,  15,   25,   35,   45,   55,   65,   65, 0,  5, 10,  15,   25,   35,   45,   55,   65,   65)
  #   cM.male   =  c(0, 10, 20,  30,   50,   70,   90,  110,  130,  130, 0, 10, 20,  30,   50,   70,   90,  110,  130,  130)
  #   cM.female =  c(0, 2.5,  5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 32.5, 0, 2.5,  5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 32.5)
  #   chromosome.ids <- rep(c(23, 24), each = 10)
  #   snp.names     <- NULL
  #   founder.mafs  <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.2)
  #   map.distances <- NULL
  #   xover.min.cM            <- 25
  #   xover.min.cM.male       <- NULL
  #   xover.min.cM.female     <- NULL                             
  #   error.rate = 1e-4
  #   missing.rate = 0.001  
  #   ped <- read.table("../../Recombination Analysis G/data/pedigree_20130920_13.txt", header = T)  
  #   pedigree.type = "simple"
  
  
  require(plyr)
  require(reshape2)
  require(GenABEL)
  require(kinship2)
  
  phased.output <-  TRUE
  
  
  #~~ Format the pedigree
  
  ped <- pedigree.format(ped, pedigree.type = pedigree.type)
  
  #~~ melt pedfile to get a unique row for each gamete transfer
  
  transped <- melt(ped[,c("ANIMAL", "MOTHER", "FATHER")], id.vars = "ANIMAL")
  transped$variable <- as.character(transped$variable)
  
  #~~ assign pedfile IDs to cohort and merge with transped
  
  cohorts <- data.frame(ANIMAL = ped[,1],
                        Cohort = kindepth(ped[,1], ped[,3], ped[,2]))
  
  suppressMessages(transped <- join(transped, cohorts))
  
  #~~ Redefine columns
  
  names(transped) <- c("Offspring.ID", "Parent.ID.SEX", "Parent.ID", "Cohort")
  transped$Key    <- paste(transped$Parent.ID, transped$Offspring.ID, sep = "_")
  
  #~~ Recode founder gametes to 0 e.g. if one parent is unknown
  
  transped$Cohort[which(transped$Parent.ID == 0)] <- 0
  
  
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
  
  message("Assuming 1cM is equivalent to r = 0.01")
  
  #~~ Tidy up the SNP names and founder MAFs etc. if not defined
  
  if(is.null(snp.names))  snp.names <- paste0("SNP", 1:length(cM.male))
  if(is.null(founder.mafs)){
    message("Founder MAFs have not been defined - assuming MAF = 0.5")
    founder.mafs <- rep(0.5, times = length(cM.male))
  }
  if(is.null(map.distances)){
    message("Base pair distances have not been defined - Dummy distances used")
    map.distances <- as.integer(cM.female*1e6)
  }
  
  if(is.null(chromosome.ids)){
    message("Assuming all markers are on the same chromosome. Errors may occur when undefined")
    chromosome.ids <- rep(1, length(cM.female))
  }
  
  
  #~~ RUN THE SIMULATIONS

  
  #~~ Convert to R
  

  r.male   <- diff(cM.male/100)
  r.female <- diff(cM.female/100)
  
  r.male[which(r.male < 0)] <- 0.5
  r.female[which(r.female < 0)] <- 0.5
  
  
  #~~ Convert to r
  
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
  names(haplo.list) <- unique(as.character(transped$Offspring.ID))
  
  #~~ Sample the founder haplotypes
  
  system.time(for(i in which(transped$Cohort == 0)){
    if(transped$Parent.ID.SEX[i] == "MOTHER") haplo.list[as.character(transped$Offspring.ID[i])][[1]]$MOTHER <- (runif(length(founder.mafs)) < founder.mafs) + 0L
    if(transped$Parent.ID.SEX[i] == "FATHER") haplo.list[as.character(transped$Offspring.ID[i])][[1]]$FATHER <- (runif(length(founder.mafs)) < founder.mafs) + 0L
  })
  
  
  #~~ Sample the non-founder haplotypes by cohort. This loops through cohorts sequentially as
  #   parental haplotypes must exist before sampling.  
  
  for(cohort in 1:max(transped$Cohort)){
    
    print(paste("Simulating haplotypes for cohort", cohort, "of", max(transped$Cohort)))
    
    for(j in which(transped$Cohort == cohort)){
      
      #~~ Pull out recombination information from template.list
      
      rec.pos <- template.list[transped$Key[j]][[1]]
      
      #~~ If there are no recombination events, sample one of the parental haplotypes at random
      
      if(length(rec.pos) == 0){
        haplo.list[as.character(transped$Offspring.ID[j])][[1]][transped$Parent.ID.SEX[j]] <- haplo.list[as.character(transped$Parent.ID[j])][[1]][sample.int(2, 1)]
      }
      
      #~~ If there are recombination events, sample the order of the parental haplotypes,
      #   exchange haplotypes and then transmit haplotype to offspring
      
      if(length(rec.pos) > 0){
        
        parental.haplotypes <- haplo.list[as.character(transped$Parent.ID[j])][[1]][sample.int(2, 2, replace = F)]
        
        start.pos <- c(1, rec.pos[1:(length(rec.pos))] + 1)
        stop.pos <- c(rec.pos, length(r.female) + 1)
        
        fragments <- list()
        
        for(k in 1:length(start.pos)){
          if(k %% 2 != 0) fragments[[k]] <- parental.haplotypes[1][[1]][start.pos[k]:stop.pos[k]]
          if(k %% 2 == 0) fragments[[k]] <- parental.haplotypes[2][[1]][start.pos[k]:stop.pos[k]]  
        }
        
        haplo.list[as.character(transped$Offspring.ID[j])][[1]][transped$Parent.ID.SEX[j]] <- list(unlist(fragments))
        
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
  
  colnames(x) <- rep(snp.names, each = 2)
  
  x <- x + 1
  
  #~~ Create errors in the dataset
  
  geno.count <- nrow(x) * (ncol(x)/2)
  error.count <- as.integer(geno.count * error.rate)
  missing.count <- as.integer(geno.count * missing.rate)
  
  if(error.count > 0){
    
    errortemp <- data.frame(row = sample(1:nrow(x), replace = T, size = error.count),
                            col = sample(1:ncol(x), replace = T, size = error.count))
    
    
    errortemp$col2 <- ifelse(errortemp$col %% 2 == 0, errortemp$col - 1, errortemp$col + 1)
    
    for(i in 1:nrow(errortemp)){
      tempgeno <- paste(x[errortemp$row[i],c(errortemp$col[i], errortemp$col2[i])], collapse = "")
      if(tempgeno == "11")            newgeno <- sample(c("12", "22"), size = 1) 
      if(tempgeno %in% c("12", "21")) newgeno <- sample(c("11", "22"), size = 1) 
      if(tempgeno == "22")            newgeno <- sample(c("11", "12"), size = 1)
      x[errortemp$row[i],c(errortemp$col[i], errortemp$col2[i])] <- as.vector(strsplit(newgeno, split = "")[[1]])
      
    }
  }
  
  
  if(missing.count > 0){
    
    missingtemp <- data.frame(row = sample(1:nrow(x), replace = T,  size = missing.count),
                              col = sample(1:ncol(x), replace = T,  size = missing.count))
    
    
    missingtemp$col2 <- ifelse(missingtemp$col %% 2 == 0, missingtemp$col - 1, missingtemp$col + 1)
    
    missingtemp <- data.frame(rbind(as.matrix(missingtemp[,c(1, 2)]),
                                    as.matrix(missingtemp[,c(1, 3)])))
    
    for(i in 1:nrow(missingtemp)) x[missingtemp$row[i],missingtemp$col[i]] <- 0
  }
  
  #~~ Create the return object in PLINK format
  
  newx <- data.frame(FAMILY = 1, ANIMAL = row.names(x))
  suppressMessages(newx <- join(newx, ped))
  newx$SEX <- ifelse(newx$ANIMAL %in% newx$FATHER, 1, 0)
  newx$PHENO <- -9
  
  x <- cbind(newx, x)
  
  
  mapobj <- data.frame(Chromosome = chromosome.ids,
                       SNP.Name = snp.names,
                       bpPosition = map.distances)
  
  genabel.phenotab <- x[,c("ANIMAL", "SEX")]
  names(genabel.phenotab) <- c("id", "sex")
  
  #~~ generate a temporary file name
  
  tempfile <- paste(sample(letters, replace = T, size = 10), collapse = "")
  
  write.table(x,      paste0(tempfile, ".ped"), row.names = F, col.names = F, quote = F)
  write.table(mapobj, paste0(tempfile, ".map"), row.names = F, col.names = F, quote = F)
  write.table(genabel.phenotab, paste0(tempfile, ".pheno"), row.names = F, quote = F)
  
  
  convert.snp.ped(pedfile = paste0(tempfile, ".ped"), 
                  mapfile = paste0(tempfile, ".map"),
                  outfile = paste0(tempfile, ".genabel"),
                  traits = 1, strand = "u",  mapHasHeaderLine = F, bcast = F)
  
  
  genabel.geno <- load.gwaa.data(genofile = paste0(tempfile, ".genabel"), 
                                 phenofile = paste0(tempfile, ".pheno"))

  
  system("cmd", input = paste0("del ", tempfile, "*"), show.output.on.console = F)
  
  # return object
  
  return(list(genabel.object = genabel.geno,
              templates = template.list))
}

