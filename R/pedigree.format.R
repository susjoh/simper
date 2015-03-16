#' Format Pedigree
#'
#' This function formats the pedigree for downstream analysis.
#' @param ped Pedigree object in "simple" format (Three columns for ANIMAL, 
#' MOTHER and FATHER) or in "plink" format (Five to Six columns for FAMILY,
#' ANIMAL, FATHER, MOTHER, SEX and Phenotype, where the phenotype column is
#' optional). The simple argument can recognise the order of parents if
#' they are named sensibly. Run simple.ped.name.rules() for an example. 
#' @param pedigree.type Defaults to "simple", can also accept "plink" which
#' is equivalent for for first 5 to 6 columns of a PLINK .ped file.
#' @keywords
#' @export
#'
#' 


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
    if(!all(names(ped) %in% c("ID", "ANIMAL", "MUM", "MOM", 
                             "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE"))) stop(simple.ped.name.rules())

    names(ped)[which(names(ped) %in% c("ID", "ANIMAL"))] <- "ANIMAL"
    names(ped)[which(names(ped) %in% c("MUM", "MOM", "MOTHER", "DAM"))] <- "MOTHER"
    names(ped)[which(names(ped) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"
    
    ped <- ped[,c("ANIMAL", "MOTHER", "FATHER")]
    
    for(i in 1:3) ped[which(is.na(ped[,i])),i] <- 0
    
    if(any(!ped$MOTHER %in% ped$ANIMAL)){
    ped <- rbind(data.frame(ANIMAL = ped$MOTHER[which(!ped$MOTHER %in% ped$ANIMAL)],
                            MOTHER = 0, FATHER = 0),
                 ped)
    }
    
    if(any(!ped$FATHER %in% ped$ANIMAL)){
      ped <- rbind(data.frame(ANIMAL = ped$FATHER[which(!ped$FATHER %in% ped$ANIMAL)],
                            MOTHER = 0, FATHER = 0),
                   ped)
    }
    
    ped <- subset(ped, ANIMAL != 0)
    ped <- droplevels(unique(ped))
    
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
      stop("Number of columns does not match those expected of PLINK format.")
    }
    
    # re-code missing values
    
    for(i in 2:4) ped[which(is.na(ped[,i])),i] <- 0
    for(i in 2:4) ped[which(ped[,i] == -9),i]  <- 0
    
  }
  
  ped
  
}


#' Extract a vector of founder individuals from a pedigree
#'
#' This function extracts a vector of founder IDs from a pedigree object. 
#' A founder is an individual that does not have parents within the pedigree.
#' @param ped Pedigree object in "simple" format (Three columns for ANIMAL, 
#' MOTHER and FATHER) or in "plink" format (Five to Six columns for FAMILY,
#' ANIMAL, FATHER, MOTHER, SEX and Phenotype, where the phenotype column is
#' optional). The simple argument can recognise the order of parents if
#' they are named sensibly. Run simple.ped.name.rules() for an example. 
#' @param pedigree.type Defaults to "simple", can also accept "plink" which
#' is equivalent for for first 5 to 6 columns of a PLINK .ped file.
#' @keywords
#' @export
#'
#' 


founderIDs <- function(pedigree, pedigree.type = "simple"){
  
  ped <- pedigree.format(pedigree, pedigree.type = pedigree.type)
  
  cohorts <- data.frame(ANIMAL = ped[,1],
                        Cohort = kindepth(ped[,1], ped[,3], ped[,2]))
  
  cohorts$ANIMAL[which(cohorts$Cohort == 0)]
}



#' Format Pedigree for genedrop analysis
#'
#' This function formats the pedigree for downstream analysis.
#' @param ped Pedigree object. Run simple.ped.name.rules() for an example. 
#' @keywords
#' @export
#'
#' 

pedigree.format.genedrop <- function(ped){
  names(ped) <- toupper(names(ped))
  if(!names(ped)[1] %in% c("ID", "ANIMAL"))                 stop(simple.ped.name.rules())
  if(!names(ped)[2] %in% c("MUM", "MOM", "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE")) stop(simple.ped.name.rules())
  if(!names(ped)[3] %in% c("MUM", "MOM", "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE")) stop(simple.ped.name.rules())
  
  names(ped)[which(names(ped) %in% c("ID", "ANIMAL"))] <- "ANIMAL"
  names(ped)[which(names(ped) %in% c("MUM", "MOM", "MOTHER", "DAM"))] <- "MOTHER"
  names(ped)[which(names(ped) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"
  
  ped <- ped[,c("ANIMAL", "MOTHER", "FATHER")]
  ped$Cohort <- cohort
  
  for(i in which(names(ped) %in% c("ANIMAL", "MOTHER", "FATHER"))) ped[which(is.na(ped[,i])),i] <- 0
  
  #~~ ped must have a line for all mothers and fathers in ID
  
  if(any(!ped$MOTHER %in% ped$ANIMAL & !ped$MOTHER == 0)) stop("MOTHER ids missing from ID")
  if(any(!ped$FATHER %in% ped$ANIMAL & !ped$FATHER == 0)) stop("FATHER ids missing from ID")
  
  ped
}
