# #~~ Sample founders
# 
# transped$Allele <- NA
# transped$Allele[which(transped$Cohort == 0)] <- replicate(length(which(transped$Cohort == 0)), (runif(length(founder.allele.freq)) < founder.allele.freq) + 0L)
# 
# #~~ sample cohort by cohort
# 
# for(cohort in 1:max(transped$Cohort)){
#   
#   print(paste("Simulating haplotypes for cohort", cohort, "of", max(transped$Cohort)))
#   
#   for(j in which(transped$Cohort == cohort)){
#     
#     transped$Allele <- transped[which(transped$Offspring.ID == transped[j,]$Parent.ID),]$Allele[sample.int(2, 1)]
#     
#   }
# }
# 
# 
# 


