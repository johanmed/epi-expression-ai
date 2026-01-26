#Script to align strains in genotype and phenotype file for GEMMA

# Pass names of BXD lines
new_names <- c("marker", "allele1", "allele2", "BXD1", "BXD2", "BXD5", "BXD6", "BXD8", "BXD9", "BXD11", "BXD12", "BXD13", "BXD14", "BXD15", "BXD16", "BXD18", "BXD19", "BXD20", "BXD21", "BXD22", "BXD23", "BXD24", "BXD24a", "BXD25", "BXD27", "BXD28", "BXD29", "BXD30", "BXD31", "BXD32", "BXD33", "BXD34", "BXD35", "BXD36", "BXD37", "BXD38", "BXD39", "BXD40", "BXD41", "BXD42", "BXD43", "BXD44", "BXD45", "BXD48", "BXD48a", "BXD49", "BXD50", "BXD51", "BXD52", "BXD53", "BXD54", "BXD55", "BXD56", "BXD59", "BXD60", "BXD61", "BXD62", "BXD63", "BXD64", "BXD65", "BXD65a", "BXD65b", "BXD66", "BXD67", "BXD68", "BXD69", "BXD70", "BXD71", "BXD72", "BXD73", "BXD73a", "BXD73b", "BXD74", "BXD75", "BXD76", "BXD77", "BXD78", "BXD79", "BXD81", "BXD83", "BXD84", "BXD85", "BXD86", "BXD87", "BXD88", "BXD89", "BXD90", "BXD91", "BXD93", "BXD94", "BXD95", "BXD98", "BXD99", "BXD100", "BXD101", "BXD102", "BXD104", "BXD105", "BXD106", "BXD107", "BXD108", "BXD109", "BXD110", "BXD111", "BXD112", "BXD113", "BXD114", "BXD115", "BXD116", "BXD117", "BXD119", "BXD120", "BXD121", "BXD122", "BXD123", "BXD124", "BXD125", "BXD126", "BXD127", "BXD128", "BXD128a", "BXD130", "BXD131", "BXD132", "BXD133", "BXD134", "BXD135", "BXD136", "BXD137", "BXD138", "BXD139", "BXD141", "BXD142", "BXD144", "BXD145", "BXD146", "BXD147", "BXD148", "BXD149", "BXD150", "BXD151", "BXD152", "BXD153", "BXD154", "BXD155", "BXD156", "BXD157", "BXD160", "BXD161", "BXD162", "BXD165", "BXD168", "BXD169", "BXD170", "BXD171", "BXD172", "BXD173", "BXD174", "BXD175", "BXD176", "BXD177", "BXD178", "BXD180", "BXD181", "BXD183", "BXD184", "BXD186", "BXD187", "BXD188", "BXD189", "BXD190", "BXD191", "BXD192", "BXD193", "BXD194", "BXD195", "BXD196", "BXD197", "BXD198", "BXD199", "BXD200", "BXD201", "BXD202", "BXD203", "BXD204", "BXD205", "BXD206", "BXD207", "BXD208", "BXD209", "BXD210", "BXD211", "BXD212", "BXD213", "BXD214", "BXD215", "BXD216", "BXD217", "BXD218", "BXD219", "BXD220", "C57BL/6JxBXD037F1", "BXD001xBXD065aF1", "BXD009xBXD170F1", "BXD009xBXD172F1", "BXD012xBXD002F1", "BXD012xBXD021F1", "BXD020xBXD012F1", "BXD021xBXD002F1", "BXD024xBXD034F1", "BXD032xBXD005F1", "BXD032xBXD028F1", "BXD032xBXD65bF1", "BXD034xBXD024F1", "BXD034xBXD073F1", "BXD055xBXD074F1", "BXD055xBXD65bF1", "BXD061xBXD071F1", "BXD062xBXD077F1", "BXD065xBXD077F1", "BXD069xBXD090F1", "BXD071xBXD061F1", "BXD073bxBXD065F1", "BXD073bxBXD077F1", "BXD073xBXD034F1", "BXD073xBXD065F1", "BXD073xBXD077F1", "BXD074xBXD055F1", "BXD077xBXD062F1", "BXD083xBXD045F1", "BXD087xBXD100F1", "BXD065bxBXD055F1", "BXD102xBXD077F1", "BXD102xBXD73bF1", "BXD170xBXD172F1", "BXD172xBXD197F1", "BXD197xBXD009F1", "BXD197xBXD170F1")

# Load genotype file and name appropriately
geno <- read.table(file="data/BXD.8_geno.txt",
                   col.names=new_names,
                   sep=",",
                   strip.white = TRUE)

# Load phenotype file and extract names of BXD lines
pheno <- read.csv("~/interest_methylation_phenofile_strains.csv", header=TRUE)
strains <- pheno$X

# Order genotypes based on strains extracted from phenotype file
new_geno <- geno[, c("marker", "allele1", "allele2", strains)]

# Save new genotype file
write.table(new_geno,
            file="~/new_BXD.8_geno.txt",
            sep=", ",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
