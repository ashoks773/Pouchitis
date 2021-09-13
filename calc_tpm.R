#Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.


#--- Load filtered file which was Generated in Step1_RNA_Seq
TABLE <- read.csv(file = "ProCoding_Gene_counts_mitoremove_Riboremove_Filtered.txt", sep = "\t", header = TRUE, row.names = 1)
TABLE_FILTERED <- subset(TABLE, MappabilityLength != 0) #--- 12 Genes have mappability length 0 so they need to be removed.
len <- TABLE_FILTERED$MappabilityLength

#--- Load function to calculate TPM values
r_tpm <- function(dfr,len)
{
  dfr1 <- sweep(dfr,MARGIN=1,(len/10^4),`/`)
  scf <- colSums(dfr1)/(10^6)
  return(sweep(dfr1,2,scf,`/`))
}

TABLE_FILTERED_TPM <- r_tpm(TABLE_FILTERED[,5:98], len)
TABLE_FILTERED_TPM <- cbind(TABLE_FILTERED[,2:4], TABLE_FILTERED_TPM) #-- No need to Include cumulative counts which was not for TPM
write.table(TABLE_FILTERED_TPM, file="AllSamples_Filtered_TPM.txt", sep="\t")

#--Other function can be used
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#----Tried but didnot work
len_normalized <- len/1000
TABLE_RPK <- sapply(TABLE_FILTERED[,c(5:98)], function(x){(x/len_normalized)})
PerM_ScalingF <- colSums(TABLE_RPK, na.rm = TRUE)/1000000
TABLE_TPM <- TABLE_RPK/PerM_ScalingF
FINAL_TABLE_TPM <- cbind(TABLE_FILTERED[,c(2:4)], TABLE_TPM)
