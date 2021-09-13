setwd("~/Box/David_Casero_Lab/Pouchitis_Project")
#- We will only load antisense counts - which are higher than Sense counts
Pouchitis_rnaCnt = read.table("POUCHITISrna_antisense.All.cnt")
dim (Pouchitis_rnaCnt)
Pouchitis_samples<-read.table("POUCHITISrna_antisenseSamples.txt")
dim (Pouchitis_samples)
Gencode_33_Selected_Geneid<-read.table("countsannot_GRCh38.primary.Selected.Geneid.txt")
dim(Gencode_33_Selected_Geneid)

Gencode_33_Selected_MappSS<-read.table("gencode.v33.Selected.ReadsPerGene.out.MappSS.txt")
Gencode_33_Selected_MappUS<-read.table("gencode.v33.Selected.ReadsPerGene.out.MappUS.txt")
Gencode_33_Selected_Biotype<-read.table("gencode.v33.annotation.Selected.biotype.txt")
Gencode_33_Selected_Genename<-read.table("gencode.v33.annotation.Selected.genename.txt")

Pouchitis_rnaCnt1<-cbind(as.data.frame(Gencode_33_Selected_Genename[,1]),Pouchitis_rnaCnt)

Pouchitis_rnaTPM = Pouchitis_rnaCnt/Gencode_33_Selected_MappSS[,1]*1000
Pouchitis_rnaTPM1<-apply(Pouchitis_rnaTPM,2,as.numeric)
Pouchitis_rnaTPM1[is.nan(Pouchitis_rnaTPM1)]<-0
#j<-which(is.infinite(IGF2BP1ripTPM1))
Pouchitis_rnaTPM1[is.infinite(Pouchitis_rnaTPM1)]<-0
Pouchitis_rnaTPM<-Pouchitis_rnaTPM1

sumrow<-apply(Pouchitis_rnaTPM1,2,sum)
Pouchitis_rnaTPM2 = Pouchitis_rnaTPM1/sumrow*1000000
colnames(Pouchitis_rnaTPM2)<-Pouchitis_samples[,1]
row.names(Pouchitis_rnaTPM2)<-Gencode_33_Selected_Genename[,1]
rn<-nrow(Pouchitis_rnaTPM2)
m<-unique(sample(rn, 1000))
rand_data<-Pouchitis_rnaTPM2[m,]
d <- dist(t(rand_data), method = "euclidean")
for(i in 1:9999){
  m<-unique(sample(rn, 1000))
  rand_data<-Pouchitis_rnaTPM2[m,]
  d <-d+ dist(t(rand_data), method = "euclidean")	
  
}

d1<-d/10000
fit <- hclust(sqrt(d1), method="ward.D2") 
jpeg("Random_1k_Genes.jpg", height = 8, width = 16, units = 'in', res = 600)
plot(fit)
dev.off ()

######################
#----- Perform Masking
######################
allbiotypes = unique(Gencode_33_Selected_Biotype)

#allbiotypes<-unique(Gencode_33_Selected_Biotype)
allbiotypeslength <- rep(NA,nrow(allbiotypes))
allbiotypescountsPouchitis_rna <- matrix(0,nrow(allbiotypes),ncol(Pouchitis_rnaTPM2))
allbiotypescountspercentsPouchitis_rna = matrix(0,nrow(allbiotypes),ncol(Pouchitis_rnaTPM2))

for (i in 1:nrow(allbiotypes)){
  
  temp = is.element(Gencode_33_Selected_Biotype[,1],allbiotypes[i,1])
  
  allbiotypeslength[i] = length(which(temp==T))
  
  if (allbiotypeslength[i]>1){
    ind = is.element(Gencode_33_Selected_Biotype[,1],allbiotypes[i,1])	
    subdata<-Pouchitis_rnaCnt[ind,]
    allbiotypescountsPouchitis_rna[i,] = apply(subdata,2,sum) 
    allbiotypescountspercentsPouchitis_rna[i,] = allbiotypescountsPouchitis_rna[i,]/apply(Pouchitis_rnaCnt,2,sum)*100; 
    
    rnasubdata<-Pouchitis_rnaCnt[ind,]
    allbiotypescountsPouchitis_rna[i,] = apply(rnasubdata,2,sum) 
    allbiotypescountspercentsPouchitis_rna[i,] = allbiotypescountsPouchitis_rna[i,]/apply(Pouchitis_rnaCnt,2,sum)*100; 
    
  }
}
allbiotypescountspercentsPouchitis_rna[c(2,3),]
allbiotypescountspercentsPouchitis_rna[c(2,3),]

row.names(allbiotypescountspercentsPouchitis_rna)<-allbiotypes[,1]
colnames(allbiotypescountspercentsPouchitis_rna)<-Pouchitis_samples[,1]

row.names(allbiotypescountspercentsPouchitis_rna)<-allbiotypes[,1]
colnames(allbiotypescountspercentsPouchitis_rna)<-Pouchitis_samples[,1]

proteincodingindx = is.element(Gencode_33_Selected_Biotype[,1],allbiotypes[3,1])
lincrnaindx = is.element(Gencode_33_Selected_Biotype[,1],allbiotypes[2,1])
biotypeindx = which(proteincodingindx | lincrnaindx)
length(biotypeindx)
#proteincoding_data<-IGF2BP1ripTPM2[proteincodingindx,]
#install.packages("stringr")              # Install stringr package
library("stringr")

indMT<-str_detect(Gencode_33_Selected_Genename[,1],'^MT-')
indH1<-str_detect(Gencode_33_Selected_Genename[,1], '^H1')
indH2<-str_detect(Gencode_33_Selected_Genename[,1], '^H2')
indH3<-str_detect(Gencode_33_Selected_Genename[,1], '^H3')
indH4<-str_detect(Gencode_33_Selected_Genename[,1], '^H4')
indRPL<-str_detect(Gencode_33_Selected_Genename[,1], '^RPL')
indRPS<-str_detect(Gencode_33_Selected_Genename[,1], '^RPS')

#additgeneid<-Gencode_33_Selected_Geneid[ind,1]
MTgeneid<-Gencode_33_Selected_Geneid[indMT,1]
H1geneid<-Gencode_33_Selected_Geneid[indH1,1]
H2geneid<-Gencode_33_Selected_Geneid[indH2,1]
H3geneid<-Gencode_33_Selected_Geneid[indH3,1]
H4geneid<-Gencode_33_Selected_Geneid[indH4,1]
RPLgeneid<-Gencode_33_Selected_Geneid[indRPL,1]
RPSgeneid<-Gencode_33_Selected_Geneid[indRPS,1]

#additionalgenes<-c(MTgeneid,H1geneid,H2geneid,H3geneid,H4geneid,RPLgeneid,RPSgeneid)
#additional_Geneid_GMask<-unique(additionalgenes)
#length(additional_Geneid_GMask)

additionalgenes<-indMT|indH1|indH2|indH3|indH4|indRPL|indRPS
length(additionalgenes)

nonadditionalgenes<-seq(length(Gencode_33_Selected_Genename[,1]))
nonadditionalgenes<-nonadditionalgenes[!additionalgenes]

mappableindx<-which(Gencode_33_Selected_MappSS[,1]>50)
length(mappableindx)
#protein_Geneid_GMask = Gencode_33_Selected_Geneid[biotypeindx,1]
#length(protein_Geneid_GMask)
#xx<-cbind(Gencode_33_Selected_Geneid,Gencode_33_Selected_MappSS)
#yy<-subset(xx,xx[,2]>50)
#mappableindx<-is.element(Gencode_33_Selected_Geneid[,1], yy[,1])
#mappable_Geneid_GMask = Gencode_33_Selected_Geneid[mappableindx,1]

finalIndexGeneric = intersect(biotypeindx,intersect(nonadditionalgenes,mappableindx));

#--Genes should have at least one count on an average- for 77 samples total count sholud be more than 77.
sm<-apply(Pouchitis_rnaCnt,1,sum)
datindex <- which(sm>77)

finalIndex<-intersect(finalIndexGeneric,datindex)
length(finalIndex)

Pouchitis_rnaCnt_GMask <- Pouchitis_rnaCnt[finalIndex,]

Gencode_33_Selected_Geneid_GMask <- Gencode_33_Selected_Geneid[finalIndex,1]
Gencode_33_Selected_Genename_GMask <- Gencode_33_Selected_Genename[finalIndex,1]
Gencode_33_Selected_MappSS_GMask <- Gencode_33_Selected_MappSS[finalIndex,1]
Gencode_33_Selected_MappUS_GMask <- Gencode_33_Selected_MappUS[finalIndex,1]
Gencode_33_Selected_Biotype_GMask<-Gencode_33_Selected_Biotype[finalIndex,1]

#--- Add Column Names and Add Gene Names- Final Filtered Count Table- Save this Table for Future Reference
colnames(Pouchitis_rnaCnt_GMask)<-Pouchitis_samples[,1]
Pouchitis_rnaCnt_GMask<-cbind(as.data.frame(Gencode_33_Selected_Genename_GMask),Pouchitis_rnaCnt_GMask)

#--- Calculate TPM
#Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#Divide the RPK values by the “per million” scaling factor. This gives you TPM.

r_tpm <- function(dfr,len)
{
  dfr1 <- sweep(dfr,MARGIN=1,(len/10^4),`/`)
  scf <- colSums(dfr1)/(10^6)
  return(sweep(dfr1,2,scf,`/`))
}

#--- Input is Counts and SS gene length
Pouchitis_rnaTPM_GMask <- r_tpm(Pouchitis_rnaCnt_GMask[,2:78], Gencode_33_Selected_MappSS_GMask)
Pouchitis_rnaTPM_GMask<-cbind(as.data.frame(Gencode_33_Selected_Genename_GMask),Pouchitis_rnaTPM_GMask)

#---- Save both Filtered Counts and TPM files to Share with Collaborators
write.table(Pouchitis_rnaCnt_GMask, file="Filtered_Counts/Pouchitis_rnaCnt_GMask.txt", sep = "\t")
write.table(Pouchitis_rnaTPM_GMask, file="Filtered_Counts/Pouchitis_rnaTPM_GMask.txt", sep = "\t")


#----------------------
# - DownStream Analysis
#----------------------
Pouchitis_rnaTPM2<-Pouchitis_rnaTPM_GMask[,-1]
dim(Pouchitis_rnaTPM2)
rn<-nrow(Pouchitis_rnaTPM2)
m<-unique(sample(rn, 1000))
rand_data<-Pouchitis_rnaTPM2[m,]
d <- dist(t(rand_data), method = "euclidean")
for(i in 1:9999){
  m<-unique(sample(rn, 1000))
  rand_data<-Pouchitis_rnaTPM2[m,]
  d <-d+ dist(t(rand_data), method = "euclidean")	
  
}
d1<-d/100000
fit <- hclust(sqrt(d1), method="ward.D2") 
jpeg("Random_1k_Genes_Gmask.jpg", height = 8, width = 16, units = 'in', res = 600)
plot(fit)
dev.off ()

#--- Load metadata
sample_meta <- read.csv(file = "SraRunTable.txt", header = T, row.names = 1, sep = ",") # Check metadata sample names and names in the count file should be in the same order

###########
#---- DESeq
############
library(DESeq2)
library(IHW)
library(ggplot2)

PouchitisrnaCntGMask = Pouchitis_rnaCnt_GMask[,2:78]
dim(PouchitisrnaCntGMask)
PouchitisrnaCntGMask1<-round(PouchitisrnaCntGMask)
Pouchitis_sampleName = as.matrix(Pouchitis_samples)

PouchitisrnaCntGMask_Factor <- DESeqDataSetFromMatrix(PouchitisrnaCntGMask1, colData=sample_meta,design= ~prognosis)
PouchitisrnaCntGMask_Factor <- DESeq(PouchitisrnaCntGMask_Factor)

PouchitisrnaCntGMask_Factor_vsd <- varianceStabilizingTransformation(PouchitisrnaCntGMask_Factor,blind=FALSE)
#write.csv(assay(PouchitisrnaCntGMask_Factor_vsd),file="Filtered_Counts/PouchitisrnaCntGMask_Factor_vsd_DESEq2.csv")

PouchitisrnaCntGMask_Factor_vsd_df <- data.frame(assay(PouchitisrnaCntGMask_Factor_vsd))
#-- Add Gene Names
PouchitisrnaCntGMask_Factor_vsd_df<-cbind(as.data.frame(Gencode_33_Selected_Genename_GMask),PouchitisrnaCntGMask_Factor_vsd_df)
write.table(PouchitisrnaCntGMask_Factor_vsd_df, file="Filtered_Counts/PouchitisrnaCntGMask_Factor_vsd.txt", sep = "\t")

#---------------------------------
#-- Differentially Expressed Genes
#---------------------------------
Genenames<-Gencode_33_Selected_Genename_GMask
PouchitisrnaCntGMask_Factor_diff <- results(PouchitisrnaCntGMask_Factor,contrast=c("prognosis", "Pouchitis", "Healthy"),filterFun=ihw)
rownames(PouchitisrnaCntGMask_Factor_diff) <- Genenames

#################################
#---- Variance Partition Analysis
##################################
#-- On VSD counts
#PouchitisrnaCntGMask_Factor_vsd <- read.csv(file = "Filtered_Counts/PouchitisrnaCntGMask_Factor_vsd.txt", header = T, row.names = 1, sep = "\t")
model <- ~ (1|biopsytime) + (1|Ethnicity) + (1|Diagnosis) + (1|LibraryLayout) + (1|prognosis) + (1|sex)
varPart_model <- fitExtractVarPartModel(PouchitisrnaCntGMask_Factor_vsd[,2:78], model, sample_meta)
dim(varPart_model)
save(varPart_model, file="varPart_model.Rdata")
load ("varPart_model.Rdata")
rownames(varPart_model) = make.names(PouchitisrnaCntGMask_Factor_vsd_df$Gencode_33_Selected_Genename_GMask, unique=TRUE)
varPart_model_df <- data.frame(varPart_model)
write.table(varPart_model_df, file="varPart_model.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2 <- sortCols( varPart_model )
# Violin plot
jpeg("variance_partision.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()
#sort genes based on variance explained by Condition (#Here it is a treatment)
head(varPart_model[order(varPart_model$prognosis, decreasing=TRUE),])

#-----
library(RColorBrewer)
library("pheatmap")

prognosis <- sample_meta$prognosis
prognosis <- data.frame(prognosis)
LibraryLayout <- sample_meta$LibraryLayout
LibraryLayout <- data.frame(LibraryLayout)
biopsytime <- sample_meta$biopsytime
biopsytime <- data.frame(biopsytime)
sample_meta_new <- cbind(prognosis, LibraryLayout, biopsytime)
rownames(sample_meta_new)<- rownames(sample_meta)

sampleDists <- dist(t(assay(PouchitisrnaCntGMask_Factor_vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(sample_meta)
colnames(sampleDistMatrix) <- rownames(sample_meta)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

jpeg("sample_Dist_vsd_heatmap.jpg", height = 7, width = 10, units = 'in', res = 600)
thisheat <- pheatmap(sampleDistMatrix,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     col=colors, annotation_row= sample_meta_new, show_rownames =FALSE,show_colnames = TRUE,legend = TRUE,legend.cex = .05)
dev.off ()

###################
#---- PCA Analysis
###################
Genenames<-Gencode_33_Selected_Genename_GMask

pcarip <- prcomp(t(assay(PouchitisrnaCntGMask_Factor_vsd )))

#percent variance for each component
percentVar <- round(100*pcarip$sdev^2/sum(pcarip$sdev^2))

#retrieve the samples coordinates, and the loadings for each gene in each component
aloadrip <- abs(pcarip$rotation)
ripaloadrelative <- sweep(aloadrip, 2, colSums(aloadrip), "/")
rownames(ripaloadrelative) <- Genenames
#create data frame for plotting
pcaALL <- pcarip$x
pcaR<- data.frame(pcaALL,sample_meta)

#plot with sample labels
library(ggrepel)
jpeg("pca_seqmethod.jpg", height = 7, width = 7, units = 'in', res = 600)
ggplot(pcaR, aes(PC1, PC2, color= prognosis, shape=LibraryLayout)) + geom_point(alpha=0.6,stroke = 3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()  +theme_bw()
dev.off ()




