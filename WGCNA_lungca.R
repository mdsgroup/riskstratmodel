###############
#Preprocessing#
##############

#setwd("./Documents/R_temp/WGCNA_lungca_codes")

# Dataset: Using RMA normalized array data of lung adenocarcinoma 

load("./Rdata/final_normdata.Rdata")
load("./Rdata/final_adc subset.Rdata")
load("./Rdata/final_survival data.RData")

library(WGCNA)
library(ggplot2)
# sample name correction
colnames(GSE50081.rma) <- substr(colnames(GSE50081.rma), 1, 10)
colnames(GSE50081.rma)
colnames(GSE19188.rma) <- substr(colnames(GSE19188.rma), 1, 9)
colnames(GSE19188.rma)
colnames(GSE31546.rma) <- substr(colnames(GSE31546.rma), 1, 9)
colnames(GSE31546.rma)
colnames(GSE31210.rma) <- substr(colnames(GSE31210.rma), 1, 9)
colnames(GSE31210.rma)
colnames(GSE37745.rma) <- substr(colnames(GSE37745.rma), 1, 10)
colnames(GSE37745.rma)
colnames(GSE10245.rma) = substr(colnames(GSE10245.rma), 1, 9)
colnames(GSE10245.rma)
colnames(GSE33532.rma) = substr(colnames(GSE33532.rma), 1, 9)
colnames(GSE33532.rma)
colnames(GSE28571.rma) = substr(colnames(GSE28571.rma), 1, 9)
colnames(GSE28571.rma)
colnames(GSE27716.rma) = substr(colnames(GSE27716.rma), 1, 9)
colnames(GSE27716.rma)
colnames(GSE12667.rma) = substr(colnames(GSE12667.rma), 1, 9)
colnames(GSE12667.rma)

# geneFilter
library(genefilter)
library(hgu133plus2.db)

GSE50081.exp.adc = exprs(featureFilter(GSE50081.rma[,GSE50081.adc]))
GSE19188.exp.adc = exprs(featureFilter(GSE19188.rma[,GSE19188.adc]))
GSE31546.exp.adc = exprs(featureFilter(GSE31546.rma[,GSE31546.adc]))
GSE37745.exp.adc = exprs(featureFilter(GSE37745.rma[,GSE37745.adc]))
GSE10245.exp.adc = exprs(featureFilter(GSE10245.rma[,GSE10245.adc]))
GSE33532.exp.adc = exprs(featureFilter(GSE33532.rma[,GSE33532.adc]))
GSE28571.exp.adc = exprs(featureFilter(GSE28571.rma[,GSE28571.adc]))
GSE27716.exp.adc = exprs(featureFilter(GSE27716.rma[,GSE27716.adc]))
GSE12667.exp.adc = exprs(featureFilter(GSE12667.rma[,GSE12667.adc]))

rownames(GSE50081.exp.adc) = unlist(mget(rownames(GSE50081.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE19188.exp.adc) = unlist(mget(rownames(GSE19188.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE31546.exp.adc) = unlist(mget(rownames(GSE31546.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE37745.exp.adc) = unlist(mget(rownames(GSE37745.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE10245.exp.adc) = unlist(mget(rownames(GSE10245.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE33532.exp.adc) = unlist(mget(rownames(GSE33532.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE28571.exp.adc) = unlist(mget(rownames(GSE28571.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE27716.exp.adc) = unlist(mget(rownames(GSE27716.exp.adc), env = hgu133plus2SYMBOL))
rownames(GSE12667.exp.adc) = unlist(mget(rownames(GSE12667.exp.adc), env = hgu133plus2SYMBOL))

f1 <- pOverA(0.25, log2(100)) 
f2 <- function(x) {IQR(x) > 0.5}
ff <- filterfun(f1)#, f2

# apply filterfunction to each expression matrix
GSE50081.exp.adc <- GSE50081.exp.adc[genefilter(GSE50081.exp.adc, ff),]
GSE19188.exp.adc <- GSE19188.exp.adc[genefilter(GSE19188.exp.adc, ff),]
GSE31546.exp.adc <- GSE31546.exp.adc[genefilter(GSE31546.exp.adc, ff),]
GSE37745.exp.adc <- GSE37745.exp.adc[genefilter(GSE37745.exp.adc, ff),]
GSE10245.exp.adc <- GSE10245.exp.adc[genefilter(GSE10245.exp.adc, ff),]
GSE33532.exp.adc <- GSE33532.exp.adc[genefilter(GSE33532.exp.adc, ff),]
GSE28571.exp.adc <- GSE28571.exp.adc[genefilter(GSE28571.exp.adc, ff),]
GSE27716.exp.adc <- GSE27716.exp.adc[genefilter(GSE27716.exp.adc, ff),]
GSE12667.exp.adc <- GSE12667.exp.adc[genefilter(GSE12667.exp.adc, ff),]

# intersect probeset
int1 <- intersect(rownames(GSE50081.exp.adc), rownames(GSE19188.exp.adc))
int2 <- intersect(int1, rownames(GSE31546.exp.adc))
int4 <- intersect(int2, rownames(GSE37745.exp.adc))
int5 <- intersect(int4, rownames(GSE10245.exp.adc))
int6 <- intersect(int5, rownames(GSE33532.exp.adc))
int7 <- intersect(int6, rownames(GSE28571.exp.adc))
int8 <- intersect(int7, rownames(GSE27716.exp.adc))
int9 <- intersect(int8, rownames(GSE12667.exp.adc))

# select genes
GSE50081.exp.adc <- GSE50081.exp.adc[int9,]
GSE19188.exp.adc <- GSE19188.exp.adc[int9,]
GSE31546.exp.adc <- GSE31546.exp.adc[int9,]
GSE37745.exp.adc <- GSE37745.exp.adc[int9,]
GSE10245.exp.adc <- GSE10245.exp.adc[int9,]
GSE33532.exp.adc <- GSE33532.exp.adc[int9,]
GSE28571.exp.adc <- GSE28571.exp.adc[int9,]
GSE27716.exp.adc <- GSE27716.exp.adc[int9,]
GSE12667.exp.adc <- GSE12667.exp.adc[int9,]

# merge expression matrix (samples with or without survival data)
GSE.exp0 = cbind(GSE50081.exp.adc, GSE19188.exp.adc, GSE31546.exp.adc,
                 GSE37745.exp.adc, GSE10245.exp.adc, GSE33532.exp.adc, GSE28571.exp.adc,
                 GSE27716.exp.adc, GSE12667.exp.adc)
dim(GSE.exp0)

gene_all = rownames(GSE.exp0) 

# Remove batch effect using combat
library(sva)

batch = c(rep(1, ncol(GSE50081.exp.adc)), rep(2, ncol(GSE19188.exp.adc)),
          rep(3, ncol(GSE31546.exp.adc)), 
          rep(5, ncol(GSE37745.exp.adc)), rep(6, ncol(GSE10245.exp.adc)),
          rep(7, ncol(GSE33532.exp.adc)), rep(8, ncol(GSE28571.exp.adc)),
          rep(9, ncol(GSE27716.exp.adc)), rep(10, ncol(GSE12667.exp.adc)))
GSE.combat = ComBat(dat = GSE.exp0, batch = batch, par.prior = T, prior.plots = F)
rownames(GSE.combat) = gene_all # change probe set ids to gene names

# IAC filtering
sizeGrWindow(5,10)
par(mfrow = c(1,2))
IAC <- cor(GSE.combat, use = "p") # cauclating IACs for all pairs of samples
hist(IAC, sub = paste("Mean =", format(mean(IAC[upper.tri(IAC)]), digits = 3)))

meanIAC <- apply(IAC, 2, mean)
sdCorr <- sd(meanIAC)
numbersd <- (meanIAC - mean(meanIAC)) / sdCorr
plot(numbersd)
abline(h = -2, col = "red", lwd = 1)

sdout <- -2
outliers <- colnames(GSE.combat)[numbersd < sdout]
show(outliers)
GSE.filt <- GSE.combat[,numbersd > sdout]
dim(GSE.filt)

IAC2 <- cor(GSE.filt, use="p")
hist(IAC2, sub = paste("Mean =", format(mean(IAC2[upper.tri(IAC2)]), digits = 3)))

meanIAC2 <- apply(IAC2, 2, mean)
sdCorr2 <- sd(meanIAC2)
numbersd2 <- (meanIAC2 - mean(meanIAC2)) / sdCorr2
plot(numbersd2)

#########
# WGCNA #
########

# checking data for excessive missing values and identification of outlier microarray samples
gsg <- goodSamplesGenes(GSE.filt, verbose = 3)ddbs
ghd
gsg$allOK # if TRUe, all genes have passtd the cuts

# transpose Expression dta
GSE.filt<- t(GSE.filt)


## Construction of gene network and identification of modules
#choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# call the network topology analysis function
sft <- pickSoftThreshold(GSE.filt, powerVector = powers, verbose = 5)
# plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
# this line corresponds to using an R^2 cut-off of h = 0.95
abline(h = 0.95, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")


### Two step

softPower <- 6
adjacency <- adjacency(GSE.filt, power = softPower)

## turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

## clustering using TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12, 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = F, hang = 0.04)
minModuleSize <- 30
# module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = F,
                             minClusterSize = minModuleSize)
table(dynamicMods)
# convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


## Merging of modules whose expression profiles are very similar
#calculate eigengenes
MEList <- moduleEigengenes(GSE.filt, colors = dynamicColors)
MEs <- MEList$eigengenes

#Calculate Principal component coefficients
module_gene_pc=list()
module_gene_sqrlatent=list()
module_gene_score=list()

for (i in colnames(MEs))
{
  module_gene = GSE.filt[, dynamicColors==gsub("ME","",i)]
  module_gene_svd = svd(t(scale(module_gene)))
  module_gene_pc[[i]] = module_gene_svd$u[,1] # for pc1 rotation
  module_gene_sqrlatent[[i]] = module_gene_svd$d[1] #for pc1 latent, square root
  module_gene_score[[i]] = scale(module_gene) %*% module_gene_pc[[i]] / module_gene_sqrlatent[[i]]
}


# calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

moduleLabels <- dynamicMods
moduleColors <- dynamicColors


## Visualize the network of eigengenes
sizeGrWindow(5, 7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 2),
                      marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8, xLabelsAngle = 90)


########################
# Survival analysis   #
#######################
library(survival)

# merge surv data
GSE.pheno.v2 = rbind(GSE50081.surv, GSE19188.surv, GSE31546.surv, GSE37745.surv)

index = rownames(GSE.pheno.v2) %in% rownames(GSE.filt)
index = which(index == T)
GSE.pheno.v2 = GSE.pheno.v2[index,] # indexing d/t outlier removal

GSE.pheno.v2$surv = Surv(GSE.pheno.v2$time, GSE.pheno.v2$status == "dead") # build surv object

# subset MEs with survival data
rownames(MEs) = rownames(GSE.filt)
index.ME = rownames(MEs) %in% rownames(GSE.pheno.v2)
index.ME = which(index.ME == T)
MEs = MEs[index.ME,] # MEs with survival data
GSE.filt.surv = GSE.filt[index.ME,] # expression matrix with survival data

# Cox, univariate
res.cox.p<-vector()
res.cox.ci<-vector()

for (i in 1:ncol(MEs))
{
  res.cox1<-coxph(GSE.pheno.v2$surv ~ MEs[,i])
  res.cox1.sum<-summary(res.cox1)
  res.cox.p[i]<-as.numeric(res.cox1.sum$coefficients[,5])
  res.cox.ci[i]<-as.numeric(res.cox1.sum$concondance[1])
}
names(res.cox.p)=colnames(MEs)

# Plotting modules and p-values of univariate cox model.
MEcolors = gsub("ME","",colnames(MEs))
sizeGrWindow(10,10)
plot(-log10(res.cox.p), xlab="Modules", ylab="-log10(p-value)", col=MEcolors,pch=16, cex=1.5, xaxt="n") #cex: size, pch: no of shape  
text(-log10(res.cox.p), MEcolors,cex=0.7, pos=4, col=1) 
axis(1, labels=MEcolors, at=1:ncol(MEs), cex.axis=1, las=2)
cox.index<-which(res.cox.p<0.05) #Significant modules... 
abline(h = -log(0.05, base = 10), lwd = 1, lty = 3, col = "red")


### validataion GSE 31210######
GSE31210.exp.adc = exprs(featureFilter(GSE31210.rma))
rownames(GSE31210.exp.adc) = unlist(mget(rownames(GSE31210.exp.adc), env = hgu133plus2SYMBOL))
GSE31210.exp = GSE31210.exp.adc[gene_all,]

GSE31210.exp= t(GSE31210.exp)

#Extract module eigengenes (modules extracted by multiple GEO dataset) <- from ME coefficient of train set
MEs_val=matrix(0,nrow(GSE31210.exp),length(module_gene_pc))
colnames(MEs_val)=names(module_gene_pc)

for (i in names(module_gene_pc))
{
  module_gene=GSE31210.exp[,dynamicColors==gsub("ME","",i)]
  module_gene_score=scale(module_gene) %*% module_gene_pc[[i]] / module_gene_sqrlatent[[i]]
  MEs_val[,i]=module_gene_score
}


inds= rownames(GSE31210.surv) %in% rownames(GSE31210.exp)
GSE31210.surv = GSE31210.surv[inds,]
GSE31210.surv$surv <- Surv(GSE31210.surv$time, GSE31210.surv$status == "dead")

rownames(MEs_val) = rownames(GSE31210.exp)
index.ME = rownames(MEs_val) %in% rownames(GSE31210.surv)
index.ME = which(index.ME == T)
MEs_val = MEs_val[index.ME,] # MEs with survival data

# Cox, univariate
resval.cox.p<-vector()
resval.cox.ci<-vector()

for (i in 1:ncol(MEs_val))
{
  resval.cox1<-coxph(GSE31210.surv$surv ~ MEs_val[,i])
  resval.cox1.sum<-summary(resval.cox1)
  resval.cox.p[i]<-as.numeric(resval.cox1.sum$coefficients[,5])
  resval.cox.ci[i]<-as.numeric(resval.cox1.sum$concordance[1])
}
names(resval.cox.p)=colnames(MEs)

#Plotting modules and p-values of univariate cox model.

#Original Training & Validation 
tmp=data.frame(MEcolors[cox.index],-log10(resval.cox.p)[cox.index])
colnames(tmp)=c("Modules","p")
p1=ggplot(data=tmp, aes(x=reorder(Modules,p),y=p))+
  geom_bar(stat="identity", fill=tmp$Modules[order(tmp$p)], alpha = 0.8)+
  geom_hline(yintercept=-log10(0.05) ,color="gray20", linetype=3)+
  labs(x="Modules", y= "-log10(p-value)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5), axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_blank())
p1


####GSE30219####
colnames(GSE30219.rma) <- substr(colnames(GSE30219.rma), 1, 9)
GSE30219.exp.adc = exprs(featureFilter(GSE30219.rma))
rownames(GSE30219.exp.adc) = unlist(mget(rownames(GSE30219.exp.adc), env = hgu133plus2SYMBOL))
GSE30219.exp = GSE30219.exp.adc[gene_all,]

GSE30219.exp= t(GSE30219.exp)
#Extract module eigengenes (modules extracted by multiple GEO dataset) <- from ME coefficient of train set
MEs_val2=matrix(0,nrow(GSE30219.exp),length(module_gene_pc))
colnames(MEs_val2)=names(module_gene_pc)

for (i in names(module_gene_pc))
{
  module_gene=GSE30219.exp[,dynamicColors==gsub("ME","",i)]
  module_gene_score=scale(module_gene) %*% module_gene_pc[[i]] / module_gene_sqrlatent[[i]]
  MEs_val2[,i]=module_gene_score
}

inds= rownames(GSE30219.surv) %in% rownames(GSE30219.exp)
GSE30219.surv = GSE30219.surv[inds,]
GSE30219.surv$surv <- Surv(GSE30219.surv$time, GSE30219.surv$status == "dead")

rownames(MEs_val2) = rownames(GSE30219.exp)
index.ME2 = rownames(MEs_val2) %in% rownames(GSE30219.surv)
index.ME2 = which(index.ME2 == T)
MEs_val2 = MEs_val2[index.ME2,] # MEs with survival data

# Cox, univariate
resval2.cox.p<-vector()
resval2.cox.ci<-vector()

for (i in 1:ncol(MEs_val2))
{
  resval2.cox1<-coxph(GSE30219.surv$surv ~ MEs_val2[,i])
  resval2.cox1.sum<-summary(resval2.cox1)
  resval2.cox.p[i]<-as.numeric(resval2.cox1.sum$coefficients[,5])
  resval2.cox.ci[i]<-as.numeric(resval2.cox1.sum$concordance[1])
}
names(resval2.cox.p)=colnames(MEs)

#Plotting modules and p-values of univariate cox model.
tmp=data.frame(MEcolors[cox.index],-log10(resval2.cox.p)[cox.index])
colnames(tmp)=c("Modules","p")
p2=ggplot(data=tmp, aes(x=reorder(Modules,p),y=p))+
  geom_bar(stat="identity", fill=tmp$Modules[order(tmp$p)], alpha = 0.8)+
  geom_hline(yintercept=-log10(0.05) ,color="gray20", linetype=3)+
  labs(x="Modules", y= "-log10(p-value)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5), axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_blank())
p2


##########################################
#Processing for Deep learning modeling   #
##########################################

###Gene module membership
geneModuleMembership = as.data.frame(cor(GSE.filt.surv, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(GSE.filt.surv)))

gMM.module=list()
gMM.order=list()
MM.cox.p=list()
MM.GSE=list()

for (which.ME in names(cox.index))
{
  which.module=gsub("ME","",which.ME)
  gMM.module[[which.module]]=geneModuleMembership[moduleColors==which.module,which.ME]
  names(gMM.module[[which.module]])=gene_all[moduleColors==which.module]
  gMM.order[[which.module]]=order(abs(gMM.module[[which.module]]),decreasing=T)
  
  MM.cox.p[[which.module]]<-vector()
  MM.GSE[[which.module]] <- GSE.filt.surv[,moduleColors==which.module]
  for (i in 1:ncol(MM.GSE[[which.module]]))
  {
    res.cox1<-coxph(GSE.pheno.v2$surv ~ MM.GSE[[which.module]][,i])
    res.cox1.sum<-summary(res.cox1)
    MM.cox.p[[which.module]][i]<-as.numeric(res.cox1.sum$coefficients[,5])
  }
  names(MM.cox.p[[which.module]])=gene_all[moduleColors==which.module]
}

#Plot survival-related modules
sizeGrWindow(5,15)
par(mfrow = c(2,3))

for (which.ME in names(cox.index))
{
  which.module=gsub("ME","",which.ME)
  plot(abs(gMM.module[[which.module]]), -log10(MM.cox.p[[which.module]]), xlab="Gene Module Membership", ylab="-log10(p-value)", cex=0.7, pch=16, col=which.module) #,xlab="Gene Module Membership", ylab="-log10(p-value)"
  text(abs(gMM.module[[which.module]][gMM.order[[which.module]][1:10]]), -log10(MM.cox.p[[which.module]][gMM.order[[which.module]][1:10]]), 
      names(MM.cox.p[[which.module]])[gMM.order[[which.module]][1:10]],cex=0.5, pos=1, col=1, offset=0.1) 
  mtext( paste("r=",toString(round(cor(abs(gMM.module[[which.module]]), -log10(MM.cox.p[[which.module]])),digits=2)),
               "\np=", toString(cor.test(abs(gMM.module[[which.module]]), -log10(MM.cox.p[[which.module]]))$p.value)), 
         cex=0.5)
}

###Export Top genes for significant modules
topno=10 # no. of genes per module

GSE.sig=list()
dir.create("./ModuleGenes", showWarnings = FALSE)

for (which.ME in names(cox.index))
{
  which.module=gsub("ME","",which.ME)
  GSE.sig[[which.module]]=MM.GSE[[which.module]][,gMM.order[[which.module]][1:topno]]
  write.table(GSE.sig[[which.module]], file = paste("./ModuleGenes/GSE",which.module, ".csv",sep=""), sep = ",", quote = TRUE, row.names = FALSE)
}
write.table(GSE.pheno.v2, file="./ModuleGenes/GSEpheno.csv", sep=",", quote=TRUE, row.names=TRUE)
save(GSE.sig, file = "./Rdata/SignificantModules.RData")

#--> To python code..


#Validation set : GSE31210
GSE31210.sig=list()
valindx= rownames(GSE31210.exp) %in% rownames(GSE31210.surv)
GSE31210.exp.val=GSE31210.exp[valindx,]

dir.create("./ModuleGenes_validation", showWarnings = FALSE)
for (which.module in names(GSE.sig))
{
  GSE31210.sig[[which.module]] = GSE31210.exp.val[, colnames(GSE.sig[[which.module]])]
  write.table(GSE31210.sig[[which.module]], file=paste("./ModuleGenes_validation/GSE",which.module,".csv",sep=""),sep = ",", quote = TRUE, row.names = FALSE)
}
write.table(GSE31210.surv, file="./ModuleGenes_validation/GSEpheno.csv", sep=",", quote=TRUE, row.names=TRUE)
#--> TO python testset.


#Validation set : GSE30219

GSE30219.sig=list()
valindx= rownames(GSE30219.exp) %in% rownames(GSE30219.surv)
GSE30219.exp.val=GSE30219.exp[valindx,]

dir.create("./ModuleGenes_validation2", showWarnings = FALSE)
for (which.module in names(GSE.sig))
{
  GSE30219.sig[[which.module]] = GSE30219.exp.val[, colnames(GSE.sig[[which.module]])]
  write.table(GSE30219.sig[[which.module]], file=paste("./ModuleGenes_validation2/GSE",which.module,".csv",sep=""),sep = ",", quote = TRUE, row.names = FALSE)
}
write.table(GSE30219.surv, file="./ModuleGenes_validation2/GSEpheno.csv", sep=",", quote=TRUE, row.names=TRUE)
#--> TO python testset.
