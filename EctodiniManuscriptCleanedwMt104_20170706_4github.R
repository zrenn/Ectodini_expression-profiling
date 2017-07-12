## FILENAME:: EctodiniManuscriptCleanedwMT104_20170706_4github.
## FILENAME:: previously ::: EctodiniManuscriptCleaned_20131206.R
## This code contains all functions needed to 
## 	1) create the genomic mask, 
# 	2) analys gene expression 
# 	3) make the figures for the manuscrt  Renn et al 2017 in GENOME
# Raw data files are available from GEO as dataset GSE97082


################
## SET UP
source("CollapseGenesRenn_etal2017.R")
source("2010_RennJones_functions.R")
source("ExpGenCorrelationFunctions_Renn_etal2017.R")
library (limma)
################

## INITIAL GENE EXPRESSION ANALYSIS  ##############################################
#	only inteersted in 6  species contratsts
####### expression analysis with ALL expression slides
####### NO faint filter (only Flag filter)
####### all pairwise comparisons
####### this is used for the genomic masking

targets <-readTargets("targets_speciesMt104.txt")
#                                         FileName   Cy3   Cy5  Notes
#1                         2008_11_25_ 4.08 27.gpr    Mt    Al   male
#2                         2008_11_25_ 4.08 31.gpr Mt104    Al   male
#3                          2008_12_2_ 4.08 34.gpr    Xf    Xo   male
#4                         2008_11_25_ 4.08 35.gpr    Al Mt104   male
#5                          2008_12_2_ 4.08 36.gpr    Xo    Xf   male
#6                          2008_12_2_ 4.08 37.gpr    Xo    Xf   male
#7                          2008_12_2_ 4.08 38.gpr    Xf    Xo   male
#8                          2008_12_2_ 4.08 39.gpr    Xo    Xf   male
#9                          2008_12_4_ 4.08 40.gpr    Xf    Al   male
#10                         2008_12_2_ 4.08 41.gpr    Xf    Xo   male
#11                         2008_12_2_ 4.08 42.gpr    Xo    Xf   male
#12                         2008_12_2_ 4.08 43.gpr    Xf    Xo   male
#13                         2008_12_2_ 4.08 44.gpr    Xo    Xf   male
#14                         2008_12_2_ 4.08 45.gpr    Xf    Xo   male
#15                         2008_12_4_ 4.08 46.gpr    Al    Xo   male
#16                        2008_11_25_ 4.08 47.gpr    Al    Mt   male
#17                         2008_12_4_ 4.08 48.gpr    Xo    Mt   male
#18                         2008_12_5_ 4.08 49.gpr    Al    Xf   male
#19                        2008_11_26_ 4.08 50.gpr    Al    Mt   male
#20                        2008_11_26_ 4.08 51.gpr    Mt    Al   male
#21                         2008_12_9_ 4.08 52.gpr    Xo    Mt   male
#22                        2008_12_17_ 4.08 53.gpr    Al    Xo   male
#23                        2008_11_26_ 4.08 54.gpr    Al    Mt   male
#24                         2008_12_5_ 4.08 55.gpr    Al    Mt   male
#25                         2008_12_5_ 4.08 56.gpr    Mt    Al   male
#26                        2008_12_09_ 4.08 57.gpr    Mt    Xf   male
#27                        2008_12_10_ 4.08 58.gpr    Xf Mt104   male
#28                        2008_12_10_ 4.08 94.gpr    Mt    Al   male
#29 2007_12_05_ 4.08 63_ Xo116f vs Xf128f fine.gpr    Xo    Xf female
#30 2007_12_06_ 4.08 64_ Xo115f vs Xf123f fine.gpr    Xo    Xf female
#31 2007_12_06_ 4.08 65_ Al112f vs Mt110f fine.gpr    Al    Mt female
#32 2007_12_11_ 4.08 68_ Mt111f vs Al114f fine.gpr    Mt    Al female
#33 2007_12_11_ 4.08 69_ Xf126f vs Al114f fine.gpr    Xf    Al female
#34 2007_12_15_ 4.08 70_ Mt108f vs Xf126f fine.gpr    Mt    Xf female
#35 2007_12_15_ 4.08 71_ Al116f vs Xo113f fine.gpr    Al    Xo female
#36 2007_12_15_ 4.08 72_ Al112f vs Xf127f fine.gpr    Al    Xf female
#37 2007_12_15_ 4.08 73_ Xf127f vs Mt110f fine.gpr    Xf    Mt female
#38 2007_12_15_ 4.08 76_ Xo117f vs Mt111f fine.gpr    Xo    Mt female
#39 2007_12_19_ 4.08 77_ Xf128f vs Xo113f fine.gpr    Xf    Xo female
#40 2007_12_19_ 4.08 78_ Xo113f vs Xf129f fine.gpr    Xo    Xf female
#41   2008_10_17_ 4.08 89_ Al2 vs Xo117f 2b me.gpr    Al    Xo female
#42   2008_10_23_ 4.08 90_ Al2 vs Mt111f 3b me.gpr    Al    Mt female
#43   2008_10_30_ 4.08 91_ Mt110f vs Al2 2b me.gpr    Mt    Al female
#44     2008_11_5_ 4.08 92_ Mt108f vs Al1 2 me.gpr    Mt    Al female
#45 2008_11_07_ 4.08 93_ Al116f vs Mt107f 3 me.gpr    Al    Mt female
#46        2009_1_6_ 4.08 85_ Mt109f vs Al116f.gpr    Mt    Al female
#47        2009_1_6_ 4.08 87_ Mt107f vs Al112f.gpr    Mt    Al female
#48        2009_01_07_ 4.08 24_ Xf126 vs Xo116.gpr    Xf    Xo female
#49          2009_1_7_ 4.08 25_ Xf123 vs Xo117.gpr    Xf    Xo female
#50        2009_1_7_ 4.08 79_ Xo114f vs Xf127f.gpr    Xo    Xf female
#51        2009_1_7_ 4.08 80_ Xf127f vs Xo115f.gpr    Xf    Xo female
#52        2009_1_7_ 4.08 88_ Xo113f vs Mt109f.gpr    Xo    Mt female
#53        2009_01_08_ 4.08 81_ Xo117 vs Xf126.gpr    Xo    Xf female
#54      2009_01_08_ 4.08 82_ Al114f vs Mt108f.gpr    Al    Mt female
#55          2009_1_9_ 4.08 33_ Xf129 vs Xo114.gpr    Xf    Xo female
RG <- read.maimages(targets$FileName,columns=list(R="F635 Median", G="F532 Median", Rb="B635 Median", Gb="B532 Median"),other.columns=c("Flags","B635 SD","B532 SD","F635 % Sat.","F532 % Sat."))
names(RG$other)<-c("Flags", "RbSD", "GbSD","RSat" ,"GSat" )
RG$genes<-readGAL("Fishchip4.03_annotatedGAL_20090730HEM.gal")
RG$printer <- getLayout(RG$genes)
###############
RG$R[RG$other$Flags < 0]<-NA  
RG$G[RG$other$Flags < 0]<-NA
sum(is.na(RG$R)) # 216422
spottypes <- readSpotTypes("SpotTypes.txt")
RG$genes$Status <- controlStatus(spottypes, RG)
batch1<-grep("old", RG$genes$Status);
batch<-grep("old", RG$genes$Status);
MA.qb<-normalizeWithinArraysBatch(RG, method="bPrinttip", batch=batch1)  
# in house function based on normalizeWithinArray (Limma) in loaded script "2010_RennJones_functions.r"
design<- modelMatrix(targets, ref = "Xf") 
contrasts<-makeContrasts(AlXo=Al-Xo, MtXf=Mt,MtXo= Mt-Xo,AlXf=Al,XfXo=-Xo, AlMt=Al - Mt, levels= design)
fit2 <-  Contrasts.fit(MA.qb, design=design, contrasts=contrasts)
fit.bayes <-EBayes(fit2)
# in house function based on ebayes (limma) calculates exact p value rather than estimate 
# in loaded script "2010_RennJones_functions.r"
totable<-cbind(fit.bayes$genes$ID, fit.bayes$coefficients, fit.bayes$p.value)
colnames(totable)<-c("ID","ExpRatio","pvalue")
write.table(totable, file="ResultsEctodiniExpBySpeciesMt104_20100714.txt",sep="\t")
## this file is used for making gneomic mask
save(list=ls(), file="EctodiniExpBySpeciesMt104_20100714.RData")  





allEL<-read.table(file="ResultsEctodiniExpBySpeciesMt104_20100714.txt",sep="\t", header=T)                                

AlXoEL<-allEL[,c(1,2,8)]
MtXfEL<-allEL[,c(1,3,9)]
MtXoEL<-allEL[,c(1,4,10)]
AlXfEL<-allEL[,c(1,5,11)]
XfXoEL<-allEL[,c(1,6,12)]
AlMtEL<-allEL[,c(1,7,13)]
                                     
writeexpa<-read.table(file="ResultsEctodiniGenFlagFilter_20100701.txt",sep="\t", header=T)
AlXo <-writeexpa[,c(4,13,19)]
head(AlXo)
colnames(AlXo)<-c("ID","GenRatio","pvalue")
MtXf <-writeexpa[,c(4,14,20)]
head(MtXf)
colnames(MtXf)<-c("ID","GenRatio","pvalue")
 MtXo <-writeexpa[,c(4,15,21)]
head( MtXo)
colnames( MtXo)<-c("ID","GenRatio","pvalue")
AlXf <-writeexpa[,c(4,16,22)]
head(AlXf)
colnames(AlXf)<-c("ID","GenRatio","pvalue")
XfXo  <-writeexpa[,c(4,17,23)]
head(XfXo )
colnames(XfXo )<-c("ID","GenRatio","pvalue")
AlMt <-writeexpa[,c(4,18,24)]
head(AlMt)
colnames(AlMt)<-c("ID","GenRatio","pvalue")

nAlXo<-na.omit(AlXo)
nMtXf<-na.omit(MtXf)
nMtXo<-na.omit(MtXo)
nAlXf<-na.omit(AlXf)
nXfXo<-na.omit(XfXo)
nAlMt<-na.omit(AlMt)

nAlXoEL<-na.omit(AlXoEL)
nMtXfEL<-na.omit(MtXfEL)
nMtXoEL<-na.omit(MtXoEL)
nAlXfEL<-na.omit(AlXfEL)
nXfXoEL<-na.omit(XfXoEL)
nAlMtEL<-na.omit(AlMtEL)
colnames(nAlXoEL)<-c("ID","ExpRatio","pvalue")
colnames(nMtXfEL)<-c("ID","ExpRatio","pvalue")
colnames(nMtXoEL)<-c("ID","ExpRatio","pvalue")
colnames(nAlXfEL)<-c("ID","ExpRatio","pvalue")
colnames(nXfXoEL)<-c("ID","ExpRatio","pvalue")
colnames(nAlMtEL)<-c("ID","ExpRatio","pvalue")


####################################
####### FUNCTION: "plot.groupsize" 
##  Check to see the effect of including more genes in the SigExpression group
# x= max group size of genes in "significant" group
# stp= the step size 
# ehyb= expression data: must have columns "ID", "ExpRatio", "pvalue"
# ghyb= genomic data: must have columns "ID", "GenRatio", "pvalue"

plot.groupsize<-function(x, stp, ehyb, ghyb, main, ... ){
xs<-x/stp  # the number of steps to plot
mat<-matrix(nrow= xs, ncol= 5)  # matrix that the data will
gh<-merge(ehyb,ghyb, by="ID")
gh2<-gh[order(gh$pvalue.x),]    # order the expression data by pvalue
for(i in 1:xs){
n<-stp*i
m1<-gh2[(1:n),]
m2<-gh2[-(1:n),]
c1<-cor.test(m1$ExpRatio, m1$GenRatio)  
c2<-cor.test(m2$ExpRatio, m2$GenRatio)
mat[i,1]<-n
mat[i,2]<-as.numeric(c1[3])    #pvalue
mat[i,3]<-as.numeric(c1[4])    #estimate
mat[i,4]<-as.numeric(c2[3])    #pvalue  ns
mat[i,5]<-as.numeric(c2[4])    #estimate  ns
}
par(oma=c(0,0,0,2))
plot(mat[,1],mat[,3], pch=16, col="red", ylim=c(-(max(abs(c(mat[,3],mat[,5])))),max(abs(c(mat[,3],mat[,5])))),
  xlab="# sig genes included", ylab="exp/gen correlation")
points(mat[,1],mat[,5], pch=16, col="black", cex=.6) 
par(new=T)  
plot(mat[,1],mat[,2], pch=17,col="purple", yaxt="n", main=main,  ylim=c(0,.1),
  xlab="# sig genes included", ylab="") 
axis(side=4)
mtext("p-value",side=4, line=2.5) 
points(mat[,1],mat[,4], pch=17,col="green", cex=.6) 
abline(h=.05, col="blue") 
legend(x="topright", legend=c("R sig","pvalue sig", "R ns","pvalue ns"),col=c("red","purple","black","green"), pch=c(16, 17), bty="n", pt.cex=c(1,1,.6,.6))
  }
#plot.groupsize(x= 1000, stp=10 , ehyb=oexpAlMt, ghyb=ogenAlMt, main="AlMt sig gene group size")

#################################################################################
####### FUNCTION: "maskplot"
# only apply the genomic hyb mask to those genes that have a pos exp and gen ratio
# or a neg exp and gen ratio
maskplot<-function(x, stp, ehyb, ghyb, main, sig.size, ... ){
ghyb$Ratioab<-abs(ghyb$GenRatio) # to measure magnitude of genomic hyb bias
m1<-merge(ehyb,ghyb, by="ID")
# the following gives a features with gen and exp ratios with opposite bias a "0" gen ratio value
# this avoids "false postive" masking
m1$Ratioab<-ifelse((m1$ExpRatio > 0 &  m1$GenRatio < 0) | (m1$ExpRatio < 0 &  m1$GenRatio > 0), 0, m1$Ratioab) 
m1.o<-m1[order(m1$Ratioab,decreasing = TRUE),] # in decreasing order of gen ratio magnitude
xs<-x/stp # the number of steps to plot
mat<-matrix(nrow= xs, ncol= 5) # will put the correlation results into this new matrix
for(i in 0:xs){
n<-stp*i # number of genes masked (from entire dataset- not known how many are masked from sig group)
m1.del<-m1.o[-(1:n),]
m2<-m1.del[order(m1.del$pvalue.x),]  # order the expression data by pvalue
m.s<-m2[(1:sig.size),]
m.ns<-m2[-(1:sig.size),]
c1<-cor.test(m.s$ExpRatio, m.s$GenRatio)  
c2<-cor.test(m.ns$ExpRatio, m.ns$GenRatio)
mat[i,1]<-n
mat[i,2]<-as.numeric(c1[3])    #pvalue
mat[i,3]<-as.numeric(c1[4])    #estimate
mat[i,4]<-as.numeric(c2[3])    #pvalue  ns
mat[i,5]<-as.numeric(c2[4])    #estimate  ns
}
par(oma=c(0,0,0,2))
plot(mat[,1],mat[,3], pch=16, col="red",  ylim=c(-(max(abs(c(mat[,3],mat[,5])))),max(abs(c(mat[,3],mat[,5])))),
  xlab="# genes masked", ylab="exp/gen correlation")
points(mat[,1],mat[,5], pch=16, col="black", cex=.7) 
par(new=T)  
plot(mat[,1],mat[,2], pch=17,col="purple", yaxt="n", main=main, ylim=c(0,.1),
  xlab="", ylab="", cex=.9) 
axis(side=4)
mtext("p-value",side=4, line=2.5) 
points(mat[,1],mat[,4], pch=17,col="green", cex=.6) 
abline(h=.05, col="blue") 
legend(x="topright", legend=c("sig correlation", "sig p-value", "ns correlation", "ns p-value"),col=c("red","purple","black","green"), pch=c(16, 17), bty="n", pt.cex=c(1,.9,.7,.6))
  }
  
##########################
############# TRIALS
plot.groupsize(x= 5000, stp=10 , ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt sig gene group size")
plot.groupsize(x= 500, stp=4 , ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt sig gene group size") #opt=189
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=2000, stp=10, ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt # genes masked, sig=189", sig.size=189) 
maskplot(x=600, stp=2, ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt # genes masked, sig=189", sig.size=189) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)    #mask 399

plot.groupsize(x= 5000, stp=50 , ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo sig gene group size")
plot.groupsize(x= 500, stp=4 , ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo sig gene group size")
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)  #106
maskplot(x=2000, stp=10, ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo # genes masked, sig=106", sig.size=106) 
maskplot(x=500, stp=2, ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo # genes masked, sig=106", sig.size=106) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)   #mask 384

plot.groupsize(x= 5000, stp=50 , ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf sig gene group size") 
plot.groupsize(x= 200, stp=4 , ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf sig gene group size") #opt=88
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=4, ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf # genes masked, sig=88", sig.size=88) 
maskplot(x=600, stp=2, ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf # genes masked, sig=88", sig.size=88) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)      #mask 547

plot.groupsize(x= 5000, stp=50 , ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo sig gene group size")
plot.groupsize(x= 400, stp=4, ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo sig gene group size") #opt= 119
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=5, ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo # genes masked, sig=119", sig.size=119) 
maskplot(x=600, stp=2, ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo # genes masked, sig=119", sig.size=119) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2) # mask 460

plot.groupsize(x= 5000, stp=10 , ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf sig gene group size")
plot.groupsize(x= 800, stp=5 , ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf sig gene group size") #opt=749
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=5, ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf # genes masked, sig=749", sig.size=749) 
maskplot(x=400, stp=2, ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf # genes masked, sig=749", sig.size=749) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)   #mask 194

plot.groupsize(x= 5000, stp=50 , ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo sig gene group size")
plot.groupsize(x= 1000, stp=4, ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo sig gene group size")# opt=304
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=5, ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo # genes masked, sig=304", sig.size=304)
maskplot(x=200, stp=2, ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo # genes masked, sig=304", sig.size=304)
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2) # mask 36

#########################################
##################  Subsetting genes to be masked into an object for each mask
maskorder<-function(ehyb, ghyb, ... ){
ghyb$Ratioab<-abs(ghyb$GenRatio) # to measure magnitude of genomic hyb bias
m1<-merge(ehyb,ghyb, by="ID")
# the following gives a features with gen and exp ratios with opposite bias a "0" gen ratio value
# this avoids "false postive" masking
m1$Ratioab<-ifelse((m1$ExpRatio > 0 &  m1$GenRatio < 0) | (m1$ExpRatio < 0 &  m1$GenRatio > 0), 0, m1$Ratioab) 
m1.o<-m1[order(m1$Ratioab,decreasing = TRUE),] # in decreasing order of gen ratio magnitude
m1.o}

n2AlXo<- maskorder(nAlXoEL,nAlXo)
n2MtXf<- maskorder(nMtXfEL,nMtXf)        
n2MtXo<- maskorder(nMtXoEL,nMtXo)
n2AlXf<- maskorder(nAlXfEL,nAlXf)
n2XfXo<- maskorder(nXfXoEL,nXfXo)
n2AlMt<- maskorder(nAlMtEL,nAlMt)

mAlXo<-n2AlXo[1:348,] 
mMtXf<-n2MtXf[1:547,] 
mMtXo<-n2MtXo[1:460,] 
mAlXf<-n2AlXf[1:194,] 
mXfXo<-n2XfXo[1:36,] 
mAlMt<-n2AlMt[1:399,] 

m2AlXo<-mAlXo[,1] 
m2MtXf<-mMtXf[,1] 
m2MtXo<-mMtXo[,1] 
m2AlXf<-mAlXf[,1] 
m2XfXo<-mXfXo[,1] 
m2AlMt<-mAlMt[,1]
allmasklist<-c(mAlXo[,1],mMtXf[,1],mMtXo[,1],mAlXf[,1],mXfXo[,1],mAlMt[,1])   #2722 genes
uallmasklist<-unique(allmasklist) #1319 total genes
write.table(m2AlXo, "AlXoMask348_20100714.txt")
write.table(m2MtXf, "MtXfMask547_20100714.txt")
write.table(m2MtXo, "MtXoMask460_20100714.txt")
write.table(m2AlXf, "AlXfMask194_20100714.txt")
write.table(m2XfXo, "XfXoMask36_20100714.txt")
write.table(m2AlMt, "AlMtMask399_20100714.txt")

save(list=ls(), file="EctodiniExpGenCorrAllMt104_20100714.RData")








#################################################################################
####### FUNCTION: "maskplot"
# only apply the genomic hyb mask to those genes that have a pos exp and gen ratio
# or a neg exp and gen ratio
maskplotdiff<-function(x, stp, ehyb, ghyb, main, sig.size, ... ){
ghyb$Ratioab<-abs(ghyb$GenRatio) # to measure magnitude of genomic hyb bias
m1<-merge(ehyb,ghyb, by="ID")
# the following gives a features with gen and exp ratios with opposite bias a "0" gen ratio value
# this avoids "false postive" masking
m1$Ratioab<-ifelse((m1$ExpRatio > 0 &  m1$GenRatio < 0) | (m1$ExpRatio < 0 &  m1$GenRatio > 0), 0, m1$Ratioab) 
m1.o<-m1[order(m1$Ratioab,decreasing = TRUE),] # in decreasing order of gen ratio magnitude
xs<-x/stp # the number of steps to plot
mat<-matrix(nrow= xs, ncol= 5) # will put the correlation results into this new matrix
for(i in 0:xs){
n<-stp*i # number of genes masked (from entire dataset- not known how many are masked from sig group)
m1.del<-m1.o[-(1:n),]
m2<-m1.del[order(m1.del$pvalue.x),]  # order the expression data by pvalue
m.s<-m2[(1:sig.size),]
m.ns<-m2[-(1:sig.size),]
c1<-cor.test(m.s$ExpRatio, m.s$GenRatio)  
c2<-cor.test(m.ns$ExpRatio, m.ns$GenRatio)
mat[i,1]<-n
mat[i,2]<-as.numeric(c1[3])    #pvalue
mat[i,3]<-as.numeric(c1[4])    #estimate
mat[i,4]<-as.numeric(c2[3])    #pvalue  ns
mat[i,5]<-as.numeric(c2[4])    #estimate  ns
}
mat2<-mat
mat2[1,2] <- 0
mat2[1,3] <- 0
mat2[1,4] <- 0
mat2[1,5] <- 0
for (i in 2:length(mat[,1])){
mat2[i,2] <- mat[i,2] - mat[i-1,2]
mat2[i,3] <- mat[i,3] - mat[i-1,3]
mat2[i,4] <- mat[i,4] - mat[i-1,4]
mat2[i,5] <- mat[i,5] - mat[i-1,5]
}
par(oma=c(0,0,0,2))
plot(mat2[,1],mat2[,3], pch=16, col="red", xlab="# genes masked", ylab="exp/gen correlation")
par(new=T)  
plot(mat2[,1],mat2[,2], pch=17,col="purple", yaxt="n", main=main, xlab="", ylab="", cex=.9) 
axis(side=4)
mtext("p-value",side=4, line=2.5) 
abline(h=.05, col="blue") 
legend(x="topright", legend=c("sig correlation", "sig p-value"),col=c("red","purple"), pch=c(16, 17), bty="n", pt.cex=c(1,.9))
  }
  
############
maskplotdiff(x=600, stp=1, ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt # genes masked, sig=189", sig.size=189) 



allEL<-read.table(file="ResultsEctodiniExpBySpeciesMt104_20100714.txt",sep="\t", header=T)                                

AlXoEL<-allEL[,c(1,2,8)]
MtXfEL<-allEL[,c(1,3,9)]
MtXoEL<-allEL[,c(1,4,10)]
AlXfEL<-allEL[,c(1,5,11)]
XfXoEL<-allEL[,c(1,6,12)]
AlMtEL<-allEL[,c(1,7,13)]
                                     
writeexpa<-read.table(file="ResultsEctodiniGenFlagFilter_20100701.txt",sep="\t", header=T)
AlXo <-writeexpa[,c(4,13,19)]
head(AlXo)
colnames(AlXo)<-c("ID","GenRatio","pvalue")
MtXf <-writeexpa[,c(4,14,20)]
head(MtXf)
colnames(MtXf)<-c("ID","GenRatio","pvalue")
 MtXo <-writeexpa[,c(4,15,21)]
head( MtXo)
colnames( MtXo)<-c("ID","GenRatio","pvalue")
AlXf <-writeexpa[,c(4,16,22)]
head(AlXf)
colnames(AlXf)<-c("ID","GenRatio","pvalue")
XfXo  <-writeexpa[,c(4,17,23)]
head(XfXo )
colnames(XfXo )<-c("ID","GenRatio","pvalue")
AlMt <-writeexpa[,c(4,18,24)]
head(AlMt)
colnames(AlMt)<-c("ID","GenRatio","pvalue")

nAlXo<-na.omit(AlXo)
nMtXf<-na.omit(MtXf)
nMtXo<-na.omit(MtXo)
nAlXf<-na.omit(AlXf)
nXfXo<-na.omit(XfXo)
nAlMt<-na.omit(AlMt)

nAlXoEL<-na.omit(AlXoEL)
nMtXfEL<-na.omit(MtXfEL)
nMtXoEL<-na.omit(MtXoEL)
nAlXfEL<-na.omit(AlXfEL)
nXfXoEL<-na.omit(XfXoEL)
nAlMtEL<-na.omit(AlMtEL)
colnames(nAlXoEL)<-c("ID","ExpRatio","pvalue")
colnames(nMtXfEL)<-c("ID","ExpRatio","pvalue")
colnames(nMtXoEL)<-c("ID","ExpRatio","pvalue")
colnames(nAlXfEL)<-c("ID","ExpRatio","pvalue")
colnames(nXfXoEL)<-c("ID","ExpRatio","pvalue")
colnames(nAlMtEL)<-c("ID","ExpRatio","pvalue")


####################################
####### FUNCTION: "plot.groupsize" 
##  Check to see the effect of including more genes in the SigExpression group
# x= max group size of genes in "significant" group
# stp= the step size 
# ehyb= expression data: must have columns "ID", "ExpRatio", "pvalue"
# ghyb= genomic data: must have columns "ID", "GenRatio", "pvalue"

plot.groupsize<-function(x, stp, ehyb, ghyb, main, ... ){
xs<-x/stp  # the number of steps to plot
mat<-matrix(nrow= xs, ncol= 5)  # matrix that the data will
gh<-merge(ehyb,ghyb, by="ID")
gh2<-gh[order(gh$pvalue.x),]    # order the expression data by pvalue
for(i in 1:xs){
n<-stp*i
m1<-gh2[(1:n),]
m2<-gh2[-(1:n),]
c1<-cor.test(m1$ExpRatio, m1$GenRatio)  
c2<-cor.test(m2$ExpRatio, m2$GenRatio)
mat[i,1]<-n
mat[i,2]<-as.numeric(c1[3])    #pvalue
mat[i,3]<-as.numeric(c1[4])    #estimate
mat[i,4]<-as.numeric(c2[3])    #pvalue  ns
mat[i,5]<-as.numeric(c2[4])    #estimate  ns
}
par(oma=c(0,0,0,2))
plot(mat[,1],mat[,3], pch=16, col="red", ylim=c(-(max(abs(c(mat[,3],mat[,5])))),max(abs(c(mat[,3],mat[,5])))),
  xlab="# sig genes included", ylab="exp/gen correlation")
points(mat[,1],mat[,5], pch=16, col="black", cex=.6) 
par(new=T)  
plot(mat[,1],mat[,2], pch=17,col="purple", yaxt="n", main=main,  ylim=c(0,.1),
  xlab="# sig genes included", ylab="") 
axis(side=4)
mtext("p-value",side=4, line=2.5) 
points(mat[,1],mat[,4], pch=17,col="green", cex=.6) 
abline(h=.05, col="blue") 
legend(x="topright", legend=c("R sig","pvalue sig", "R ns","pvalue ns"),col=c("red","purple","black","green"), pch=c(16, 17), bty="n", pt.cex=c(1,1,.6,.6))
  }
#plot.groupsize(x= 1000, stp=10 , ehyb=oexpAlMt, ghyb=ogenAlMt, main="AlMt sig gene group size")

#################################################################################
####### FUNCTION: "maskplot"
# only apply the genomic hyb mask to those genes that have a pos exp and gen ratio
# or a neg exp and gen ratio
maskplot<-function(x, stp, ehyb, ghyb, main, sig.size, ... ){
ghyb$Ratioab<-abs(ghyb$GenRatio) # to measure magnitude of genomic hyb bias
m1<-merge(ehyb,ghyb, by="ID")
# the following gives a features with gen and exp ratios with opposite bias a "0" gen ratio value
# this avoids "false postive" masking
m1$Ratioab<-ifelse((m1$ExpRatio > 0 &  m1$GenRatio < 0) | (m1$ExpRatio < 0 &  m1$GenRatio > 0), 0, m1$Ratioab) 
m1.o<-m1[order(m1$Ratioab,decreasing = TRUE),] # in decreasing order of gen ratio magnitude
xs<-x/stp # the number of steps to plot
mat<-matrix(nrow= xs, ncol= 5) # will put the correlation results into this new matrix
for(i in 0:xs){
n<-stp*i # number of genes masked (from entire dataset- not known how many are masked from sig group)
m1.del<-m1.o[-(1:n),]
m2<-m1.del[order(m1.del$pvalue.x),]  # order the expression data by pvalue
m.s<-m2[(1:sig.size),]
m.ns<-m2[-(1:sig.size),]
c1<-cor.test(m.s$ExpRatio, m.s$GenRatio)  
c2<-cor.test(m.ns$ExpRatio, m.ns$GenRatio)
mat[i,1]<-n
mat[i,2]<-as.numeric(c1[3])    #pvalue
mat[i,3]<-as.numeric(c1[4])    #estimate
mat[i,4]<-as.numeric(c2[3])    #pvalue  ns
mat[i,5]<-as.numeric(c2[4])    #estimate  ns
}
par(oma=c(0,0,0,2))
plot(mat[,1],mat[,3], pch=16, col="red",  ylim=c(-(max(abs(c(mat[,3],mat[,5])))),max(abs(c(mat[,3],mat[,5])))),
  xlab="# genes masked", ylab="exp/gen correlation")
points(mat[,1],mat[,5], pch=16, col="black", cex=.7) 
par(new=T)  
plot(mat[,1],mat[,2], pch=17,col="purple", yaxt="n", main=main, ylim=c(0,.1),
  xlab="", ylab="", cex=.9) 
axis(side=4)
mtext("p-value",side=4, line=2.5) 
points(mat[,1],mat[,4], pch=17,col="green", cex=.6) 
abline(h=.05, col="blue") 
legend(x="topright", legend=c("sig correlation", "sig p-value", "ns correlation", "ns p-value"),col=c("red","purple","black","green"), pch=c(16, 17), bty="n", pt.cex=c(1,.9,.7,.6))
  }
  
##########################
############# TRIALS
plot.groupsize(x= 5000, stp=10 , ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt sig gene group size")
plot.groupsize(x= 500, stp=4 , ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt sig gene group size") #opt=189
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=2000, stp=10, ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt # genes masked, sig=189", sig.size=189) 
maskplot(x=600, stp=2, ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt # genes masked, sig=189", sig.size=189) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)    #mask 399

plot.groupsize(x= 5000, stp=50 , ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo sig gene group size")
plot.groupsize(x= 500, stp=4 , ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo sig gene group size")
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)  #106
maskplot(x=2000, stp=10, ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo # genes masked, sig=106", sig.size=106) 
maskplot(x=500, stp=2, ehyb=nAlXoEL, ghyb=nAlXo, main="AlXo # genes masked, sig=106", sig.size=106) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)   #mask 384

plot.groupsize(x= 5000, stp=50 , ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf sig gene group size") 
plot.groupsize(x= 200, stp=4 , ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf sig gene group size") #opt=88
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=4, ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf # genes masked, sig=88", sig.size=88) 
maskplot(x=600, stp=2, ehyb=nMtXfEL, ghyb=nMtXf, main="MtXf # genes masked, sig=88", sig.size=88) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)      #mask 547

plot.groupsize(x= 5000, stp=50 , ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo sig gene group size")
plot.groupsize(x= 400, stp=4, ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo sig gene group size") #opt= 119
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=5, ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo # genes masked, sig=119", sig.size=119) 
maskplot(x=600, stp=2, ehyb=nMtXoEL, ghyb=nMtXo, main="MtXo # genes masked, sig=119", sig.size=119) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2) # mask 460

plot.groupsize(x= 5000, stp=10 , ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf sig gene group size")
plot.groupsize(x= 800, stp=5 , ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf sig gene group size") #opt=749
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=5, ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf # genes masked, sig=749", sig.size=749) 
maskplot(x=400, stp=2, ehyb=nAlXfEL, ghyb=nAlXf, main="AlXf # genes masked, sig=749", sig.size=749) 
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)   #mask 194

plot.groupsize(x= 5000, stp=50 , ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo sig gene group size")
plot.groupsize(x= 1000, stp=4, ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo sig gene group size")# opt=304
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2)
maskplot(x=1000, stp=5, ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo # genes masked, sig=304", sig.size=304)
maskplot(x=200, stp=2, ehyb=nXfXoEL, ghyb=nXfXo, main="XfXo # genes masked, sig=304", sig.size=304)
text(locator(1),labels=ceiling(as.numeric(locator(1)[1])), cex=2, font=2) # mask 36

#########################################
##################  Subsetting genes to be masked into an object for each mask
maskorder<-function(ehyb, ghyb, ... ){
ghyb$Ratioab<-abs(ghyb$GenRatio) # to measure magnitude of genomic hyb bias
m1<-merge(ehyb,ghyb, by="ID")
# the following gives a features with gen and exp ratios with opposite bias a "0" gen ratio value
# this avoids "false postive" masking
m1$Ratioab<-ifelse((m1$ExpRatio > 0 &  m1$GenRatio < 0) | (m1$ExpRatio < 0 &  m1$GenRatio > 0), 0, m1$Ratioab) 
m1.o<-m1[order(m1$Ratioab,decreasing = TRUE),] # in decreasing order of gen ratio magnitude
m1.o}

n2AlXo<- maskorder(nAlXoEL,nAlXo)
n2MtXf<- maskorder(nMtXfEL,nMtXf)        
n2MtXo<- maskorder(nMtXoEL,nMtXo)
n2AlXf<- maskorder(nAlXfEL,nAlXf)
n2XfXo<- maskorder(nXfXoEL,nXfXo)
n2AlMt<- maskorder(nAlMtEL,nAlMt)

mAlXo<-n2AlXo[1:348,] 
mMtXf<-n2MtXf[1:547,] 
mMtXo<-n2MtXo[1:460,] 
mAlXf<-n2AlXf[1:194,] 
mXfXo<-n2XfXo[1:36,] 
mAlMt<-n2AlMt[1:399,] 

m2AlXo<-mAlXo[,1] 
m2MtXf<-mMtXf[,1] 
m2MtXo<-mMtXo[,1] 
m2AlXf<-mAlXf[,1] 
m2XfXo<-mXfXo[,1] 
m2AlMt<-mAlMt[,1]
allmasklist<-c(mAlXo[,1],mMtXf[,1],mMtXo[,1],mAlXf[,1],mXfXo[,1],mAlMt[,1])   #2722 genes
uallmasklist<-unique(allmasklist) #1319 total genes
write.table(m2AlXo, "AlXoMask348_20100714.txt")
write.table(m2MtXf, "MtXfMask547_20100714.txt")
write.table(m2MtXo, "MtXoMask460_20100714.txt")
write.table(m2AlXf, "AlXfMask194_20100714.txt")
write.table(m2XfXo, "XfXoMask36_20100714.txt")
write.table(m2AlMt, "AlMtMask399_20100714.txt")

save(list=ls(), file="EctodiniExpGenCorrAllMt104_20100714.RData")

#################################################################################
####### FUNCTION: "maskplot"
# only apply the genomic hyb mask to those genes that have a pos exp and gen ratio
# or a neg exp and gen ratio
maskplotdiff<-function(x, stp, ehyb, ghyb, main, sig.size, ... ){
ghyb$Ratioab<-abs(ghyb$GenRatio) # to measure magnitude of genomic hyb bias
m1<-merge(ehyb,ghyb, by="ID")
# the following gives a features with gen and exp ratios with opposite bias a "0" gen ratio value
# this avoids "false postive" masking
m1$Ratioab<-ifelse((m1$ExpRatio > 0 &  m1$GenRatio < 0) | (m1$ExpRatio < 0 &  m1$GenRatio > 0), 0, m1$Ratioab) 
m1.o<-m1[order(m1$Ratioab,decreasing = TRUE),] # in decreasing order of gen ratio magnitude
xs<-x/stp # the number of steps to plot
mat<-matrix(nrow= xs, ncol= 5) # will put the correlation results into this new matrix
for(i in 0:xs){
n<-stp*i # number of genes masked (from entire dataset- not known how many are masked from sig group)
m1.del<-m1.o[-(1:n),]
m2<-m1.del[order(m1.del$pvalue.x),]  # order the expression data by pvalue
m.s<-m2[(1:sig.size),]
m.ns<-m2[-(1:sig.size),]
c1<-cor.test(m.s$ExpRatio, m.s$GenRatio)  
c2<-cor.test(m.ns$ExpRatio, m.ns$GenRatio)
mat[i,1]<-n
mat[i,2]<-as.numeric(c1[3])    #pvalue
mat[i,3]<-as.numeric(c1[4])    #estimate
mat[i,4]<-as.numeric(c2[3])    #pvalue  ns
mat[i,5]<-as.numeric(c2[4])    #estimate  ns
}
mat2<-mat
mat2[1,2] <- 0
mat2[1,3] <- 0
mat2[1,4] <- 0
mat2[1,5] <- 0
for (i in 2:length(mat[,1])){
mat2[i,2] <- mat[i,2] - mat[i-1,2]
mat2[i,3] <- mat[i,3] - mat[i-1,3]
mat2[i,4] <- mat[i,4] - mat[i-1,4]
mat2[i,5] <- mat[i,5] - mat[i-1,5]
}
par(oma=c(0,0,0,2))
plot(mat2[,1],mat2[,3], pch=16, col="red", xlab="# genes masked", ylab="exp/gen correlation")
par(new=T)  
plot(mat2[,1],mat2[,2], pch=17,col="purple", yaxt="n", main=main, xlab="", ylab="", cex=.9) 
axis(side=4)
mtext("p-value",side=4, line=2.5) 
abline(h=.05, col="blue") 
legend(x="topright", legend=c("sig correlation", "sig p-value"),col=c("red","purple"), pch=c(16, 17), bty="n", pt.cex=c(1,.9))
  }
  
############
maskplotdiff(x=600, stp=1, ehyb=nAlMtEL, ghyb=nAlMt, main="AlMt # genes masked, sig=189", sig.size=189) 




#################################################################################
##### repeat Gene Expression Analysis with contrasts of interest applying the gneomic masks from above

#################################################################################


alxo<-read.table("AlXoMask348_20100714.txt")
mtxf<-read.table("MtXfMask547_20100714.txt")
mtxo<-read.table("MtXoMask460_20100714.txt")
alxf<-read.table("AlXfMask194_20100714.txt")
xfxo<-read.table("XfXoMask36_20100714.txt") 
almt<-read.table("AlMtMask399_20100714.txt")
#
full<-c(as.character(alxo[,1]),as.character(mtxf[,1]),as.character(mtxo[,1]),as.character(alxf[,1]),as.character(xfxo[,1]),as.character(almt[,1])) 
fullu<-unique(full) 
length(fullu) #1319 genes 

##############
## READ IN FEMALE
targetsF <-readTargets("targetsF.txt")
targetsF<-targetsF[-27,] # this array is bad
#
RG <- read.maimages(targetsF$FileName,columns=list(R="F635 Median", G="F532 Median", Rb="B635 Median", Gb="B532 Median"),other.columns=c("Flags","B635 SD","B532 SD","F635 % Sat.","F532 % Sat.","Dia."))
names(RG$other)<-c("Flags", "RbSD", "GbSD","RSat" ,"GSat","Dia" )
##### USE 2011 .gal file but note that collapse genes runs on DFCI_TC_2009 so need to change column header.
############### need this to link properly with current annotations 
# 
RG$genes<-readGAL("Fishchip4.03_annotatedGAL_20121012_SCPR.gal", quote ="")
colnames(RG$genes) <- colnames(RG$genes)[c(1:7,9,8,10:21)] # switches 8 and 9 so it uses 2011TC in collapse
RG$printer <- getLayout(RG$genes)
# filter
RG$R[RG$other$Flags < 0]<-NA  
RG$G[RG$other$Flags < 0]<-NA
sum(is.na(RG$R)) #  86225 
RG$R[RG$other$Dia <= 60]<-NA  
RG$G[RG$other$Dia <= 60]<-NA
sum(is.na(RG$R)) # 86251  ### must have been at least one array where the small spots were not flagged
RG$R[RG$genes$Spot == "empty"]<-NA  
RG$G[RG$genes$Spot == "empty"]<-NA
sum(is.na(RG$R)) #86251 #
RG$R[RG$genes$Spot == "control"]<-NA  
RG$G[RG$genes$Spot == "control"]<-NA
sum(is.na(RG$R)) #  91374  ## numbers all match Heather's 
sum(is.na(RG$G))
RG$R[RG$genes$TCorGBorHH == "KickedOut"]<-NA  # these are 2011 TCs
RG$G[RG$genes$TCorGBorHH == "KickedOut"]<-NA  
sum(is.na(RG$R)) # 97112
# filter faint  ### 
RG$R[(RG$R < (RG$Rb + 2*RG$other$RbSD)) & (RG$G < (RG$Gb + 2*RG$other$GbSD))]<-0
RG$G[RG$R == 0]<-NA
RG$R[RG$R == 0]<-NA 
sum(is.na(RG$R)) # 117851
cbind(apply(is.na(RG$R), 2, sum))
# 2007_12_05_ 4.08 63_ Xo116f vs Xf128f fine 3960
# 2007_12_06_ 4.08 64_ Xo115f vs Xf123f fine 5132
# 2007_12_06_ 4.08 65_ Al112f vs Mt110f fine 5379
# 2007_12_11_ 4.08 68_ Mt111f vs Al114f fine 4069
# 2007_12_11_ 4.08 69_ Xf126f vs Al114f fine 3752
# 2007_12_15_ 4.08 70_ Mt108f vs Xf126f fine 3805
# 2007_12_15_ 4.08 71_ Al116f vs Xo113f fine 4481
# 2007_12_15_ 4.08 72_ Al112f vs Xf127f fine 4908
# 2007_12_15_ 4.08 73_ Xf127f vs Mt110f fine 3507
# 2007_12_15_ 4.08 76_ Xo117f vs Mt111f fine 4026
# 2007_12_19_ 4.08 77_ Xf128f vs Xo113f fine 4091
# 2007_12_19_ 4.08 78_ Xo113f vs Xf129f fine 4418
# 2008_10_17_ 4.08 89_ Al2 vs Xo117f 2b me   4446
# 2008_10_23_ 4.08 90_ Al2 vs Mt111f 3b me   4125
# 2008_10_30_ 4.08 91_ Mt110f vs Al2 2b me   3831
# 2008_11_5_ 4.08 92_ Mt108f vs Al1 2 me     4185
# 2008_11_07_ 4.08 93_ Al116f vs Mt107f 3 me 4302
# 2009_1_6_ 4.08 85_ Mt109f vs Al116f        5500
# 2009_1_6_ 4.08 87_ Mt107f vs Al112f        6693
# 2009_01_07_ 4.08 24_ Xf126 vs Xo116        4135
# 2009_1_7_ 4.08 25_ Xf123 vs Xo117          5815
# 2009_1_7_ 4.08 79_ Xo114f vs Xf127f        3693
# 2009_1_7_ 4.08 80_ Xf127f vs Xo115f        3813
# 2009_1_7_ 4.08 88_ Xo113f vs Mt109f        4113
# 2009_01_08_ 4.08 81_ Xo117 vs Xf126        4817
# 2009_01_08_ 4.08 82_ Al114f vs Mt108f      6855
# This is 500 - 1000 more features eliminated for each array
spottypes <- readSpotTypes("SpotTypes.txt")
RG$genes$Status <- controlStatus(spottypes, RG) # Found 4504 old 
batch1<-grep("old", RG$genes$Status);
batch<-grep("old", RG$genes$Status);
## make MA without normalization to mask before normalizing
## played with MA.RG etc and these functions somehow add NAs?
MA.none<-normalizeWithinArrays(RG, method = "none", bc.method = "none", weights = NULL) 
sum(is.na(MA.none$M)) # now  117851 
### same as in RG (before it added one NA and more in other analyses,
###                   must have had to do with faint somehow?)

#########################
############ MASK features
## Heather did masking asfter normalization.  seems like it should be done before
## you wouldn't think it would make that much difference
MA.mask.none <-MA.none
mask<-c(fullu) # made above
for( i in 1:length(mask)){
MA.mask.none$M[which(MA.none$genes$ID==mask[i]),]<-NA}
sum(is.na(MA.mask.none$M)) #  148317

### Now normalize masked data
MA.masked.b<-normalizeWithinArraysBatch(MA.mask.none, method="batch", batch=batch1, bc.method = "minimum") 
# changed from auguest analysis to not be ptip (too few features)
sum(is.na(MA.masked.b$M)) # [1] 148317
############ Collapsing
collapsedF<-CollapseGenes(RG, MA.masked.b)
dim(collapsedF) # 17118   26
#########################
#######  setup ### 
designF<- modelMatrix(targetsF, ref = "Xo116f")
contrastsF<-makeContrasts( 
AlvsMt= ((Al1 +Al112f +Al114f +Al116f +Al2) - (Mt107f +Mt108f +Mt109f +Mt110f +Mt111f))/5,
AlvsXf= ((Al1 +Al112f +Al114f +Al116f +Al2) - (Xf123f +Xf126f +Xf127f +Xf128f +Xf129f))/5,
AlvsXo= ((Al1 +Al112f +Al114f +Al116f +Al2) - (Xo113f +Xo114f +Xo115f +Xo117f))/5,
MtvsXf= ((Mt107f +Mt108f +Mt109f +Mt110f +Mt111f) - (Xf123f +Xf126f +Xf127f +Xf128f +Xf129f))/5 ,
MtvsXo= ((Mt107f +Mt108f +Mt109f +Mt110f +Mt111f) - (Xo113f +Xo114f +Xo115f +Xo117f))/5,
XfvsXo= ((Xf123f +Xf126f +Xf127f +Xf128f +Xf129f) - (Xo113f +Xo114f +Xo115f +Xo117f))/5,
XfXovsAlMt=((Xf123f +Xf126f +Xf127f +Xf128f +Xf129f +Xo113f +Xo114f +Xo115f +Xo117f) -
  (Al1 +Al112f +Al114f +Al116f +Al2 +Mt107f +Mt108f +Mt109f +Mt110f +Mt111f))/10,
XfAlvsXoMt=((Al1 +Al112f +Al114f +Al116f +Al2 + Xf123f +Xf126f +Xf127f +Xf128f +Xf129f) -
  (Xo113f +Xo114f +Xo115f +Xo117f +Mt107f +Mt108f +Mt109f +Mt110f +Mt111f))/10,
AlvsAll=(Al1 +Al112f +Al114f +Al116f +Al2)/5 -
(Xo113f +Xo114f +Xo115f +Xo117f +Mt107f +Mt108f +Mt109f +Mt110f +Mt111f + Xf123f +Xf126f +Xf127f +Xf128f +Xf129f)/15,
XfvsAll=(Xf123f +Xf126f +Xf127f +Xf128f +Xf129f)/5 -
  (Al1 +Al112f +Al114f +Al116f +Al2 + Xo113f +Xo114f +Xo115f +Xo117f +Mt107f +Mt108f +Mt109f +Mt110f +Mt111f)/15,
MtvsAll=(Mt107f +Mt108f +Mt109f +Mt110f +Mt111f)/5 -
  (Al1 +Al112f +Al114f +Al116f +Al2 + Xo113f +Xo114f +Xo115f +Xo117f +Xf123f +Xf126f +Xf127f +Xf128f +Xf129f)/15,
XovsAll=(Xo113f +Xo114f +Xo115f +Xo117f)/5 -
 (Al1 +Al112f +Al114f +Al116f +Al2 +Xf123f +Xf126f +Xf127f +Xf128f +Xf129f+Mt107f +Mt108f +Mt109f +Mt110f +Mt111f)/15, 
XfXoDiffAlMt=((Xf123f +Xf126f +Xf127f +Xf128f +Xf129f) - (Xo113f +Xo114f +Xo115f +Xo117f))/5 -
  ((Al1 +Al112f +Al114f +Al116f +Al2) - (Mt107f +Mt108f +Mt109f +Mt110f +Mt111f))/5,  levels=designF)

#########################
#### RUN FIT  for exact P-values 
fit2F <-  Contrasts.fit(collapsedF, design=designF, contrasts=contrastsF)
fit2Find <-  Contrasts.fit(collapsedF, design=designF)
fit.bayesF <-EBayes(fit2F)
fit.bayesFind <-EBayes(fit2Find)
resultsfdr05F<-decideTests(fit.bayesF, method = "separate", adjust.method="fdr", p.value=0.05)
summary (resultsfdr05F)
   # AlvsMt AlvsXf AlvsXo MtvsXf MtvsXo XfvsXo XfXovsAlMt XfAlvsXoMt AlvsAll XfvsAll MtvsAll XovsAll XfXoDiffAlMt
# -1    181      2      2    488    450      6         91         36       6       6     732      61          180
# 0    9284   9877   9880   9228   9117   9878       9611       9698    9861    9872    8721    9769         9603
# 1     423      9      6    172    321      4        186        154      21      10     435      58          105

181+9284+423 # 9888 = the number of gene analyzed 
results005F<-decideTests(fit.bayesF, method = "separate", adjust.method="none", p.value=0.005)
summary (results005F)
   # AlvsMt AlvsXf AlvsXo MtvsXf MtvsXo XfvsXo XfXovsAlMt XfAlvsXoMt AlvsAll XfvsAll MtvsAll XovsAll XfXoDiffAlMt
# -1    218     85     89    563    484     52        193         90      48      86     705     203          269
# 0    9173   9726   9684   9115   9055   9745       9375       9488    9747    9621    8773    9491         9415
# 1     497     77    115    210    349     91        320        310      93     181     410     194          204

which(results005F[,6]==1 & results005F[,4]==-1 & results005F[,2]==-1)  # 3 Xf up in Female

results05F<-decideTests(fit.bayesF, method = "separate", adjust.method="none", p.value=0.05)
summary (results05F)
   # AlvsMt AlvsXf AlvsXo MtvsXf MtvsXo XfvsXo XfXovsAlMt XfAlvsXoMt AlvsAll XfvsAll MtvsAll XovsAll XfXoDiffAlMt
# -1    683    434    367   1211   1035    230        623        455     272     369    1326     642          662
# 0    8152   9125   9057   7947   7964   9174       8414       8482    9250    8696    7566    8628         8577
# 1    1053    329    464    730    889    484        851        951     366     823     996     618          649

resultsfdr1F<-decideTests(fit.bayesF, method = "separate", adjust.method="fdr", p.value=0.1)
summary (resultsfdr1F)
   # AlvsMt AlvsXf AlvsXo MtvsXf MtvsXo XfvsXo XfXovsAlMt XfAlvsXoMt AlvsAll XfvsAll MtvsAll XovsAll XfXoDiffAlMt
# -1    283     10     30    732    608     13        197         64       7      21     950     168          265
# 0    9015   9861   9822   8839   8794   9864       9362       9560    9849    9832    8292    9564         9428
# 1     590     17     36    317    486     11        329        264      32      35     646     156          195
results005Find<-decideTests(fit.bayesFind, method = "separate", adjust.method="none", p.value=0.005)
summary(results005Find) # reference is Xo116f
     # Al1 Al112f Al114f Al116f   Al2 Mt107f Mt108f Mt109f Mt110f Mt111f Xf123f Xf126f Xf127f Xf128f Xf129f Xo113f Xo114f Xo115f Xo117f
# -1    11     87     28     69    12     89    134     60     71     12     23     32     28     36     17     21      9     10     14
# 0  12250  12441  12959  12959 13117  12265  12585  12629  13036  13011  12390  12844  13063  12800  11993  12175  12426  12316  12124
# 1     34    103     34     32    35     46    106     21     20     19     35     57     17     26     27     14      8      6     28

############
############
### MALE ###
############
############

targetsM <-readTargets("targetsM.txt")
targetsM                 # FileName   Cy3   Cy5 Notes
# 1  2008_11_25_ 4.08 27.gpr Mt105 Al106  male
# 2  2008_11_25_ 4.08 31.gpr Mt104 Al110  male
# 3   2008_12_2_ 4.08 34.gpr Xf121 Xo101  male
# 4  2008_11_25_ 4.08 35.gpr Al109 Mt104  male
# 5   2008_12_2_ 4.08 36.gpr Xo101 Xf117  male
# 6   2008_12_2_ 4.08 37.gpr Xo111 Xf118  male
# 7   2008_12_2_ 4.08 38.gpr Xf120 Xo111  male
# 8   2008_12_2_ 4.08 39.gpr Xo110 Xf120  male
# 9   2008_12_4_ 4.08 40.gpr Xf120 Al106  male
# 10  2008_12_2_ 4.08 41.gpr Xf118 Xo109  male
# 11  2008_12_2_ 4.08 42.gpr Xo109 Xf121  male
# 12  2008_12_2_ 4.08 43.gpr Xf117 Xo107  male
# 13  2008_12_2_ 4.08 44.gpr Xo107 Xf122  male
# 14  2008_12_2_ 4.08 45.gpr Xf122 Xo110  male
# 15  2008_12_4_ 4.08 46.gpr Al110 Xo110  male
# 16 2008_11_25_ 4.08 47.gpr Al110 Mt105  male
# 17  2008_12_4_ 4.08 48.gpr Xo110 Mt105  male
# 18  2008_12_5_ 4.08 49.gpr Al109 Xf117  male
# 19 2008_11_26_ 4.08 50.gpr Al103 Mt103  male
# 20 2008_11_26_ 4.08 51.gpr Mt103 Al107  male
# 21  2008_12_9_ 4.08 52.gpr Xo109 Mt103  male
# 22 2008_12_17_ 4.08 53.gpr Al107 Xo109  male
# 23 2008_11_26_ 4.08 54.gpr Al106 Mt102  male
# 24  2008_12_5_ 4.08 55.gpr Al107 Mt101  male
# 25  2008_12_5_ 4.08 56.gpr Mt101 Al109  male
# 26 2008_12_09_ 4.08 57.gpr Mt102 Xf120  male
# 27 2008_12_10_ 4.08 58.gpr Xf117 Mt104  male
# 28 2008_12_10_ 4.08 94.gpr Mt102 Al103  male

RGM <- read.maimages(targetsM$FileName,columns=list(R="F635 Median", G="F532 Median", Rb="B635 Median", Gb="B532 Median"),other.columns=c("Flags","B635 SD","B532 SD","F635 % Sat.","F532 % Sat."))
names(RGM$other)<-c("Flags", "RbSD", "GbSD","RSat" ,"GSat" )
####### 
RGM$genes<-readGAL("Fishchip4.03_annotatedGAL_20121012_SCPR.gal", quote ="")
colnames(RGM$genes) <- colnames(RGM$genes)[c(1:7,9,8,10:21)] # switches 8 and 9
RGM$printer <- getLayout(RGM$genes)
#######
RGM$R[RGM$other$Flags < 0]<-NA  
RGM$G[RGM$other$Flags < 0]<-NA
sum(is.na(RGM$R)) # 119884
RGM$R[RGM$other$Dia <= 60]<-NA  
RGM$G[RGM$other$Dia <= 60]<-NA
sum(is.na(RGM$R)) # 119884
RGM$R[RGM$genes$Spot == "empty"]<-NA  
RGM$G[RGM$genes$Spot == "empty"]<-NA
sum(is.na(RGM$R)) # 119884
RGM$R[RGM$genes$Spot == "control"]<-NA  
RGM$G[RGM$genes$Spot == "control"]<-NA
sum(is.na(RGM$R)) #  124965
RGM$R[RGM$genes$TCorGBorHH == "KickedOut"]<-NA  
RGM$G[RGM$genes$TCorGBorHH == "KickedOut"]<-NA  
sum(is.na(RGM$R)) # [1] 130798
# filter faint  ### 
RGM$R[(RGM$R < (RGM$Rb + 2*RGM$other$RbSD)) & (RGM$G < (RGM$Gb + 2*RGM$other$GbSD))]<-0
RGM$G[RGM$R == 0]<-NA
RGM$R[RGM$R == 0]<-NA 
sum(is.na(RGM$R)) # [1] 117881 # shockingly similar to female
cbind(apply(is.na(RGM$R), 2, sum))
                     # [,1]
# 2008_11_25_ 4.08 27 15009
# 2008_11_25_ 4.08 31 13566
# 2008_12_2_ 4.08 34   4882
# 2008_11_25_ 4.08 35 12625
# 2008_12_2_ 4.08 36   5733
# 2008_12_2_ 4.08 37   4149
# 2008_12_2_ 4.08 38   4477
# 2008_12_2_ 4.08 39   4611
# 2008_12_4_ 4.08 40   4404
# 2008_12_2_ 4.08 41   3977
# 2008_12_2_ 4.08 42   5049
# 2008_12_2_ 4.08 43   3891
# 2008_12_2_ 4.08 44   4664
# 2008_12_2_ 4.08 45   7078
# 2008_12_4_ 4.08 46   4333
# 2008_11_25_ 4.08 47  4405
# 2008_12_4_ 4.08 48   4900
# 2008_12_5_ 4.08 49   3708
# 2008_11_26_ 4.08 50  5100
# 2008_11_26_ 4.08 51  6396
# 2008_12_9_ 4.08 52   4147
# 2008_12_17_ 4.08 53  4386
# 2008_11_26_ 4.08 54  4303
# 2008_12_5_ 4.08 55   7376
# 2008_12_5_ 4.08 56   4936
# 2008_12_09_ 4.08 57  4828
# 2008_12_10_ 4.08 58 15180
# 2008_12_10_ 4.08 94  9768 #### holy shit, some are so faint.... that the whole array is nearly eliminated
cbind(apply(!is.na(RGM$R), 2, sum))
                     # [,1]
# 2008_11_25_ 4.08 27  4191
# 2008_11_25_ 4.08 31  5634
# 2008_12_2_ 4.08 34  14318
# 2008_11_25_ 4.08 35  6575
# 2008_12_2_ 4.08 36  13467
# 2008_12_2_ 4.08 37  15051
# 2008_12_2_ 4.08 38  14723
# 2008_12_2_ 4.08 39  14589
# 2008_12_4_ 4.08 40  14796
# 2008_12_2_ 4.08 41  15223
# 2008_12_2_ 4.08 42  14151
# 2008_12_2_ 4.08 43  15309
# 2008_12_2_ 4.08 44  14536
# 2008_12_2_ 4.08 45  12122
# 2008_12_4_ 4.08 46  14867
# 2008_11_25_ 4.08 47 14795
# 2008_12_4_ 4.08 48  14300
# 2008_12_5_ 4.08 49  15492
# 2008_11_26_ 4.08 50 14100
# 2008_11_26_ 4.08 51 12804
# 2008_12_9_ 4.08 52  15053
# 2008_12_17_ 4.08 53 14814
# 2008_11_26_ 4.08 54 14897
# 2008_12_5_ 4.08 55  11824
# 2008_12_5_ 4.08 56  14264
# 2008_12_09_ 4.08 57 14372
# 2008_12_10_ 4.08 58  4020
# 2008_12_10_ 4.08 94  9432  # these are good spots after masking

##### get infor for batch normalization
RGM$genes$Status <- controlStatus(spottypes, RGM)
batch1<-grep("old", RGM$genes$Status);  # can't remember why I need both??
batch<-grep("old", RGM$genes$Status);

#########################
## make MA without normalization to mask before normalizing
MAM.none<-normalizeWithinArrays(RGM, method = "none", bc.method = "none", weights = NULL) # make MA
dim(MAM.none) # [1] 19200    28
#########################
########## Masking out genes
sum(is.na(MAM.none$M)) # [1]  177881 # same as before as I'd expect
MAM.mask.none <-MAM.none
mask<-c(fullu) # made above
for( i in 1:length(mask)){
MAM.mask.none$M[which(MAM.none$genes$ID==mask[i]),]<-NA}
sum(is.na(MAM.mask.none$M))# now  206932

#########################
## Now normalize masked data
MAM.masked.b<-normalizeWithinArraysBatch(MAM.mask.none, method="batch", batch=batch1)
sum(is.na(MAM.masked.b$M)) # [1]  206932


####20170312 add a step here to write a files for GEO submission)
## GEO platform GPL6416 does not have empty features. use $genes$ID but will have to eliminate all "empty"
writeGEO.table <- cbind(MA.mask.none$genes$ID, MA.masked.b$M, MAM.masked.b$M)
colnames(writeGEO.table)[1] <- "ID_REF"
dim(writeGEO.table) # [1] 19200    55
# will have to add .gpr in excel
writeGEO.tableE <- writeGEO.table[-c(which(writeGEO_F.table[,1] == "empty")), ]
dim(writeGEO.tableE) # [1] 17722    55
write.table(writeGEO.tableE,file="EctodiniF_RNA_GE0_normMvalues_20170312.tab",  row.names = FALSE, sep="\t")


#########################
################ Collapsing
collapsedM<-CollapseGenes(RGM, MAM.masked.b)
dim(collapsedM)  # [1] 17118    28

#########################
################ set up
designM<- modelMatrix(targetsM, ref = "Al103") 
contrastsM<-makeContrasts( 
AlvsMt=(Al106 +Al107 +Al109 +Al110)/5 - (Mt104 +Mt101 +Mt102 +Mt103 +Mt105)/5,
AlvsXf=(Al106 +Al107 +Al109 +Al110)/5 - (Xf117 +Xf118 +Xf120 +Xf121 +Xf122)/5,
AlvsXo=(Al106 +Al107 +Al109 +Al110)/5 - (Xo101 +Xo107 +Xo109 +Xo110 +Xo111 )/5,
MtvsXf=(Mt104 +Mt101 +Mt102 +Mt103+Mt105)/5 - (Xf117 +Xf118 +Xf120 +Xf121 +Xf122)/5,
MtvsXo=(Mt104 +Mt101 +Mt102 +Mt103+Mt105)/5 - (Xo101 +Xo107 +Xo109 +Xo110 +Xo111 )/5,
XfvsXo=((Xf117 +Xf118 +Xf120 +Xf121 +Xf122) - (Xo101 +Xo107 +Xo109 +Xo110 +Xo111 ))/5,
XfXovsAlMt=((Xf117 +Xf118 +Xf120 +Xf121 +Xf122 +Xo101 +Xo107 +Xo109 +Xo110 +Xo111 ) -
  (Al106 +Al107 +Al109 +Al110 +Mt104 +Mt101 +Mt102 +Mt103+Mt105))/10,
XfAlvsXoMt=((Xf117 +Xf118 +Xf120 +Xf121 +Xf122 +Al106 +Al107 +Al109 +Al110) -
  (Xo101 +Xo107 +Xo109 +Xo110 +Xo111 +Mt104 +Mt101 +Mt102 +Mt103+Mt105))/10,
AlvsAll= (Al106 +Al107 +Al109 +Al110)/5 - (Xo101 +Xo107 +Xo109 +Xo110 +Xo111 +Mt104 +Mt101 +Mt102 +Mt103+Mt105 + Xf117 +Xf118 +Xf120 +Xf121 +Xf122)/15,
XfvsAll= (Xf117 +Xf118 +Xf120 +Xf121 +Xf122)/5 - (Mt104 +Mt101 +Mt102 +Mt103+Mt105 + Xo101 +Xo107 +Xo109 +Xo110 +Xo111 + Al106 +Al107 +Al109 +Al110)/15,
MtvsAll= (Mt104 +Mt101 +Mt102 +Mt103+Mt105)/5 - (Xf117 +Xf118 +Xf120 +Xf121 +Xf122+Xo101 +Xo107 +Xo109 +Xo110 +Xo111 + + Al106 +Al107 +Al109 +Al110)/15,
XovsAll= (Xo101 +Xo107 +Xo109 +Xo110 +Xo111)/5 - (Mt104 +Mt101 +Mt102 +Mt103 +Mt105+Xf117 +Xf118 +Xf120 +Xf121 + Al106 +Al107 +Al109 +Al110)/15,
XfXoDiffAlMt= ((Xf117 +Xf118 +Xf120 +Xf121 +Xf122) - (Xo101 +Xo107 +Xo109 +Xo110 +Xo111 ))/5 - ((Al106 +Al107 +Al109 +Al110)/5 - (Mt104 +Mt101 +Mt102 +Mt103 +Mt105)/5),
  levels=designM)

#########################
############### RUN the Fit  
fit2M <-  Contrasts.fit(collapsedM, design=designM, contrasts=contrastsM)
fit.bayesM <-EBayes(fit2M)

fit2Mind <-  Contrasts.fit(collapsedM, design=designM)
fit.bayesMind <-EBayes(fit2Mind)

resultsfdr05M<-decideTests(fit.bayesM, method = "separate", adjust.method="fdr", p.value=0.05)
summary (resultsfdr05M)
   # AlvsMt AlvsXf AlvsXo MtvsXf MtvsXo XfvsXo XfXovsAlMt XfAlvsXoMt AlvsAll XfvsAll MtvsAll XovsAll XfXoDiffAlMt
# -1     32      0      0     91      4     18          1         20       0      10     105       0           14
# 0    5592   5696   5696   5556   5683   5665       5695       5639    5696    5661    5543    5696         5676
# 1      72      0      0     49      9     13          0         37       0      25      48       0            6

32+5592+72 # only 5695 genes are analyzed.... this is even more horrific
### without filtering faint 265+10402+461 # 11128 suvived it was 6907 without Mt104
#########################  ## how are there none for AlvsXo ### that's what Heather had too
results005M<-decideTests(fit.bayesM, method = "separate", adjust.method="none", p.value=0.005)
summary (results005M)
   # AlvsMt AlvsXf AlvsXo MtvsXf MtvsXo XfvsXo XfXovsAlMt XfAlvsXoMt AlvsAll XfvsAll MtvsAll XovsAll XfXoDiffAlMt
# -1     67     62     22    172    103     47         71         48      35      99     173      19           88
# 0    5473   5561   5656   5426   5526   5594       5529       5511    5619    5456    5432    5654         5544
# 1     156     73     18     98     67     55         96        137      42     141      91      23           64

results05M<-decideTests(fit.bayesM, method = "separate", adjust.method="none", p.value=0.05)
summary (results05M)   
   # AlvsMt AlvsXf AlvsXo MtvsXf MtvsXo XfvsXo XfXovsAlMt XfAlvsXoMt AlvsAll XfvsAll MtvsAll XovsAll XfXoDiffAlMt
# -1    241    377    143    614    330    268        378        189     199     526     524     212          441
# 0    4958   4837   5320   4692   5120   5028       4949       4962    5153    4516    4895    5334         5000
# 1     497    482    233    390    246    400        369        545     344     654     277     150          255


results005Mind<-decideTests(fit.bayesMind, method = "separate", adjust.method="none", p.value=0.005)
summary (results005Mind)
   # Al106 Al107 Al109 Al110 Mt101 Mt102 Mt103 Mt104 Mt105 Xf117 Xf118 Xf120 Xf121 Xf122 Xo101 Xo107 Xo109 Xo110 Xo111
# -1    21    34    46    10    38    21    13    75    11    43     9    34    38    22    27    26     9    14     7
# 0  12727 12485 12762 12695 11792 12702 12702  6426 12339 12608 12707 12607 12247 12299 11751 11617 12442 11800 11346
# 1      8     9    71     6    22    11    16     8     5    57     4    12   126    84    62    30    11     3     0


head(fit.bayesMind$coeff) ## there are NAs in some and not others so clustering would be driven by missingness
# should do PCA on individuals for only those that pass in both male and female (will there still be NAs?)


#####################################
#################### Redo the genomic analysis in collapsing the features for a TC to make figures
targetsG <-readTargets("targetsgenomic.txt")
#                                  FileName Cy3 Cy5 Notes
#1  p406s47HEM20090410_Xo3_Xf5_20090629.gpr  Xo  Xf    48
#2  p406s48HEM20090415_Xf3_Al5_20090629.gpr  Xf  Al    48
#3  p406s50HEM20090415_Xf3_Mt5_20090629.gpr  Xf  Mt    48
#4  p406s51HEM20090416_Al3_Mt5_20090629.gpr  Al  Mt    48
#5  p406s52HEM20090416_Al3_Xo5_20090629.gpr  Al  Xo    48
#6 p406s53HEM20090416_Mt3_Xo5b_20090629.gpr  Mt  Xo    48
RGG <- read.maimages(targetsG$FileName,columns=list(R="F635 Median", G="F532 Median", Rb="B635 Median", Gb="B532 Median"),other.columns=c("Flags","B635 SD","B532 SD","F635 % Sat.","F532 % Sat."))
names(RG$other)<-c("Flags", "RbSD", "GbSD","RSat" ,"GSat" )
#######  this new gal with Sec
# RGG$genes<-readGAL("Fishchip4.03_annotatedGAL_20090730HEM.gal")
##### USE 2011 .gal file but note that collapse genes runs on DFCI_TC_2009 so need to change column header.
############### need this to link properly with current annotations 
# 
RGG$genes<-readGAL("Fishchip4.03_annotatedGAL_20121012_SCPR.gal", quote ="")
colnames(RGG$genes) <- colnames(RGG$genes)[c(1:7,9,8,10:21)] # switches 8 and 9 so it uses 2011TC in collapse
RGG$printer <- getLayout(RGG$genes)
# filter
RGG$R[RGG$other$Flags < 0]<-NA  
RGG$G[RGG$other$Flags < 0]<-NA
sum(is.na(RGG$R)) #  29118 
RGG$R[RGG$genes$Spot == "empty"]<-NA  
RGG$G[RGG$genes$Spot == "empty"]<-NA
sum(is.na(RGG$R)) #29118 #
RGG$R[RGG$genes$Spot == "control"]<-NA  
RGG$G[RGG$genes$Spot == "control"]<-NA
sum(is.na(RGG$R)) #  29835  ## 
RGG$R[RGG$genes$TCoRGGBorHH == "KickedOut"]<-NA  # these are 2011 TCs
RGG$G[RGG$genes$TCoRGGBorHH == "KickedOut"]<-NA  
sum(is.na(RGG$R)) # 29835 ## there should be a couple hundred
cbind(apply(is.na(RGG$R), 2, sum))                                     
# [,1]
# p406s47HEM20090410_Xo3_Xf5_20090629  9450
# p406s48HEM20090415_Xf3_Al5_20090629  3109
# p406s50HEM20090415_Xf3_Mt5_20090629  3986
# p406s51HEM20090416_Al3_Mt5_20090629  4152
# p406s52HEM20090416_Al3_Xo5_20090629  4338
# p406s53HEM20090416_Mt3_Xo5b_20090629 4350
#########################
## make MA without normalization to mask before normalizing
MAG.none<-normalizeWithinArrays(RGG, method = "none", bc.method = "none", weights = NULL) # make MA
dim(MAG.none) # [1] 19200    6
#########################
########## Masking out genes
sum(is.na(MAG.none$M)) # [1]  177881 # same as before as I'd expect
MAG.mask.none <-MAG.none
mask<-c(fullu) # made above
for( i in 1:length(mask)){
MAG.mask.none$M[which(MAG.none$genes$ID==mask[i]),]<-NA}
sum(is.na(MAG.mask.none$M))# now  36424
#########################
## Now normalize masked data
MAG.masked.b<-normalizeWithinArraysBatch(MAG.mask.none, method="batch", batch=batch1)
sum(is.na(MAG.masked.b$M)) # [1]  36424
#########################
################ Collapsing
collapsedG<-CollapseGenes(RGG, MAG.masked.b)
dim(collapsedG)  # [1] 17118    6
#########################
################ set up
designG<- modelMatrix(targetsG, ref = "Xf")
contrastsG<-makeContrasts(AlMt=Al-Mt, AlXf=Al,AlXo=Al-Xo, MtXf=Mt, MtXo= Mt-Xo, XfXo=-Xo,  levels= designG)
fitG <-  Contrasts.fit(collapsedG, design=designG,contrasts=contrastsG)
fitG.bayes <-EBayes(fitG)
resultsGn005<-decideTests(fitG.bayes, method = "separate", adjust.method="none", p.value=0.005)
summary(resultsGn005) # But these don't really matter
    # AlMt  AlXf  AlXo  MtXf  MtXo  XfXo
# -1    24    35    16    78    39    28
# 0  12146 12249 12250 12205 12260 12228
# 1    145    31    49    32    16    59
colnames(fitG.bayes$coefficients) # "AlMt" "AlXf" "AlXo" "MtXf" "MtXo" "XfXo"
writeG<-cbind(fitG.bayes$genes, fitG.bayes$coefficients,  fitG.bayes$p.value, fitG.bayes$lods)
colnames(writeG)<-c(colnames(fitG.bayes$genes),c(paste("coef.g",c("AlvsMt", "AlvsXf" ,"AlvsXo", "MtvsXf", "MtvsXo", "XfvsXo"), sep=".")),c(paste("P.g",c("AlvsMt", "AlvsXf" ,"AlvsXo", "MtvsXf", "MtvsXo", "XfvsXo"), sep=".")),c(paste("lods.g",c("AlvsMt", "AlvsXf" ,"AlvsXo", "MtvsXf", "MtvsXo", "XfvsXo"), sep="."))) 
head(writeG)
write.table(writeG, file="21031010_ResultsEctodiniGenFlagFilterLODScollapsed.txt",sep="\t")

## counting  NEED TO FIRST LIMIT TO THOSE THAT SURVIVE BOTH LOOPS
#############################################################################
# MALE
## make a table to merge and keep results around

writeM.table <- cbind(fit.bayesM$genes[, c(1:4,8,9,13:15,21), ], fit.bayesM$lods, fit.bayesM$coeff, fit.bayesM$p.value, results005M)
#writeM.table <- cbind(fit.bayesM$genes[, c(1:4,8,9,13:15,21), ], fit.bayesM$lods, fit.bayesM$coeff, fit.bayesM$p.value, resultsfdr05M)
#writeM.table <- cbind(fit.bayesM$genes[, c(1:4,8,9,13:15,21), ], fit.bayesM$lods, fit.bayesM$coeff, fit.bayesM$p.value, results05M)
cbind(colnames(fit.bayesM$genes)) # [c(1:7,9,8,13:15,21)])
colnameM<-c(
colnames(fit.bayesM$genes[c(1:4,9,8,13:15,21)]),
c(paste("lods.M",colnames(fit.bayesM$lods), sep=".")), 
c(paste("coeff.M",colnames(fit.bayesM$coeff), sep=".")), 
c(paste("p.value.M",colnames(fit.bayesM$p.value), sep=".")), 
c(paste("n005.M",colnames(results005M), sep=".")))  ## switch colnames back on TC 2009 TC 2011
#
colnames(writeM.table)<-colnameM
# write.table(writeM.table,file="EctodiniMaleCollapse005_20140206filterfaintwMT104.tab",  row.names = FALSE, sep="\t")
# FEMALE
## make a table to merge and keep results around
writeF.table <- cbind(fit.bayesF$genes[, c(1:4,8,9,13:15,21), ], fit.bayesF$lods, fit.bayesF$coeff, fit.bayesF$p.value, results005F)
#writeF.table <- cbind(fit.bayesF$genes[, c(1:4,8,9,13:15,21), ], fit.bayesF$lods, fit.bayesF$coeff, fit.bayesF$p.value, resultsfdr05F)
#writeF.table <- cbind(fit.bayesF$genes[, c(1:4,8,9,13:15,21), ], fit.bayesF$lods, fit.bayesF$coeff, fit.bayesF$p.value, results05F)
colnameF<-c(
colnames(fit.bayesF$genes[c(1:4,9,8,13:15,21)]),
c(paste("lods.F",colnames(fit.bayesF$lods), sep=".")), 
c(paste("coeff.F",colnames(fit.bayesF$coeff), sep=".")), 
c(paste("p.value.F",colnames(fit.bayesF$p.value), sep=".")), 
c(paste("n005.F",colnames(results005F), sep="."))) ## switch colnames back on TC 2009 TC 2011
#
colnames(writeF.table)<-colnameF
#which(writeF.table[,55]==1 & writeF.table[,53]==-1 & writeF.table[,51]==-1)  ## 3 
dim(writeF.table) #17118    62
cbind(colnames(writeF.table)) # just to check
# write.table(writeF.table,file="EctodiniFemaleCollapse005_20130903fileterfaint.tab",  row.names = FALSE, sep="\t")

## MALE INDIVIDUAL
writeMind.table <- cbind(fit.bayesMind$genes[, c(1:4,8,9,13:15,21), ], fit.bayesMind$lods, fit.bayesMind$coeff, fit.bayesMind$p.value, fit.bayesMind)
colnames(writeMind.table)  ### don't really need to change em, ID# indicates male
write.table(writeMind.table,file="EctodiniMaleINDCollapse005_20140213fileterfaintWmt104.tab",  row.names = FALSE, sep="\t")
 ### NOTE SUZY THIS TABLE IS SCREWY CUZ YOU HAVE fit.bayesMind at the end so everything is in there twice. same for female
## FEMALE INDIVIDUAL
writeFind.table <- cbind(fit.bayesFind$genes[, c(1:4,8,9,13:15,21), ], fit.bayesFind$lods, fit.bayesFind$coeff, fit.bayesFind$p.value, fit.bayesFind)
write.table(writeFind.table,file="EctodiniFemaleINDCollapse005_20140213fileterfaintWmt104.tab",  row.names = FALSE, sep="\t")


#################################################################################







#################################################################################
#MERGE MALE AND FEMALE by cbind and keep only one set of gene names
## check that ID is the same
which(writeM.table[,4]!=writeF.table[,4]) # integer(0) 
BOTH<- cbind (writeM.table[ ,c(1:62)], writeF.table[ , c(11:62)])
cbind(colnames(BOTH)) # looks good
dim(BOTH) #[1] 17118   114
length(colnames(BOTH)) # 114
colnames(BOTH)
write.table(BOTH,file="EctodiniBOTHCollapse005_20140206filterfaintwMt104.tab",  row.names = FALSE, sep="\t", na = "NA", col.names = TRUE)
data <- read.table(file="EctodiniBOTHCollapse005_20140206filterfaintwMt104.tab",   header = TRUE)

# write.table(BOTH,file="EctodiniBOTHCollapse05_20131009filterfaint.tab",  row.names = FALSE, sep="\t", na = "NA", col.names = TRUE)
# data <- read.table(file="EctodiniBOTHCollapse05_20131009filterfaint.tab",   header = TRUE)

dim(data) # [1] 17118   114
cbind(colnames(data))


#################################################################################
## LIMIT to survived in both ::::: what is the point in doing this?
dataB <- data[which(!is.na(data$coeff.F.AlvsMt) &!is.na(data$coeff.M.AlvsMt )), ]
length(which(dataB$n005.F.AlvsMt != 0))  # 848 this checks your looking at 005
dim(dataB) # [1] 5232   114  # number survived in both
 which(dataB[,107]==1 & dataB[,103]==-1 & dataB[,105]==-1)  ## all 3 survived Xf up (again a check)
cbind(colnames(dataB))
 # [11,] "lods.M.AlvsMt"         	 [24,] "coeff.M.AlvsMt"        	 [37,] "p.value.M.AlvsMt"      	 [50,] "n005.M.AlvsMt"         
 # [12,] "lods.M.AlvsXf"         	 [25,] "coeff.M.AlvsXf"        	 [38,] "p.value.M.AlvsXf"      	 [51,] "n005.M.AlvsXf"         
 # [13,] "lods.M.AlvsXo"         	 [26,] "coeff.M.AlvsXo"        	 [39,] "p.value.M.AlvsXo"      	 [52,] "n005.M.AlvsXo"         
 # [14,] "lods.M.MtvsXf"         	 [27,] "coeff.M.MtvsXf"        	 [40,] "p.value.M.MtvsXf"      	 [53,] "n005.M.MtvsXf"         
 # [15,] "lods.M.MtvsXo"         	 [28,] "coeff.M.MtvsXo"        	 [41,] "p.value.M.MtvsXo"      	 [54,] "n005.M.MtvsXo"         
 # [16,] "lods.M.XfvsXo"         	 [29,] "coeff.M.XfvsXo"        	 [42,] "p.value.M.XfvsXo"      	 [55,] "n005.M.XfvsXo"         
 # [17,] "lods.M.XfXovsAlMt"     	 [30,] "coeff.M.XfXovsAlMt"    	 [43,] "p.value.M.XfXovsAlMt"  	 [56,] "n005.M.XfXovsAlMt"     
 # [18,] "lods.M.XfAlvsXoMt"     	 [31,] "coeff.M.XfAlvsXoMt"    	 [44,] "p.value.M.XfAlvsXoMt"  	 [57,] "n005.M.XfAlvsXoMt"     
 # [19,] "lods.M.AlvsAll"        	 [32,] "coeff.M.AlvsAll"       	 [45,] "p.value.M.AlvsAll"     	 [58,] "n005.M.AlvsAll"        
 # [20,] "lods.M.XfvsAll"        	 [33,] "coeff.M.XfvsAll"       	 [46,] "p.value.M.XfvsAll"     	 [59,] "n005.M.XfvsAll"        
 # [21,] "lods.M.MtvsAll"        	 [34,] "coeff.M.MtvsAll"       	 [47,] "p.value.M.MtvsAll"     	 [60,] "n005.M.MtvsAll"        
 # [22,] "lods.M.XovsAll"        	 [35,] "coeff.M.XovsAll"       	 [48,] "p.value.M.XovsAll"     	 [61,] "n005.M.XovsAll"        
 # [23,] "lods.M.XfXoDiffAlMt"   	 [36,] "coeff.M.XfXoDiffAlMt"  	 [49,] "p.value.M.XfXoDiffAlMt"	 [62,] "n005.M.XfXoDiffAlMt"   

 # [63,] "lods.F.AlvsMt"         	 [76,] "coeff.F.AlvsMt"        	 [89,] "p.value.F.AlvsMt"      	[102,] "n005.F.AlvsMt"         
 # [64,] "lods.F.AlvsXf"         	 [77,] "coeff.F.AlvsXf"        	 [90,] "p.value.F.AlvsXf"      	[103,] "n005.F.AlvsXf"         
 # [65,] "lods.F.AlvsXo"         	 [78,] "coeff.F.AlvsXo"        	 [91,] "p.value.F.AlvsXo"      	[104,] "n005.F.AlvsXo"         
 # [66,] "lods.F.MtvsXf"         	 [79,] "coeff.F.MtvsXf"        	 [92,] "p.value.F.MtvsXf"      	[105,] "n005.F.MtvsXf"         
 # [67,] "lods.F.MtvsXo"         	 [80,] "coeff.F.MtvsXo"        	 [93,] "p.value.F.MtvsXo"      	[106,] "n005.F.MtvsXo"         
 # [68,] "lods.F.XfvsXo"         	 [81,] "coeff.F.XfvsXo"        	 [94,] "p.value.F.XfvsXo"      	[107,] "n005.F.XfvsXo"         
 # [69,] "lods.F.XfXovsAlMt"     	 [82,] "coeff.F.XfXovsAlMt"    	 [95,] "p.value.F.XfXovsAlMt"  	[108,] "n005.F.XfXovsAlMt"     
 # [70,] "lods.F.XfAlvsXoMt"     	 [83,] "coeff.F.XfAlvsXoMt"    	 [96,] "p.value.F.XfAlvsXoMt"  	[109,] "n005.F.XfAlvsXoMt"     
 # [71,] "lods.F.AlvsAll"        	 [84,] "coeff.F.AlvsAll"       	 [97,] "p.value.F.AlvsAll"     	[110,] "n005.F.AlvsAll"        
 # [72,] "lods.F.XfvsAll"        	 [85,] "coeff.F.XfvsAll"       	 [98,] "p.value.F.XfvsAll"     	[111,] "n005.F.XfvsAll"        
 # [73,] "lods.F.MtvsAll"        	 [86,] "coeff.F.MtvsAll"       	 [99,] "p.value.F.MtvsAll"     	[112,] "n005.F.MtvsAll"        
 # [74,] "lods.F.XovsAll"        	 [87,] "coeff.F.XovsAll"       	[100,] "p.value.F.XovsAll"     	[113,] "n005.F.XovsAll"        
 # [75,] "lods.F.XfXoDiffAlMt"   	 [88,] "coeff.F.XfXoDiffAlMt"  	[101,] "p.value.F.XfXoDiffAlMt"	[114,] "n005.F.XfXoDiffAlMt" 

 ### USE THE FOLLOWING IF YOU WANT ONLY THOSE THAT SURVIVED IN BOTH (1st # is 005 no FDR no faint filter)(2nd number is fdr 05 w/o faint filter) (3rd is 005 no FDR w/ faint filter) (4th # is 05 no FDR w/ faint filter)
length(which(dataB$n005.F.AlvsMt ==1)) #   583   # 522  # 412 #749
length(which(dataB$n005.F.AlvsXf ==1)) #      93 # 11 # 54	#195
length(which(dataB$n005.F.AlvsXo ==1)) #      168 # 30 # 88	#332
length(which(dataB$n005.F.MtvsXf ==1)) #      277 # 239 # 169 #490
length(which(dataB$n005.F.MtvsXo  ==1)) #     467 # 456 # 296 #626
length(which(dataB$n005.F.XfvsXo ==1)) #      99 # 20 # 70 #331
length(which(dataB$n005.F.XfXovsAlMt==1)) #   368 #238 # 274	#624
length(which(dataB$n005.F.XfAlvsXoMt ==1)) #  373 # 162 # 260	#677
length(which(dataB$n005.F.AlvsAll ==1)) #     127 # 32 # 67	#225
length(which(dataB$n005.F.XfvsAll==1)) #      212 # 20 # 156	#583
length(which(dataB$n005.F.MtvsAll ==1)) #     518 # 572 # 327 #670
length(which(dataB$n005.F.XovsAll ==1)) #     219 # 112 # 166  #446
length(which(dataB$n005.F.XfXoDiffAlMt ==1)) #230  # 160 # 170 #474
# (#s here are jus 005 no FDR w/ faint filter)
length(which(dataB$n005.F.AlvsMt ==1 & dataB$n005.F.XfAlvsXoMt ==1)) #225 (1/2 of the AlMt regulation)
length(which(dataB$n005.F.XfvsXo ==1 & dataB$n005.F.XfAlvsXoMt ==1)) #22 (this is about 1/3rd of the XfXo regulation)
																	# but see only 7 are in both

length(which(dataB$n005.F.AlvsMt ==-1)) #      265	# 223 #181 #465
length(which(dataB$n005.F.AlvsXf ==-1)) #      108 #	2 # 70 #317
length(which(dataB$n005.F.AlvsXo ==-1)) #      128 # 23  # 71 #259
length(which(dataB$n005.F.MtvsXf ==-1)) #      633 # 575 # 457  # 846
length(which(dataB$n005.F.MtvsXo  ==-1)) #     574 # 558 # 399  # 731
length(which(dataB$n005.F.XfvsXo ==-1)) #      68 # 16 # 33  # 121
length(which(dataB$n005.F.XfXovsAlMt==-1)) #   244 # 122  # 157 #443
length(which(dataB$n005.F.XfAlvsXoMt ==-1)) #  102 # 37 # 71  # 289
length(which(dataB$n005.F.AlvsAll ==-1)) #     84 # 5 # 42 #192
length(which(dataB$n005.F.XfvsAll==-1)) #      95 # 9 # 62 # 222
length(which(dataB$n005.F.MtvsAll ==-1)) #     821 # 879  ### check the 005 number, how can it be less than the 05 fdr? # 565 #914
length(which(dataB$n005.F.XovsAll ==-1)) #     264 # 122 # 173 #471
length(which(dataB$n005.F.XfXoDiffAlMt ==-1)) #301 # 203 # 206 #451
(1st # is 005 no FDR no faint filter)(2nd number is fdr 05 w/o faint filter) (3rd is 005 no FDR w/ faint filter) (4th # is 05 no FDR w/ faint filter) (## 5TH NUMBER IS WITH MT104 INCLUDED)
length(which(dataB$n005.M.AlvsMt ==1)) #      525 # 446 # 126 # 437  ## 146
length(which(dataB$n005.M.AlvsXf ==1)) #      275 # 106 # 70 # 426	## 70
length(which(dataB$n005.M.AlvsXo ==1)) #      139 # 0 # 16 # 212	## 16	
length(which(dataB$n005.M.MtvsXf ==1)) #      389 # 377 # 109 #421	## 95
length(which(dataB$n005.M.MtvsXo  ==1)) #     294 # 227  # 67 # 253	## 66
length(which(dataB$n005.M.XfvsXo ==1)) #      134 # 39 # 47 # 363	## 47
length(which(dataB$n005.M.XfXovsAlMt==1)) #   347 # 202 # 75 #329	## 92
length(which(dataB$n005.M.XfAlvsXoMt ==1)) #  473 # 395 # 163 # 673	## 131
length(which(dataB$n005.M.AlvsAll ==1)) #     195 # 13 # 36 # 277	## 38
length(which(dataB$n005.M.XfvsAll==1)) #      466 # 359 # 127 # 629	## 133
length(which(dataB$n005.M.MtvsAll ==1)) #     380 # 362  # 91 #315	## 88
length(which(dataB$n005.M.XovsAll ==1)) #     189 # 3 # 15 # 123	## 22
length(which(dataB$n005.M.XfXoDiffAlMt ==1)) #224 # 157 # 41 #176	## 54
 
length(which(dataB$n005.M.AlvsMt ==-1)) # 	298    # 254 # 57 # 172	## 61
length(which(dataB$n005.M.AlvsXf ==-1)) #    274   # 118 # 57 # 352	## 57
length(which(dataB$n005.M.AlvsXo ==-1)) #      151 # 0  # 19 # 131 	## 19
length(which(dataB$n005.M.MtvsXf ==-1)) #      626 # 616 # 179 # 648## 	163
length(which(dataB$n005.M.MtvsXo  ==-1)) #     470 # 370 # 102 # 326## 	99
length(which(dataB$n005.M.XfvsXo ==-1)) #      245 # 135 # 44 # 236	## 44
length(which(dataB$n005.M.XfXovsAlMt==-1)) #   258 # 122 # 74 #423 	## 69
length(which(dataB$n005.M.XfAlvsXoMt ==-1)) #  288 # 230 # 84 #509 	## 46
46 + 131 # = 177 (WAS 247) total mating strategty
length(which(dataB$n005.M.AlvsAll ==-1)) #     171 # 22 # 32 # 152	## 31 
length(which(dataB$n005.M.XfvsAll==-1)) #      332 # 264  # 102 # 522 ##  97		
length(which(dataB$n005.M.MtvsAll ==-1)) #     635 # 622 # 170 # 537 ##  166		
length(which(dataB$n005.M.XovsAll ==-1)) #     111 # 1  # 19 # 192	## 19	
length(which(dataB$n005.M.XfXoDiffAlMt ==-1)) #366 # 243 # 69 #253	## 81	
 
#MALEn005 w/o faint filter  (2nd number is fdr 05 w/o faint filter) (3rd is w/ faint filter) (4th # is 05 no FDR w/ faint filter) (## 5TH NUMBER IS WITH MT104 INCLUDED)
length(which(dataB$n005.M.AlvsMt==1 & dataB$n005.M.XfvsXo==1)) #16 Monogamy	# 7 # 12 #84 ## 7
length(which(dataB$n005.M.AlvsMt==-1 & dataB$n005.M.XfvsXo==-1)) #26 polygamy 	# 18 # 7 #19	## 6
length(which(dataB$n005.M.AlvsMt==1 & dataB$n005.M.XfvsXo==-1)) #28 twisted Al/Xo # 15 # 6 # 22	## 8
length(which(dataB$n005.M.AlvsMt==-1 & dataB$n005.M.XfvsXo==1)) #12 twisted Mt/Xf # 4 # 0 #14	## 2
length(which(dataB$n005.M.AlvsMt==1 & dataB$n005.M.AlvsXf==1 & dataB$n005.M.AlvsXo==1 )) # 13 Al up # 0 # 1 #18 ## 1
length(which(dataB$n005.M.AlvsMt==-1 & dataB$n005.M.AlvsXf==-1 & dataB$n005.M.AlvsXo==-1 )) # 11 Al down # 0 # 2 # 19	## 2
length(which(dataB$n005.M.MtvsXf==1 & dataB$n005.M.MtvsXo==1 & dataB$n005.M.AlvsMt==-1 )) # 104 Mt up # 89 # 24 #74	## 25
length(which(dataB$n005.M.MtvsXf==-1 & dataB$n005.M.MtvsXo==-1 & dataB$n005.M.AlvsMt==1 )) # 232 Mt down # 191 # 46 #194 ## 50
length(which(dataB$n005.M.XfvsXo==1 & dataB$n005.M.AlvsXf==-1 & dataB$n005.M.MtvsXf==-1 )) # 34 Xf up # 11 #22 #151	## 20
length(which(dataB$n005.M.XfvsXo==-1 & dataB$n005.M.AlvsXf==1 & dataB$n005.M.MtvsXf==1 )) # 34 Xf down # 20 # 7 #87	## 5
length(which(dataB$n005.M.XfvsXo==-1 & dataB$n005.M.AlvsXo==-1 & dataB$n005.M.MtvsXo==-1 )) # 8 Xo up # 0 # 0 #5	## 0
length(which(dataB$n005.M.XfvsXo==1 & dataB$n005.M.AlvsXo==1 & dataB$n005.M.MtvsXo==1 )) # 1 Xo down # 0 # 0 # 1	## 0

#FEMALE
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.F.XfvsXo==1))# 11 Monogamy # 2 # 7 #49
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.F.XfvsXo==-1)) # 2 polygamy # 0 # 2 #13
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.F.XfvsXo==-1)) # 18 twisted AL/Xo # 4 # 7 # 35
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.F.XfvsXo==1))# 11 twisted Mt/Xf # 1 # 8 #72
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.F.AlvsXf==1 & dataB$n005.F.AlvsXo==1 )) # 21 Al up # 3 # 14 #43
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.F.AlvsXf==-1 & dataB$n005.F.AlvsXo==-1 )) # 6 Al down # 0 # 7 # 31
length(which(dataB$n005.F.MtvsXf==1 & dataB$n005.F.MtvsXo==1 & dataB$n005.F.AlvsMt==-1 )) # 102 Mt up # 86 # 75 # 227
length(which(dataB$n005.F.MtvsXf==-1 & dataB$n005.F.MtvsXo==-1 & dataB$n005.F.AlvsMt==1 )) # 289 Mt down  # 268 # 232 #446
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.F.AlvsXf==-1 & dataB$n005.F.MtvsXf==-1 )) # 5 Xf up # 0 # 3  # 42
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.F.AlvsXf==1 & dataB$n005.F.MtvsXf==1 )) # 4 Xf down # 1 # 3 #14
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.F.AlvsXo==-1 & dataB$n005.F.MtvsXo==-1 )) # 14 Xo up # 4 # 2 # 21
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.F.AlvsXo==1 & dataB$n005.F.MtvsXo==1 )) # 23 Xo down # 3 # 8 #89

# BOTH
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.F.XfvsXo==1 & dataB$n005.M.AlvsMt==1 & dataB$n005.M.XfvsXo==1))# 0 Monogamy # 0 # 0 # 3 ## 0
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.F.XfvsXo==-1 & dataB$n005.M.AlvsMt==-1 & dataB$n005.M.XfvsXo==-1)) # 1 polygamy # 0 # 0 # 0 ## 0
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.F.XfvsXo==-1 & dataB$n005.M.AlvsMt==1 & dataB$n005.M.XfvsXo==-1)) # 4 twisted AL/Xo # 0 # 1 # 3 ## 1
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.F.XfvsXo==1 & dataB$n005.M.AlvsMt==-1 & dataB$n005.M.XfvsXo==1))# 0 twisted Mt/Xf # 0 # 0 # 1 ## 0
length(which(dataB$n005.F.XfAlvsXoMt==1 & dataB$n005.M.XfAlvsXoMt==1)) # monogamy 180 # 86 # 71 # 269	## 66
length(which(dataB$n005.F.XfAlvsXoMt==-1 & dataB$n005.M.XfAlvsXoMt==-1)) # polygamy 47 # 23 # 17 # 93	## 15
length(which(dataB$n005.F.XfAlvsXoMt==1 & dataB$n005.M.XfAlvsXoMt==-1)) # twisted 3 # 1 # 2 # 13	## 1
length(which(dataB$n005.F.XfAlvsXoMt==-1 & dataB$n005.M.XfAlvsXoMt==1)) # twisted 2 #2 # 0 # 9	## 0
length(which(dataB$n005.F.XfXovsAlMt==1 & dataB$n005.M.XfXovsAlMt==1)) # X-lineage up 162 #81 # 35 # 183	## 46
length(which(dataB$n005.F.XfXovsAlMt==-1 & dataB$n005.M.XfXovsAlMt==-1)) # X-lineage down 95 # 36 # 30 # 156	## 33
length(which(dataB$n005.F.XfXovsAlMt==1 & dataB$n005.M.XfXovsAlMt==-1)) # lineage - reversed 0 # 0 # 0 # 5	## 0
length(which(dataB$n005.F.XfXovsAlMt==-1 & dataB$n005.M.XfXovsAlMt==1))# lineage - reversed 0  # 0 # 0 #2	## 0
length(which(dataB$n005.F.XfXovsAlMt==1 & dataB$coeff.M.XfXovsAlMt<0)) #non-sig reversal 7 # 4 # 34 #136	## 19
length(which(dataB$n005.F.XfXovsAlMt==-1 & dataB$coeff.M.XfXovsAlMt>0)) # non-sig reversal 11 # 4 # 24 #77	## 12
# (numbers for nex two are just 005 no FDR w/ faint filter) (then 05 no fdr w/ faint filter)	
length(which(dataB$n005.M.XfXovsAlMt==1 & dataB$coeff.F.XfXovsAlMt<0)) #non-sig reversal # 4 # 47	## 5
length(which(dataB$n005.M.XfXovsAlMt==-1 & dataB$coeff.F.XfXovsAlMt>0)) # non-sig reversal # 1 #64	## 0

mean(c(dataB$coeff.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt==1 & dataB$coeff.M.XfXovsAlMt<0)])) # -0.09571475 # -0.1322 # -0.2118 -0.2449	#-0.1447318
mean(dataB$coeff.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt==-1 & dataB$coeff.M.XfXovsAlMt>0)]) #  0.07347005 # 0.0608 # 0.20106 0.2210 ## 1.03069345
mean(abs(c(dataB$coeff.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt==1 & dataB$coeff.M.XfXovsAlMt<0)], dataB$coeff.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt==-1 & dataB$coeff.M.XfXovsAlMt>0)]))) # 0.0821  AT log2 this is 1.05X difference ## 0.09614 0.20736 # 0.23628 ## 0.2075199

## BOTH  other way of counting "species specific" (numbers just for w/faint filter 005 no FDR) (4th # is 05 no FDR w/ faint filter) (## 5TH NUMBER IS WITH MT104 INCLUDED)
length(which(dataB$n005.M.AlvsAll ==1 & dataB$n005.F.AlvsAll ==1)) # Al up 9 #64	## 11
length(which(dataB$n005.M.XfvsAll ==1 & dataB$n005.F.XfvsAll ==1)) # Xf up 34   # 220	## 35
length(which(dataB$n005.M.MtvsAll ==1 & dataB$n005.F.MtvsAll ==1)) # Mt up  67   # 211	## 66
length(which(dataB$n005.M.XovsAll ==1 & dataB$n005.F.XovsAll ==1)) # Xo up  5 # 65	## 7
length(which(dataB$n005.M.AlvsAll ==-1 & dataB$n005.F.AlvsAll ==-1)) # Al down 2  # 27	## 3
length(which(dataB$n005.M.XfvsAll ==-1 & dataB$n005.F.XfvsAll ==-1)) # Xf down  14  # 83	## 15
length(which(dataB$n005.M.MtvsAll ==-1 & dataB$n005.F.MtvsAll ==-1)) # Mt down  130   # 366	## 127
length(which(dataB$n005.M.XovsAll ==-1 & dataB$n005.F.XovsAll ==-1)) # Xo down  10 #110	## 10


length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.F.AlvsXf==1 & dataB$n005.F.AlvsXo==1  & dataB$n005.M.AlvsMt==1 & dataB$n005.M.AlvsXf==1 & dataB$n005.M.AlvsXo==1)) # 2 Al up # 0 # 0 # 2	## 0
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.F.AlvsXf==-1 & dataB$n005.F.AlvsXo==-1 & dataB$n005.M.AlvsMt==-1 & dataB$n005.M.AlvsXf==-1 & dataB$n005.M.AlvsXo==-1)) # 1 Al down # 0 # 0 # 1	## 0
length(which(dataB$n005.F.MtvsXf==1 & dataB$n005.F.MtvsXo==1 & dataB$n005.F.AlvsMt==-1 & dataB$n005.M.MtvsXf==1 & dataB$n005.M.MtvsXo==1 & dataB$n005.M.AlvsMt==-1)) # 39 Mt up shared  # 29 # 10 # 48		## 11
length(which(dataB$n005.F.MtvsXf==-1 & dataB$n005.F.MtvsXo==-1 & dataB$n005.F.AlvsMt==1 & dataB$n005.M.MtvsXf==-1 & dataB$n005.M.MtvsXo==-1 & dataB$n005.M.AlvsMt==1)) # 117 Mt down shared # 97 # 33 #136	## 34
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.F.AlvsXf==-1 & dataB$n005.F.MtvsXf==-1 & dataB$n005.M.XfvsXo==1 & dataB$n005.M.AlvsXf==-1 & dataB$n005.M.MtvsXf==-1)) # 0 Xf up # 0 # 0 # 1	## 0
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.F.AlvsXf==1 & dataB$n005.F.MtvsXf==1 & dataB$n005.M.XfvsXo==-1 & dataB$n005.M.AlvsXf==1 & dataB$n005.M.MtvsXf==1)) # 1 Xf down # 0 # 0 # 2	## 0
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.F.AlvsXo==-1 & dataB$n005.F.MtvsXo==-1 & dataB$n005.M.XfvsXo==-1 & dataB$n005.M.AlvsXo==-1 & dataB$n005.M.MtvsXo==-1  )) # 1 Xo up # 0 # 0 # 0	## 0
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.F.AlvsXo==1 & dataB$n005.F.MtvsXo==1 & dataB$n005.M.XfvsXo==1 & dataB$n005.M.AlvsXo==1 & dataB$n005.M.MtvsXo==1 )) # 0 Xo down # 0 # 0 # 0	## 0
 
# PAIRWISE BOTH AND REVERSED
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.M.AlvsMt==1 )) # 275  #23 # 79 # 285	## 81
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.M.AlvsMt==-1 )) #107 # 87 # 30 # 95 	## 28
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.M.AlvsMt==-1 )) # 2 # 1 # 0 # 2	## 0
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.M.AlvsMt==1 )) # 2 # 1 # 0 # 9	## 0
length(which(dataB$n005.F.AlvsXf==1 & dataB$n005.M.AlvsXf==1 )) # 32 # 5 # 9 # 67	## 9
length(which(dataB$n005.F.AlvsXf==-1 & dataB$n005.M.AlvsXf==-1 )) # 44  # 0 # 7 #79	## 7
length(which(dataB$n005.F.AlvsXf==1 & dataB$n005.M.AlvsXf==-1 )) # 0 # 0 # 0 # 2	## 0
length(which(dataB$n005.F.AlvsXf==-1 & dataB$n005.M.AlvsXf==1 )) # 0 # 0 # 0 # 4	## 0
length(which(dataB$n005.F.AlvsXo==1 & dataB$n005.M.AlvsXo==1 )) # 35 # 0 # 5 # 94	## 5
length(which(dataB$n005.F.AlvsXo==-1 & dataB$n005.M.AlvsXo==-1 )) # 29 # 0 # 4 # 37	## 4
length(which(dataB$n005.F.AlvsXo==1 & dataB$n005.M.AlvsXo==-1 ))  # 8 # 0 # 2 # 7 	## 2
length(which(dataB$n005.F.AlvsXo==-1 & dataB$n005.M.AlvsXo==1 )) # 0 # 0 # 0# 1	0
length(which(dataB$n005.F.MtvsXf==1 & dataB$n005.M.MtvsXf==1 )) # 145 # 132 # 42 #183	## 45
length(which(dataB$n005.F.MtvsXf==-1 & dataB$n005.M.MtvsXf==-1 )) # 378 # 351 # 117 # 375	## 106
length(which(dataB$n005.F.MtvsXf==1 & dataB$n005.M.MtvsXf==-1 )) # 0 # 0 # 1 # 11	## 0
length(which(dataB$n005.F.MtvsXf==-1 & dataB$n005.M.MtvsXf==1 )) # 0 # 0 # 0 # 3	## 0
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.M.XfvsXo==1  )) # 5 # 1 # 1 # 23	## 1
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.M.XfvsXo==-1  )) # 22 # 3 # 7 # 28	## 7
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.M.XfvsXo==-1  )) #11 # 6 # 6 # 16	## 6
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.M.XfvsXo==1  )) #2 # 2 # 2 # 7	## 2


#### INTERACTION
length ( which((dataB$n005.F.XfXoDiffAlMt ==1 | dataB$n005.F.XfXoDiffAlMt ==-1)  & 
         dataB$coeff.F.XfvsXo > 0 & dataB$coeff.F.AlvsMt > 0)) # 65 interaction with both monogamy biased  #37 # 42 # 124 ## 42
length ( which((dataB$n005.F.XfXoDiffAlMt ==1 | dataB$n005.F.XfXoDiffAlMt ==-1)  &  
         dataB$coeff.F.XfvsXo < 0 & dataB$coeff.F.AlvsMt < 0))	# 16 interaction with both polygamy biased # 8 # 11 # 44 ## 11
length ( which((dataB$n005.F.XfXoDiffAlMt ==1 | dataB$n005.F.XfXoDiffAlMt ==-1)  & 
         dataB$coeff.F.XfvsXo > 0 & dataB$coeff.F.AlvsMt < 0)) # 214 interaction, monogamy biased in X only (not thresholeded) # 111 #161 # 430 ## 161
length ( which((dataB$n005.F.XfXoDiffAlMt ==1 | dataB$n005.F.XfXoDiffAlMt ==-1)  & 
         dataB$n005.F.XfvsXo == 1 & dataB$n005.F.AlvsMt ==-1 )) # 11 of theh 214 are significantly regulated in oppos direction #1 # 8 # 72  ## 8 	
length ( which((dataB$n005.F.XfXoDiffAlMt ==1 | dataB$n005.F.XfXoDiffAlMt ==-1)  & 
         dataB$coeff.F.XfvsXo < 0 & dataB$coeff.F.AlvsMt > 0)) # 236 interaction monogamy biased in non-X only (not thresholeded) #163 #182 # 327	##182
length ( which((dataB$n005.F.XfXoDiffAlMt ==1 | dataB$n005.F.XfXoDiffAlMt ==-1)  & 
         dataB$n005.F.XfvsXo == -1 & dataB$n005.F.AlvsMt ==1 )) # 18 of the 236 are significantly regulated in opposing direction # 4 # 7 # 35 ## 7
         


######################
## Plotting  

Xocol <- "darkmagenta"
Xfcol <- "goldenrod3"
Mtcol <- "blue4"
Alcol <- "darkorange3"
nscol <- "gray80"   
                 
# ####### NOT USED IN MANUSCRIPT
# ## plot similarity between sexes for  AlvsMT
# # male on female
 # plot (dataB$coeff.M.AlvsMt, dataB$lods.M.AlvsMt, pch = 17, cex =0.5, col = nscol, ylab = "Lods", xlab = "Log fold change MALE", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.M.AlvsMt[which(dataB$n005.F.AlvsMt ==1)], dataB$lods.M.AlvsMt[which(dataB$n005.F.AlvsMt ==1)], pch = 19, cex =0.5, col = Alcol, ylab = "Lods", xlab = "Log fold change")
  # points (dataB$coeff.M.AlvsMt[which(dataB$n005.F.AlvsMt ==-1)], dataB$lods.M.AlvsMt[which(dataB$n005.F.AlvsMt ==-1)], pch = 19, cex =0.5, col = Mtcol, ylab = "Lods", xlab = "Log fold change")
# legend(x=.4 ,y=12.8 ,  bty="n", pch = c(19,19), pt.lwd=2, legend = c("female A. leptura", "female M. tenuindentata"), col = c(Alcol,Mtcol), cex = 0.8)
# # female on males
 # plot (dataB$coeff.F.AlvsMt, dataB$lods.F.AlvsMt, pch = 19, cex =0.5, col = nscol, ylab = "Lods", xlab = "Log fold change FEMALE", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.AlvsMt ==1)], dataB$lods.F.AlvsMt[which(dataB$n005.M.AlvsMt ==1)], pch = 17, cex =0.7 , col = Alcol, ylab = "Lods", xlab = "Log fold change")
# points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.AlvsMt ==-1)], dataB$lods.F.AlvsMt[which(dataB$n005.M.AlvsMt ==-1)], pch = 17, cex =0.7, col = Mtcol, ylab = "Lods", xlab = "Log fold change")
# legend(x=.4 ,y=12.8 ,  bty="n", pch = c(17,17), pt.lwd=2, legend = c("male A. leptura", "male M. tenuindentata"), col = c(Mtcol,Alcol), cex = 0.8)
 
# ##plot similarity between sexes XfvsXo ##### not so similar
# # male on female
 # plot (dataB$coeff.F.XfvsXo, dataB$lods.F.XfvsXo, pch = 17, cex =0.5, col = nscol,  ylab = "Lods", xlab = "Log fold change FEMALE", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.F.XfvsXo[which(dataB$n005.M.XfvsXo ==1)], dataB$lods.F.XfvsXo[which(dataB$n005.M.XfvsXo ==1)], pch = 19, cex =0.5, col = Alcol)
  # points (dataB$coeff.F.XfvsXo[which(dataB$n005.M.XfvsXo ==-1)], dataB$lods.F.XfvsXo[which(dataB$n005.M.XfvsXo ==-1)], pch = 19, cex =0.5, col = Mtcol)
  # legend(x=.4 ,y=12.8 ,  bty="n", pch = c(19,19), pt.lwd=2, legend = c("male X. flavipinnis", "male X. ochrogenys"), col = c(Alcol,Mtcol), cex = 0.8)
# # female on male
 # plot (dataB$coeff.M.XfvsXo, dataB$lods.M.XfvsXo, pch = 19, cex =0.5, col = nscol,  ylab = "Lods", xlab = "Log fold change MALE", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.M.XfvsXo[which(dataB$n005.F.XfvsXo ==1)], dataB$lods.M.XfvsXo[which(dataB$n005.F.XfvsXo ==1)], pch = 17, cex =0.5, col = Alcol)
  # points (dataB$coeff.M.XfvsXo[which(dataB$n005.F.XfvsXo ==-1)], dataB$lods.M.XfvsXo[which(dataB$n005.F.XfvsXo ==-1)], pch = 17, cex =0.5, col = Mtcol)
  # legend(x=.4 ,y=12.8 ,  bty="n", pch = c(19,19), pt.lwd=2, legend = c("female X. flavipinnis", "female X. ochrogenys"), col = c(Alcol,Mtcol), cex = 0.8)
  
  
####### NOT USED IN MANUSCRIPT
# ##plot similarity between lineages Mon vs Poly for female
# # non-X on X for MAle
# #jpeg("FigureXwithinLineageMating.jpg", width = 8, height = 8, units = "in", quality=100, pointsize=12, res = 300)
# #par(mfrow=c(2,2))
 # plot (dataB$coeff.M.XfvsXo, dataB$lods.M.XfvsXo, pch = 19, cex =0.5, col = nscol,  ylab = "Lods", xlab = "Log fold change X-lineage Male", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt ==1)], dataB$lods.M.XfvsXo[which(dataB$n005.M.AlvsMt ==1)], pch = 17, cex =0.5, col = Alcol, ylab = "Lods")
  # points (dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt ==-1)], dataB$lods.M.XfvsXo[which(dataB$n005.M.AlvsMt ==-1)], pch = 17, cex =0.5, col = Mtcol)
  # legend(x=-0.2 ,y=12.8 ,  bty="n", pch = c(17,17), pt.lwd=2, legend = c("nonX-lineage monogamy-biased", "nonX-lineage polygamy-biased"), col = c(Alcol,Mtcol), cex = 0.7)
 # abline(-2.298,0, col =nscol, lty = 3, lwd = 2)
# #  X on NonX for MAle
 # plot (dataB$coeff.M.AlvsMt, dataB$lods.M.AlvsMt, pch = 19, cex =0.5, col = nscol,  ylab = "Lods", xlab = "Log fold change nonX-lineage Male", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo ==1)], dataB$lods.M.AlvsMt[which(dataB$n005.M.XfvsXo ==1)], pch = 17, cex =0.5, col = Alcol, ylab = "Lods")
  # points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo ==-1)], dataB$lods.M.AlvsMt[which(dataB$n005.M.XfvsXo ==-1)], pch = 17, cex =0.5, col = Mtcol)
  # legend(x=0 ,y=12.8 ,  bty="n", pch = c(17,17), pt.lwd=2, legend = c("X-lineage monogamy-biased", "X-lineage polygamy-biased"), col = c(Alcol,Mtcol), cex = 0.7) 
   # abline(-2.298,0, col =nscol, lty = 3, lwd = 2) 
# # non-X on X for FEMAle
 # plot (dataB$coeff.F.XfvsXo, dataB$lods.F.XfvsXo, pch = 19, cex =0.5, col = nscol,  ylab = "Lods", xlab = "Log fold change X-lineage Female", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.F.XfvsXo[which(dataB$n005.F.AlvsMt ==1)], dataB$lods.F.XfvsXo[which(dataB$n005.F.AlvsMt ==1)], pch = 19, cex =0.5, col = Alcol, ylab = "Lods")
  # points (dataB$coeff.F.XfvsXo[which(dataB$n005.F.AlvsMt ==-1)], dataB$lods.F.XfvsXo[which(dataB$n005.F.AlvsMt ==-1)], pch = 19, cex =0.5, col = Mtcol)
  # legend(x=-0.2 ,y=12.8 ,  bty="n", pch = c(19,19), pt.lwd=2, legend = c("nonX-lineage monogamy-biased", "nonX-lineage polygamy-biased"), col = c(Alcol,Mtcol), cex = 0.7)
 # abline(-2.298,0, col =nscol, lty = 3, lwd = 2)  
  
# #  X on NonX for FEMAle
 # plot (dataB$coeff.F.AlvsMt, dataB$lods.F.AlvsMt, pch = 19, cex =0.5, col = nscol,  ylab = "Lods", xlab = "Log fold change nonX-lineage Female", xlim = c(-2,2), ylim = c(-7,13))
 # points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.XfvsXo ==1)], dataB$lods.F.AlvsMt[which(dataB$n005.F.XfvsXo ==1)], pch = 19, cex =1, col = Alcol, ylab = "Lods")
  # points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.XfvsXo ==-1)], dataB$lods.F.AlvsMt[which(dataB$n005.F.XfvsXo ==-1)], pch = 19, cex =1, col = Mtcol)
  # legend(x=-2 ,y=12.8 ,  bty="n", pch = c(19,19), pt.lwd=2, legend = c("X-lineage monogamy-biased", "X-lineage polygamy-biased"), col = c(Alcol,Mtcol), cex = 1) 
  # abline(-2.298,0, col =nscol, lty = 3, lwd = 2)  
# #dev.off()  
  
#### Plot LINEAGE VOLCANOES  ::::: FIGURE 3
# Male X vs non-X colored for male (A)
 plot (dataB$coeff.M.XfXovsAlMt, dataB$lods.M.XfXovsAlMt, pch = 19, cex =1, cex.axis = 1.8, cex.lab = 1.5, col = nscol,  ylab = "Lods", xlab = "Log fold change X-lineage vs non-X-lineage (Male)", xlim = c(-2,2), ylim = c(-7,9), main= "Male Analysis of Lineage")
 points (dataB$coeff.M.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==1)], dataB$lods.M.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==1)], pch = 4, lwd = 3, cex =1.5, col = "black")
  points (dataB$coeff.M.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==-1)], dataB$lods.M.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==-1)], pch = 1, lwd = 3, cex =1.5, col = "black")
  legend(x="topleft",  bty="n", pch = c(4,1), pt.lwd=3, legend = c("Male X-lineage-biased", "Male non-X-lineage-biased"), col = c("black", "black"), cex = 1.5) 
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 2)
# Male X vs non-X colored for female (using polygamous colors? is this confusing) (B)
 plot (dataB$coeff.M.XfXovsAlMt, dataB$lods.M.XfXovsAlMt, pch = 19, cex =1, cex.axis = 1.8, cex.lab = 1.5,col = nscol,  ylab = "Lods",  xlab = "Log fold change X-lineage vs non-X-lineage (Male)", xlim = c(-2,2), ylim = c(-7,9), main= "Male Analysis of Lineage")
 points (dataB$coeff.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==1)], dataB$lods.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==1)], pch = 4, lwd = 3, cex =1.5, col = "black")
  points (dataB$coeff.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==-1)], dataB$lods.M.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==-1)], pch = 1, lwd = 3, cex =1.5, col = "black")
  legend(x="topleft",  bty="n", pch = c(4,1), pt.lwd=3, legend = c("Female X-lineage-biased", "Female non-X-lineage-biased"), col = c("black","black"), cex = 1.5) 
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 2)  
  
# Female X vs non-X colored for female (using polygamous colors? is this confusing)(C)
 plot (dataB$coeff.F.XfXovsAlMt, dataB$lods.F.XfXovsAlMt, pch = 19, cex =1, cex.axis = 1.8, cex.lab = 1.5,col = nscol,  ylab = "Lods",  xlab = "Log fold change X-lineage vs non-X-lineage (Female)", xlim = c(-2,2), ylim = c(-7,9), main = "Female Analysis of Lineage")
 points (dataB$coeff.F.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==1)], dataB$lods.F.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==1)], pch = 4, lwd = 3, cex =1.5, col = "black")
  points (dataB$coeff.F.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==-1)], dataB$lods.F.XfXovsAlMt[which(dataB$n005.F.XfXovsAlMt ==-1)], pch = 1, lwd = 3, cex =1.5, col = "black")
  legend(x="topleft",  bty="n", pch = c(4,1), pt.lwd=3, legend = c("Female X-lineage-biased", "Female non-X-lineage-biased"), col = c("black", "black"), cex = 1.5) 
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 2)
# Female X vs non-X colored for male (using polygamous colors? is this confusing)(D)
  plot (dataB$coeff.F.XfXovsAlMt, dataB$lods.F.XfXovsAlMt, pch = 19, cex =1, cex.axis = 1.8, cex.lab = 1.5,col = nscol,  ylab = "Lods",  xlab = "Log fold change X-lineage vs non-X-lineage (Female)", xlim = c(-2,2), ylim = c(-7,9), main = "Female Analysis of Lineage")
 points (dataB$coeff.F.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==1)], dataB$lods.F.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==1)], pch = 4, lwd = 3, cex =1.5, col = "black")
  points (dataB$coeff.F.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==-1)], dataB$lods.F.XfXovsAlMt[which(dataB$n005.M.XfXovsAlMt ==-1)], pch = 1, lwd = 3, cex =1.5, col = "black")
  legend(x="topleft",  bty="n", pch = c(4,1), pt.lwd=3, legend = c("Male X-lineage-biased", "Male non-X-lineage-biased"), col = c("black", "black"), cex = 1.5) 
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 2)   





  
 # [11,] "lods.M.AlvsMt"         	 [24,] "coeff.M.AlvsMt"        	 [37,] "p.value.M.AlvsMt"      	 [50,] "n005.M.AlvsMt"         
 # [12,] "lods.M.AlvsXf"         	 [25,] "coeff.M.AlvsXf"        	 [38,] "p.value.M.AlvsXf"      	 [51,] "n005.M.AlvsXf"         
 # [13,] "lods.M.AlvsXo"         	 [26,] "coeff.M.AlvsXo"        	 [39,] "p.value.M.AlvsXo"      	 [52,] "n005.M.AlvsXo"         
 # [14,] "lods.M.MtvsXf"         	 [27,] "coeff.M.MtvsXf"        	 [40,] "p.value.M.MtvsXf"      	 [53,] "n005.M.MtvsXf"         
 # [15,] "lods.M.MtvsXo"         	 [28,] "coeff.M.MtvsXo"        	 [41,] "p.value.M.MtvsXo"      	 [54,] "n005.M.MtvsXo"         
 # [16,] "lods.M.XfvsXo"         	 [29,] "coeff.M.XfvsXo"        	 [42,] "p.value.M.XfvsXo"      	 [55,] "n005.M.XfvsXo"         
 # [17,] "lods.M.XfXovsAlMt"     	 [30,] "coeff.M.XfXovsAlMt"    	 [43,] "p.value.M.XfXovsAlMt"  	 [56,] "n005.M.XfXovsAlMt"     
 # [18,] "lods.M.XfAlvsXoMt"     	 [31,] "coeff.M.XfAlvsXoMt"    	 [44,] "p.value.M.XfAlvsXoMt"  	 [57,] "n005.M.XfAlvsXoMt"     
 # [19,] "lods.M.AlvsAll"        	 [32,] "coeff.M.AlvsAll"       	 [45,] "p.value.M.AlvsAll"     	 [58,] "n005.M.AlvsAll"        
 # [20,] "lods.M.XfvsAll"        	 [33,] "coeff.M.XfvsAll"       	 [46,] "p.value.M.XfvsAll"     	 [59,] "n005.M.XfvsAll"        
 # [21,] "lods.M.MtvsAll"        	 [34,] "coeff.M.MtvsAll"       	 [47,] "p.value.M.MtvsAll"     	 [60,] "n005.M.MtvsAll"        
 # [22,] "lods.M.XovsAll"        	 [35,] "coeff.M.XovsAll"       	 [48,] "p.value.M.XovsAll"     	 [61,] "n005.M.XovsAll"        
 # [23,] "lods.M.XfXoDiffAlMt"   	 [36,] "coeff.M.XfXoDiffAlMt"  	 [49,] "p.value.M.XfXoDiffAlMt"	 [62,] "n005.M.XfXoDiffAlMt"   

 # [63,] "lods.F.AlvsMt"         	 [76,] "coeff.F.AlvsMt"        	 [89,] "p.value.F.AlvsMt"      	[102,] "n005.F.AlvsMt"         
 # [64,] "lods.F.AlvsXf"         	 [77,] "coeff.F.AlvsXf"        	 [90,] "p.value.F.AlvsXf"      	[103,] "n005.F.AlvsXf"         
 # [65,] "lods.F.AlvsXo"         	 [78,] "coeff.F.AlvsXo"        	 [91,] "p.value.F.AlvsXo"      	[104,] "n005.F.AlvsXo"         
 # [66,] "lods.F.MtvsXf"         	 [79,] "coeff.F.MtvsXf"        	 [92,] "p.value.F.MtvsXf"      	[105,] "n005.F.MtvsXf"         
 # [67,] "lods.F.MtvsXo"         	 [80,] "coeff.F.MtvsXo"        	 [93,] "p.value.F.MtvsXo"      	[106,] "n005.F.MtvsXo"         
 # [68,] "lods.F.XfvsXo"         	 [81,] "coeff.F.XfvsXo"        	 [94,] "p.value.F.XfvsXo"      	[107,] "n005.F.XfvsXo"         
 # [69,] "lods.F.XfXovsAlMt"     	 [82,] "coeff.F.XfXovsAlMt"    	 [95,] "p.value.F.XfXovsAlMt"  	[108,] "n005.F.XfXovsAlMt"     
 # [70,] "lods.F.XfAlvsXoMt"     	 [83,] "coeff.F.XfAlvsXoMt"    	 [96,] "p.value.F.XfAlvsXoMt"  	[109,] "n005.F.XfAlvsXoMt"     
 # [71,] "lods.F.AlvsAll"        	 [84,] "coeff.F.AlvsAll"       	 [97,] "p.value.F.AlvsAll"     	[110,] "n005.F.AlvsAll"        
 # [72,] "lods.F.XfvsAll"        	 [85,] "coeff.F.XfvsAll"       	 [98,] "p.value.F.XfvsAll"     	[111,] "n005.F.XfvsAll"        
 # [73,] "lods.F.MtvsAll"        	 [86,] "coeff.F.MtvsAll"       	 [99,] "p.value.F.MtvsAll"     	[112,] "n005.F.MtvsAll"        
 # [74,] "lods.F.XovsAll"        	 [87,] "coeff.F.XovsAll"       	[100,] "p.value.F.XovsAll"     	[113,] "n005.F.XovsAll"        
 # [75,] "lods.F.XfXoDiffAlMt"   	 [88,] "coeff.F.XfXoDiffAlMt"  	[101,] "p.value.F.XfXoDiffAlMt"	[114,] "n005.F.XfXoDiffAlMt" 



#######  XY SCATTER PLOTS OF FEMALE AND MALE MONGAMOUS AND POLYGAMOUS SPECIES
#### FEMALE interactions plot with each lineage on axis
plot (dataB$coeff.F.AlvsMt, dataB$coeff.F.XfvsXo, xlab = "Mating Strategy Log Fold Bias (non-X-lineage)", ylab = "Mating Strategy Log Fold Bias (X-lineage)", cex.lab = 0.8, main= "Female Analysis", cex.main = 0.8, pch = 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main =2, col = nscol, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), )
abline(0,0)
abline(v=0, untf = FALSE)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.XfvsXo==-1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.XfvsXo==-1)], cex = 1.8, lwd=2, pch = 4, col= Xocol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.XfvsXo==1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.XfvsXo==1)], cex = 1.8, lwd=2, pch = 4, col= Xfcol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.AlvsMt==-1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.AlvsMt==-1)], cex = 1.5, pch = 16, col= Mtcol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.AlvsMt==1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.AlvsMt==1)], cex = 1.5, pch = 16, col= Alcol)
legend(x=-1.6 ,y=1.6 ,   pch = c(16,16,4,4), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al","Mt", "Xf","Xo"), col = c(Alcol,Mtcol,Xfcol,Xocol), cex = 1.5)
text(x=c(-0.2,0.2,0.2,-0.2), cex = 1.5,y=c(1.5,1.5,-1.5,-1.5), labels = c("8","7","7","2"))

#### MALE interactions plot with each lineage on axis
plot (dataB$coeff.M.AlvsMt, dataB$coeff.M.XfvsXo, xlab = "Mating Strategy Log Fold Bias (non-X-lineage)", ylab = "Mating Strategy Log Fold Bias (X-lineage)", cex.lab = 0.8, main= "Male Analysis", cex.main = 0.8, pch = 19, cex = 1, cex.axis = 1.5, cex.main = 2,cex.lab = 1.4, col = nscol, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
abline(0,0)
abline(v=0, untf = FALSE)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==-1)], cex = 1.8, lwd=2, pch = 4, col= Xocol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==1)], cex = 1.8, lwd=2, pch = 4, col= Xfcol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt==-1)], cex = 1.5, pch = 16, col= Mtcol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt==1)], cex = 1.5, pch = 16, col= Alcol)
legend(x=-1.6 ,y=1.6 ,   pch = c(16,16,4,4), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al","Mt", "Xf","Xo"), col = c(Alcol,Mtcol,Xfcol,Xocol), cex = 1.5)
text(x=c(-0.2,0.2,0.2,-0.2), cex = 1.5,y=c(1.5,1.5,-1.5,-1.5), labels = c("2","7","8","6"))

########## 20170706 reviewers asked for bias, even if not significant between male and famle
## might work to show these scatter mapped between sex
## These figurews are so hard to look at, they will not be included in supplementary material... 
#   but if anyone ever runs this script and wants me to explain them.... I'm happy too renns@reed.edu


# Male interaction plot w female expression mapped 
plot (dataB$coeff.M.AlvsMt, dataB$coeff.M.XfvsXo, xlab = "Mating Strategy Log Fold Bias (non-X-lineage)", ylab = "Mating Strategy Log Fold Bias (X-lineage)", cex.lab = 0.8, main= "Plot Male Analysis + female mapped",  pch = 19, cex = 1, cex.axis = 1.5, cex.main = 1.2,cex.lab = 1.4, col = nscol, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
abline(0,0)
abline(v=0, untf = FALSE)

points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==-1)], cex = 1.8, lwd=2, pch = 4, col= Xocol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==1)], cex = 1.8, lwd=2, pch = 4, col= Xfcol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt==-1)], cex = 1.5, pch = 16, col= Mtcol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt==1)], cex = 1.5, pch = 16, col= Alcol)

points (dataB$coeff.M.AlvsMt[which(dataB$n005.F.XfvsXo==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.F.XfvsXo==-1)], cex = 0.8, lwd=2, pch = 4, col= "pink")
points (dataB$coeff.M.AlvsMt[which(dataB$n005.F.XfvsXo==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.F.XfvsXo==1)], cex = 0.8, lwd=2, pch = 4, col= "red")
points (dataB$coeff.M.AlvsMt[which(dataB$n005.F.AlvsMt==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.F.AlvsMt==-1)], cex = 0.5, pch = 16, col= "light blue")
points (dataB$coeff.M.AlvsMt[which(dataB$n005.F.AlvsMt==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.F.AlvsMt==1)], cex = 0.5, pch = 16, col= "yellow")
legend(x=-1.6 ,y=1.6 ,   pch = c(16,16,4,4), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al male","Mt male", "Xf male","Xo male"), col = c(Alcol,Mtcol,Xfcol,Xocol), cex = 1)
legend(x=-1.6 ,y=-1.0 ,   pch = c(16,16,4,4), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al female","Mt female", "Xf female","Xo female"), col = c("yellow","light blue","red","pink"), cex = 1)
####### this is very ugly but it works


#### FEMALE interactions plot with male expression mapped  ### just reverse of above
plot (dataB$coeff.F.AlvsMt, dataB$coeff.F.XfvsXo, xlab = "Mating Strategy Log Fold Bias (non-X-lineage)", ylab = "Mating Strategy Log Fold Bias (X-lineage)", cex.lab = 0.8, main= "Female Analysis plotted + male mapped",  pch = 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main =1.2, col = nscol, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), )
abline(0,0)
abline(v=0, untf = FALSE)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.XfvsXo==-1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.XfvsXo==-1)], cex = 1.8, lwd=2, pch = 4, col= Xocol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.XfvsXo==1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.XfvsXo==1)], cex = 1.8, lwd=2, pch = 4, col= Xfcol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.AlvsMt==-1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.AlvsMt==-1)], cex = 1.5, pch = 16, col= Mtcol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.AlvsMt==1)], dataB$coeff.F.XfvsXo[which(dataB$n005.F.AlvsMt==1)], cex = 1.5, pch = 16, col= Alcol)

points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.XfvsXo==-1)], dataB$coeff.F.XfvsXo[which(dataB$n005.M.XfvsXo==-1)], cex = 0.5, lwd=1, pch = 4, col= "pink")
points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.XfvsXo==1)], dataB$coeff.F.XfvsXo[which(dataB$n005.M.XfvsXo==1)], cex = 0.5, lwd=1, pch = 4, col= "red")
points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.AlvsMt==-1)], dataB$coeff.F.XfvsXo[which(dataB$n005.M.AlvsMt==-1)], cex = 0.5, pch = 16, col= "light blue")
points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.AlvsMt==1)], dataB$coeff.F.XfvsXo[which(dataB$n005.M.AlvsMt==1)], cex = 0.5, pch = 16, col= "yellow")

legend(x=-1.6 ,y=-1 ,   pch = c(16,16,4,4), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al female","Mt female", "Xf female","Xo female"), col = c(Alcol,Mtcol,Xfcol,Xocol), cex = 1)
legend(x=-1.6 ,y=1.6 ,   pch = c(16,16,4,4), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al male","Mt male", "Xf male","Xo male"), col = c("yellow","light blue","red","pink"), cex = 1)





#### VOLCANO PLOTS MATING STYLE :::::: FIGURE 4
plot ( dataB$coeff.F.XfAlvsXoMt  , dataB$lods.F.XfAlvsXoMt , pch= 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-7,10), xlim=c(-1.5,1.5), col=nscol, xlab = "Log Fold Expression Difference (Female)", ylab = "Log Odds" )
# colored for female monogamy/polygamy (A)
points(dataB$coeff.F.XfAlvsXoMt[which(dataB$n005.F.XfAlvsXoMt==-1)],  dataB$lods.F.XfAlvsXoMt[dataB$n005.F.XfAlvsXoMt==-1], cex = 1.5, pch = 16, col = Mtcol)
points(dataB$coeff.F.XfAlvsXoMt[dataB$n005.F.XfAlvsXoMt==1],  dataB$lods.F.XfAlvsXoMt[which(dataB$n005.F.XfAlvsXoMt==1)], cex = 1.5, pch = 16, col = Alcol)
legend(x=0 ,y=10.8 ,  pch = c(16,16), bty="n",legend = c("female monogamy", "female polygyny"), col = c(Alcol,Mtcol), cex = 1.5) 
 abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 3)
##### plot female mating strategy volcanoe with male color (B)
plot ( dataB$coeff.F.XfAlvsXoMt  , dataB$lods.F.XfAlvsXoMt , pch= 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-7,10), xlim=c(-1.5,1.5), col=nscol, xlab = "Log Fold Expression Difference (Female)", ylab = "Log Odds" )
# colored for male monogamy/polygamy
points(dataB$coeff.F.XfAlvsXoMt[which(dataB$n005.M.XfAlvsXoMt==-1)],  dataB$lods.F.XfAlvsXoMt[dataB$n005.M.XfAlvsXoMt==-1], cex = 1.5, pch = 16, col = Mtcol)
points(dataB$coeff.F.XfAlvsXoMt[dataB$n005.M.XfAlvsXoMt==1],  dataB$lods.F.XfAlvsXoMt[which(dataB$n005.M.XfAlvsXoMt==1)], cex = 1.5, pch = 16, col = Alcol)
legend(x=0 ,y=10.8 ,  pch = c(16,16), bty="n",legend = c("male monogamy", "male polygyny"), col = c(Alcol,Mtcol), cex = 1.5)
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 3)
##### plot male mating strategyvolcanoe (C)
plot ( dataB$coeff.M.XfAlvsXoMt  , dataB$lods.M.XfAlvsXoMt ,pch= 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-7,10), xlim=c(-1.5,1.5), col=nscol,xlab = "Log Fold Expression Difference (Male)", ylab = "Log Odds"  )
# colored for male monogamy/polygamy
points(dataB$coeff.M.XfAlvsXoMt[which(dataB$n005.M.XfAlvsXoMt==-1)],  dataB$lods.M.XfAlvsXoMt[dataB$n005.M.XfAlvsXoMt==-1], cex = 1.5, pch = 16, col = Mtcol)
points(dataB$coeff.M.XfAlvsXoMt[dataB$n005.M.XfAlvsXoMt==1],  dataB$lods.M.XfAlvsXoMt[which(dataB$n005.M.XfAlvsXoMt==1)], cex = 1.5, pch = 16, col = Alcol)
legend(x=0 ,y=10.8 ,  pch = c(16,16), bty="n",legend = c("male monogamy", "male polygyny"), col = c(Alcol,Mtcol), cex = 1.5)
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 3)
##### plot male mating strategy volcanoe with female color (D)
plot ( dataB$coeff.M.XfAlvsXoMt  , dataB$lods.M.XfAlvsXoMt , pch= 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-7,10), xlim=c(-1.5,1.5), col=nscol, xlab = "Log Fold Expression Difference (Male)", ylab = "Log Odds"  )
# colored for female monogamy/polygamy
points(dataB$coeff.M.XfAlvsXoMt[which(dataB$n005.F.XfAlvsXoMt==-1)],  dataB$lods.M.XfAlvsXoMt[dataB$n005.F.XfAlvsXoMt==-1], cex = 1.5, pch = 16, col = Mtcol)
points(dataB$coeff.M.XfAlvsXoMt[dataB$n005.F.XfAlvsXoMt==1],  dataB$lods.M.XfAlvsXoMt[which(dataB$n005.F.XfAlvsXoMt==1)], cex = 1.5, pch = 16, col = Alcol)
legend(x=0 ,y=10.8 ,  pch = c(16,16), bty="n",legend = c("female monogamy", "female polygyny"), col = c(Alcol,Mtcol), cex = 1.5)
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 3)
  
  
  ##################### 2017 a couple extra figures
  ######################
##
##### plot pairwise species volcanoe with mating strategy 
### male mating strategy
plot ( dataB$coeff.M.XfAlvsXoMt  , dataB$lods.M.XfAlvsXoMt , pch= 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-7,10), xlim=c(-1.5,1.5), col=nscol, xlab = "Log Fold Expression Difference (Male)", ylab = "Log Odds"  )
# colored for Al/Mt male
points(dataB$coeff.M.XfAlvsXoMt[which(dataB$n005.M.AlvsMt==-1)],  dataB$lods.M.XfAlvsXoMt[which(dataB$n005.M.AlvsMt==-1)], cex = 0.8, pch = 16, col = Mtcol)
points(dataB$coeff.M.XfAlvsXoMt[which(dataB$n005.M.AlvsMt==1)],  dataB$lods.M.XfAlvsXoMt[which(dataB$n005.M.AlvsMt==1)], cex = 0.8, pch = 16, col = Alcol)
legend(x=0 ,y=10.8 ,  pch = c(16,16), bty="n",legend = c("Al monogamy", "Mt polygamy"), col = c(Alcol,Mtcol), cex = 1.5)
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 3)
#### male mating strategy
plot ( dataB$coeff.M.XfAlvsXoMt  , dataB$lods.M.XfAlvsXoMt , pch= 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-7,10), xlim=c(-1.5,1.5), col=nscol, xlab = "Log Fold Expression Difference (Male)", ylab = "Log Odds"  )
# colored for Al/Mt male
points(dataB$coeff.M.XfAlvsXoMt[which(dataB$n005.M.XfvsXo==-1)],  dataB$lods.M.XfAlvsXoMt[dataB$n005.M.XfvsXo==-1], cex = 1.5, pch = 12, col = Xocol)
points(dataB$coeff.M.XfAlvsXoMt[dataB$n005.M.XfvsXo==1],  dataB$lods.M.XfAlvsXoMt[which(dataB$n005.M.XfvsXo==1)], cex = 1.5, pch = 12, col = Xfcol)
legend(x=0 ,y=8.8 ,  pch = c(12,12), bty="n",legend = c("Xf monogamy", "Xo polygamy"), col = c(Alcol,Mtcol), cex = 1.5)
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 3)
####### these are great, there is so much here
####### try the other way
# plot species comparison but color for mating system
plot ( dataB$coeff.M.AlvsMt  , dataB$lods.M.AlvsMt , pch= 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, ylim = c(-7,10), xlim=c(-1.5,1.5), col=nscol, xlab = "Log Fold Expression Difference (Male) nonX", ylab = "Log Odds"  )
# colored for Al/Mt male
points(dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfAlvsXoMt==-1)],  dataB$lods.M.AlvsMt[dataB$n005.M.XfAlvsXoMt==-1], cex = 1.5, pch = 16, col = Mtcol)
points(dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfAlvsXoMt==1)],  dataB$lods.M.AlvsMt[which(dataB$n005.M.XfAlvsXoMt==1)], cex = 1.5, pch = 16, col = Alcol)
legend(x=0 ,y=10.8 ,  pch = c(16,16), bty="n",legend = c("Xf.Al monogamy", "Xo.Mt polygamy"), col = c(Alcol,Mtcol), cex = 1.5)
  abline(-2.298853076409706,0, col =nscol, lty = 3, lwd = 3)

#####################
#######  XY SCATTER PLOTS OF XfXovsALMT MONGAMOUS AND POLYGAMOUS between sexes 
#### NonX lineage across sex
plot (dataB$coeff.F.AlvsMt, dataB$coeff.M.AlvsMt, xlab = "Mating Strategy Log Fold Bias (Female)", ylab = "Mating Strategy Log Fold Bias (Male)", cex.lab = 0.8, main= "NonX-lineage", cex.main = 0.8, pch = 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main =2, col = nscol, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), )
abline(0,0)
abline(v=0, untf = FALSE)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.AlvsMt==-1)], dataB$coeff.M.AlvsMt[which(dataB$n005.F.AlvsMt==-1)], cex = 0.8, lwd=2, pch = 16, col= Mtcol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.F.AlvsMt==1)], dataB$coeff.M.AlvsMt[which(dataB$n005.F.AlvsMt==1)], cex = 0.8, lwd=2, pch = 16, col= Alcol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.AlvsMt==-1)], dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==-1)], cex = 1.5, pch = 24, col= Mtcol)
points (dataB$coeff.F.AlvsMt[which(dataB$n005.M.AlvsMt==1)], dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==1)], cex = 1.5, pch = 24, col= Alcol)
legend(x=-1.6 ,y=1.6 ,   pch = c(16,16,24,24), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al.F","Mt.F", "Al.M","Mt.M"), col = c(Alcol,Mtcol,Alcol,Mtcol), cex = 1.5)
text(x=c(-0.2,0.2,0.2,-0.2), cex = 1.5,y=c(1.5,1.5,-1.5,-1.5), labels = c("0","81","0","28"))
length(which(dataB$n005.F.AlvsMt==-1 & dataB$n005.M.AlvsMt==-1)) # 28
length(which(dataB$n005.F.AlvsMt==1 & dataB$n005.M.AlvsMt==1)) # 81

#### X lineage across sex
plot (dataB$coeff.F.XfvsXo, dataB$coeff.M.XfvsXo, xlab = "Mating Strategy Log Fold Bias (Female)", ylab = "Mating Strategy Log Fold Bias (Male)", cex.lab = 0.8, main= "X-lineage", cex.main = 0.8, pch = 16, cex = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main =2, col = nscol, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), )
abline(0,0)
abline(v=0, untf = FALSE)
points (dataB$coeff.F.XfvsXo[which(dataB$n005.F.XfvsXo==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.F.XfvsXo==-1)], cex = 0.8, lwd=2, pch = 16, col= Xocol)
points (dataB$coeff.F.XfvsXo[which(dataB$n005.F.XfvsXo==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.F.XfvsXo==1)], cex = 0.8, lwd=2, pch = 16, col= Xfcol)
points (dataB$coeff.F.XfvsXo[which(dataB$n005.M.XfvsXo==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==-1)], cex = 1.5, pch = 24, col= Xocol)
points (dataB$coeff.F.XfvsXo[which(dataB$n005.M.XfvsXo==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==1)], cex = 1.5, pch = 24, col= Xfcol)
legend(x=-1.6 ,y=1.6 ,   pch = c(16,16,24,24), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Xo.F","Xf.F", "Xo.M","Xf.M"), col = c(Xfcol,Xocol,Xfcol,Xocol), cex = 1.5)
text(x=c(-0.2,0.2,0.2,-0.2), cex = 1.5,y=c(1.5,1.5,-1.5,-1.5), labels = c("2","1","6","7"))
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.M.XfvsXo==-1)) # 7
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.M.XfvsXo==1)) # 1
length(which(dataB$n005.F.XfvsXo==1 & dataB$n005.M.XfvsXo==-1)) # 6
length(which(dataB$n005.F.XfvsXo==-1 & dataB$n005.M.XfvsXo==1)) # 2



#### MALE interactions plot with each lineage on axis
plot (dataB$coeff.M.AlvsMt, dataB$coeff.M.XfvsXo, xlab = "Mating Strategy Log Fold Bias (non-X-lineage)", ylab = "Mating Strategy Log Fold Bias (X-lineage)", cex.lab = 0.8, main= "Male Analysis", cex.main = 0.8, pch = 19, cex = 1, cex.axis = 1.5, cex.main = 2,cex.lab = 1.4, col = nscol, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
abline(0,0)
abline(v=0, untf = FALSE)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==-1)], cex = 1.8, lwd=2, pch = 4, col= Xocol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.XfvsXo==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.XfvsXo==1)], cex = 1.8, lwd=2, pch = 4, col= Xfcol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==-1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt==-1)], cex = 1.5, pch = 16, col= Mtcol)
points (dataB$coeff.M.AlvsMt[which(dataB$n005.M.AlvsMt==1)], dataB$coeff.M.XfvsXo[which(dataB$n005.M.AlvsMt==1)], cex = 1.5, pch = 16, col= Alcol)
legend(x=-1.6 ,y=1.6 ,   pch = c(16,16,4,4), bty="n", pt.lwd = c(0,0,2,2) ,legend = c("Al","Mt", "Xf","Xo"), col = c(Alcol,Mtcol,Xfcol,Xocol), cex = 1.5)
text(x=c(-0.2,0.2,0.2,-0.2), cex = 1.5,y=c(1.5,1.5,-1.5,-1.5), labels = c("2","7","8","6"))
  
   
  
  

#########################################################
## FILENAME:: GenomicVolcano_ExpressionCompareFigures.R  #### MOVED TO THIS CODE
## GOAL:   on 20130725_GenomicVolcano_ExpressionCompareREcheck.R I convinced myself all was OK 
#	(20140205 found error in code but it is still OK corrected to use all gDNA5232 which was in different order than dataB)
##        remade the figures using the 005 no FDR faint filtered results
##         and plotted only those that passed in both  analyses
##         ploted only relevant for paper but checked correlation on all
## 

gDNA<- read.table(file="21031010_ResultsEctodiniGenFlagFilterLODScollapsed.txt", header = TRUE)
dim(gDNA) # [1] 17118    39
dim(dataB) # [1] 5232  114
colnames(gDNA)
gDNA$ID # has hhID and TC for full array
# merge with dataB by $ID
dataB$ID

gDNA5232 <- merge(dataB, gDNA, by = "ID")
dim(gDNA5232) # 5232  152
gDNA5232$ID[28]
dataB$ID[28]  ## they are not in the same order, but gDNA5232 has all columns


length(which(gDNA5232$coef.g.AlvsMt > 0 & (gDNA5232$n005.F.AlvsMt==1))) # E:AL G:AL 195
length(which(gDNA5232$coef.g.AlvsMt < 0 & (gDNA5232$n005.F.AlvsMt==1))) # E:AL G:Mt 203
length(which(gDNA5232$coef.g.AlvsMt > 0 & (gDNA5232$n005.F.AlvsMt==-1))) # E:Mt G:AL 98
length(which(gDNA5232$coef.g.AlvsMt < 0 & (gDNA5232$n005.F.AlvsMt==-1))) # E:Mt G:Mt 78

length(which(gDNA5232$coef.g.AlvsXf > 0 & (gDNA5232$n005.F.AlvsXf==1))) # E:Al G:Al 30
length(which(gDNA5232$coef.g.AlvsXf < 0 & (gDNA5232$n005.F.AlvsXf==1))) # E:Al G:Xf 22
length(which(gDNA5232$coef.g.AlvsXf > 0 & (gDNA5232$n005.F.AlvsXf==-1))) # E:Xf G:Al 34
length(which(gDNA5232$coef.g.AlvsXf < 0 & (gDNA5232$n005.F.AlvsXf==-1))) # E:Xf G:Xf 34

length(which(gDNA5232$coef.g.AlvsXo > 0 & (gDNA5232$n005.F.AlvsXo==1))) # E:Al G:Al 38
length(which(gDNA5232$coef.g.AlvsXo < 0 & (gDNA5232$n005.F.AlvsXo==1))) # E:Al G:X0 47
length(which(gDNA5232$coef.g.AlvsXo > 0 & (gDNA5232$n005.F.AlvsXo==-1))) # E:Xo G:Al 35
length(which(gDNA5232$coef.g.AlvsXo < 0 & (gDNA5232$n005.F.AlvsXo==-1))) # E:Xo G:Xo 33

length(which(gDNA5232$coef.g.MtvsXf > 0 & (gDNA5232$n005.F.MtvsXf==1))) # E:Mt G:Mt 90
length(which(gDNA5232$coef.g.MtvsXf < 0 & (gDNA5232$n005.F.MtvsXf==1))) # E:Mt G:Xf 76
length(which(gDNA5232$coef.g.MtvsXf > 0 & (gDNA5232$n005.F.MtvsXf==-1))) # E:Xf G:Mt 185
length(which(gDNA5232$coef.g.MtvsXf < 0 & (gDNA5232$n005.F.MtvsXf==-1))) # E:Xf G:Xf 259

length(which(gDNA5232$coef.g.MtvsXo > 0 & (gDNA5232$n005.F.MtvsXo==1))) # E:Mt G:Mt 125
length(which(gDNA5232$coef.g.MtvsXo < 0 & (gDNA5232$n005.F.MtvsXo==1))) # E:Mt G:Xo 159
length(which(gDNA5232$coef.g.MtvsXo > 0 & (gDNA5232$n005.F.MtvsXo==-1))) # E:Xo G:Mt 237
length(which(gDNA5232$coef.g.MtvsXo < 0 & (gDNA5232$n005.F.MtvsXo==-1))) # E:Xo G:Xo 156

length(which(gDNA5232$coef.g.XfvsXo > 0 & (gDNA5232$n005.F.XfvsXo==1))) # E:Xf G:Xf 38
length(which(gDNA5232$coef.g.XfvsXo < 0 & (gDNA5232$n005.F.XfvsXo==1))) # E:Xf G:Xo 30
length(which(gDNA5232$coef.g.XfvsXo > 0 & (gDNA5232$n005.F.XfvsXo==-1))) # E:Xo G:Xf 12
length(which(gDNA5232$coef.g.XfvsXo < 0 & (gDNA5232$n005.F.XfvsXo==-1))) # E:Xo G:Xo 20


######### PLOTTING GENOMIC VOLCANOES  ::::: FIGURE 2
# AlMT (colored for Female expression) (B2)
plot(gDNA5232$coef.g.AlvsMt, gDNA5232$lods.g.AlvsMt, xlim = c(-2,2), ylim = c(-6,4.5), xlab = "Genomic Hybridization Ratio", ylab = "Log Odds", col=nscol, pch=19, cex = 1, cex.lab = 1.5, cex.axis = 1.5)
points(gDNA5232$coef.g.AlvsMt[which(gDNA5232$n005.F.AlvsMt==1)], gDNA5232$lods.g.AlvsMt[which(gDNA5232$n005.F.AlvsMt==1)], col= Alcol, pch = 19, cex = 1.5 )
points(gDNA5232$coef.g.AlvsMt[which(gDNA5232$n005.F.AlvsMt==-1)], gDNA5232$lods.g.AlvsMt[which(gDNA5232$n005.F.AlvsMt==-1)], col= Mtcol, pch = 19, cex = 1.5 )
legend(x="topleft", legend=c("Al biased expression in females","Mt biased expression in females"),col=c(Alcol,Mtcol), pch=c(19, 19), bty="n", cex = 1.5)
# AlMt colored for ,ale expression (B1)
plot(gDNA5232$coef.g.AlvsMt, gDNA5232$lods.g.AlvsMt, xlim = c(-2,2),  ylim = c(-6,4.5), xlab = "Genomic Hybridization Ratio", ylab = "Log Odds", col=nscol, pch=19, cex = 1, cex.lab = 1.5, cex.axis = 1.5)
points(gDNA5232$coef.g.AlvsMt[which(gDNA5232$n005.M.AlvsMt==1)], gDNA5232$lods.g.AlvsMt[which(gDNA5232$n005.M.AlvsMt==1)], col= Alcol, pch = 19, cex = 1.5 )
points(gDNA5232$coef.g.AlvsMt[which(gDNA5232$n005.M.AlvsMt==-1)], gDNA5232$lods.g.AlvsMt[which(gDNA5232$n005.M.AlvsMt==-1)], col= Mtcol, pch = 19, cex = 1.5 )
legend(x="topleft", legend=c("Al biased expression in males","Mt biased expression in males"),col=c(Alcol,Mtcol), pch=c(19, 19), bty="n", cex = 1.5)
# XfXo (colored for Female expression)(A2)
plot(gDNA5232$coef.g.XfvsXo, gDNA5232$lods.g.XfvsXo, xlim = c(-2,2),  ylim = c(-6,4.5),  xlab = "Genomic Hybridization Ratio", ylab = "Log Odds", col=nscol, pch=19, cex = 1, cex.lab = 1.5, cex.axis = 1.5)
points(gDNA5232$coef.g.XfvsXo[which(gDNA5232$n005.F.XfvsXo==1)], gDNA5232$lods.g.XfvsXo[which(gDNA5232$n005.F.XfvsXo==1)], col= Xfcol, pch = 19, cex = 1.5 )
points(gDNA5232$coef.g.XfvsXo[which(gDNA5232$n005.F.XfvsXo==-1)], gDNA5232$lods.g.XfvsXo[which(gDNA5232$n005.F.XfvsXo==-1)], col= Xocol, pch = 19, cex = 1.5 )
legend(x="topleft", legend=c("Xf biased expression in females","Xo biased expression in females"),col=c(Xfcol,Xocol), pch=c(19, 19), bty="n", cex = 1.5)
# XfXo colored for male expression (A1)
plot(gDNA5232$coef.g.XfvsXo, gDNA5232$lods.g.XfvsXo, xlim = c(-2,2),  ylim = c(-6,4.5), xlab = "Genomic Hybridization Ratio", ylab = "Log Odds", col=nscol, pch=19, cex = 1, cex.lab = 1.5, cex.axis = 1.5)
points(gDNA5232$coef.g.XfvsXo[which(gDNA5232$n005.M.XfvsXo==1)], gDNA5232$lods.g.XfvsXo[which(gDNA5232$n005.M.XfvsXo==1)], col= Xfcol, pch = 19, cex = 1.5 )
points(gDNA5232$coef.g.XfvsXo[which(gDNA5232$n005.M.XfvsXo==-1)], gDNA5232$lods.g.XfvsXo[which(gDNA5232$n005.M.XfvsXo==-1)], col= Xocol, pch = 19, cex = 1.5 )
legend(x="topleft", legend=c("Xf biased expression in males","Xo biased expression in males"),col=c(Xfcol,Xocol), pch=c(19, 19), bty="n", cex = 1.5)

########
### 20140205 supplementary table 1 correlation data
###
?cor.test
cor.test(gDNA5232$coef.g.AlvsMt[which(gDNA5232$n005.M.AlvsMt!=0)], dataB$coeff.M.AlvsMt[which(gDNA5232$n005.M.AlvsMt!=0)], method = "pearson")  # 0.0108920 df = 197, p-value = 0.6168 ******
cor.test(gDNA5232$coef.g.AlvsMt[which(gDNA5232$n005.F.AlvsMt!=0)], dataB$coeff.F.AlvsMt[which(gDNA5232$n005.F.AlvsMt!=0)], method = "pearson")  # 0.02407462  df = 572, p-value = 0.5649
cor.test(gDNA5232$coef.g.AlvsXf[which(gDNA5232$n005.M.AlvsXf!=0)], dataB$coeff.M.AlvsXf[which(gDNA5232$n005.M.AlvsXf!=0)], method = "pearson") # -0.07062709  df = 123, p-value = 0.4338
cor.test(gDNA5232$coef.g.AlvsXf[which(gDNA5232$n005.F.AlvsXf!=0)], dataB$coeff.F.AlvsXf[which(gDNA5232$n005.F.AlvsXf!=0)], method = "pearson") # 0.1189866   df = 118, p-value = 0.1955
cor.test(gDNA5232$coef.g.AlvsXo[which(gDNA5232$n005.M.AlvsXo!=0)], dataB$coeff.M.AlvsXo[which(gDNA5232$n005.M.AlvsXo!=0)], method = "pearson") # -0.08211452  df = 32, p-value = 0.6443
cor.test(gDNA5232$coef.g.AlvsXo[which(gDNA5232$n005.F.AlvsXo!=0)], dataB$coeff.F.AlvsXo[which(gDNA5232$n005.F.AlvsXo!=0)], method = "pearson") # -0.06430023  df = 151, p-value = 0.4297
cor.test(gDNA5232$coef.g.MtvsXf[which(gDNA5232$n005.M.MtvsXf!=0)], dataB$coeff.M.MtvsXf[which(gDNA5232$n005.M.MtvsXf!=0)], method = "pearson") # 0.1338658 df = 250, p-value = 0.03367 ********
cor.test(gDNA5232$coef.g.MtvsXf[which(gDNA5232$n005.F.MtvsXf!=0)], dataB$coeff.F.MtvsXf[which(gDNA5232$n005.F.MtvsXf!=0)], method = "pearson") # 0.09415609 df = 608, p-value = 0.02002  OOOPS SIGNIFICANT  BUT NOT ONE WE CARE MUCH ABOUT
cor.test(gDNA5232$coef.g.MtvsXo[which(gDNA5232$n005.M.MtvsXo!=0)], dataB$coeff.M.MtvsXo[which(gDNA5232$n005.M.MtvsXo!=0)], method = "pearson") # 0.01083949 df = 161, p-value = 0.8908 *********
cor.test(gDNA5232$coef.g.MtvsXo[which(gDNA5232$n005.F.MtvsXo!=0)], dataB$coeff.F.MtvsXo[which(gDNA5232$n005.F.MtvsXo!=0)], method = "pearson") # -0.05160101 df = 675, p-value = 0.1799
cor.test(gDNA5232$coef.g.XfvsXo[which(gDNA5232$n005.M.XfvsXo!=0)], dataB$coeff.M.XfvsXo[which(gDNA5232$n005.M.XfvsXo!=0)], method = "pearson") # 0.197019 df = 85, p-value = 0.06739
cor.test(gDNA5232$coef.g.XfvsXo[which(gDNA5232$n005.F.XfvsXo!=0)], dataB$coeff.F.XfvsXo[which(gDNA5232$n005.F.XfvsXo!=0)], method = "pearson") # -0.04015056 df = 98, p-value = 0.6916
###### total number regulated by all 12 pairwise comparisons (6 male 6 female)
cbind(colnames(gDNA5232))  ## F pairwise = 102-107  M pairwise = 50 - 55
### fuck me this is dumb but I don't care
t1<- which(gDNA5232[,50]!=0)
t2<- which(gDNA5232[,51]!=0)
t3<- which(gDNA5232[,52]!=0)
t4<- which(gDNA5232[,53]!=0)
t5<- which(gDNA5232[,54]!=0)
t6<- which(gDNA5232[,55]!=0)
t7<- which(gDNA5232[,102]!=0)
t8<- which(gDNA5232[,103]!=0)
t9<- which(gDNA5232[,104]!=0)
t10<- which(gDNA5232[,105]!=0)
t11<- which(gDNA5232[,106]!=0)
t12<- which(gDNA5232[,107]!=0)
length(unique(c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12))) # 1284

# ### NOT USED IN MANUSCRIPT
# ################################
# ## PCA Analysis of individual M and F analyses
# ## get just those that passed M and F for analysis (i.e. dataB)
# F4PCA <- merge(dataB, writeFind.table, by = "ID") # coeffs are 181 - 199
# M4PCA <- merge(dataB, writeMind.table, by = "ID")
# colnames(F4PCA) # coeffs are 181 - 200
# # pull those out, add a column of 0s for reference make matrix  m= Al103 ; f=Xo116f

# # Plot female
# FPCAcolors <- c(rep(Alcol,5), rep(Mtcol,5), rep(Xfcol,5), rep(Xocol,4))
# Xo116f  <- rep(0,5232)
# FPCAlabelsRef <- c(colnames(F4PCA[124:142]),"Xo116f" ) 
# #### inefficient code here cuz I didn't know how to center = TRU
# F4PCAmRef <- as.matrix(cbind(F4PCA[,c(181:199)], Xo116f))
# FPRCref <- prcomp(F4PCAmRef, center=TRUE )  #### try scale = 0
# plot(FPRCref$rotation[,1], FPRCref$rotation[,2],  col = c( FPCAcolors, Xocol), pch = 19, xlab="Female PC 1", ylab= "Female PC 2", main = "Female Analysis", cex = 2, ylim = c(0.42,-0.42), xlim = c(-0.02, 0.42))
# text(FPRCref$rotation[,1], FPRCref$rotation[,2], FPCAlabelsRef, pos=1, col = c(FPCAcolors, Xocol))

# # Plot male
# MPCAcolors <- c(rep(Alcol,4), rep(Mtcol,5), rep(Xfcol,5), rep(Xocol,5))
# Al103 <- rep(0,5232)
# MPCAlabelsRef <- c("Al103",colnames(M4PCA[124:142])) 
# #### inefficient code here cuz I didn't know how to center = TRU
# M4PCAmRef <- as.matrix(cbind(M4PCA[,c(181:199)], Al103))
# MPRCref <- prcomp(M4PCAmRef, center=TRUE )  #### try scale = 0
# plot(MPRCref$rotation[,1], MPRCref$rotation[,2],  col = c( Alcol, MPCAcolors), pch = 19, xlab="Male PC 1", ylab= "Male PC 2", main = "Male Analysis", cex = 2, ylim = c(-0.45,0.5), xlim = c(-0.32, 0.02))
# text(MPRCref$rotation[,1], MPRCref$rotation[,2], MPCAlabelsRef, pos=1, col = c(Alcol, MPCAcolors))


###
###
## getting the AVT clone plotted   ### NOT USED IN MANUSCRIPT
cbind(colnames(F4PCA))  ## need 181 - 200  (200-218 is the stdev)
#############hh_Ab_HarvardCol_000008099
which(F4PCA[,1]=="hh_Ab_HarvardCol_000008099") #1093  ## in house AVT
F4PCA[1093,181:199]
plot(c(1:20),c(F4PCA[1093,c(181:199)],0),  ylab = "Log2 Expression Ratio", col = c(FPCAcolors, Xocol) , pch =19,  xaxt="n" , xlab= "", main = "hh_Ab_HarvardCol_000008099 in Females", ylim = c(-0.65, 0.85), cex=2)
axis( side = 1, at = c(1:20), labels = c(colnames(F4PCA[124:142]),"Xo116f") , las = 2 )
abline(0,0)
#
which(M4PCA[,1]=="hh_Ab_HarvardCol_000008099") #1093
M4PCA[1093,181:199] # OK
plot(c(1:20),c(0,M4PCA[1093,c(181:199)]),  ylab = "Log2 Expression Ratio", col = c(Alcol,MPCAcolors ), pch =19,  xaxt="n" , xlab= "", main = "hh_Ab_HarvardCol_000008099 in Males", ylim = c(-0.65, 0.85), cex = 2)
axis( side = 1, at = c(1:20), labels = c("Al103",colnames(M4PCA[124:142])) , las = 2 )
abline(0,0)

############# hh_Ab_StanfordCol_000005706  # ab.gene.s26.40  ## AVT
which(F4PCA[,1]=="hh_Ab_StanfordCol_000005706") #4407
F4PCA[4407,c(181:199)]
plot(c(1:20),c(F4PCA[4407,c(181:199)],0),  ylab = "Log2 Expression Ratio", col = c(FPCAcolors, Xocol) , pch =19,  xaxt="n" , xlab= "", main = "hhh_Ab_StanfordCol_000005706 in Females", ylim = c(-0.65, 0.85), cex=2)
axis( side = 1, at = c(1:20), labels = c(colnames(F4PCA[124:142]),"Xo116f") , las = 2 )
abline(0,0)

plot(c(1:20),c(0,M4PCA[4407,c(181:199)]),  ylab = "Log2 Expression Ratio", col = c(Alcol,MPCAcolors ), pch =19,  xaxt="n" , xlab= "", main = "hh_Ab_StanfordCol_000005706 in males", ylim = c(-0.65, 0.85), cex = 2)
axis( side = 1, at = c(1:20), labels = c("Al103",colnames(M4PCA[124:142])) , las = 2 )
abline(0,0)


#TC1685 = another AVT (these are not compressed because they dodn't have overlapping sequences ?)
which(writeFind.table[,7]=="TC1685") #15938
cbind(colnames(writeFind.table))  ## 68:86 are coeffs (87-105 are stdev)
writeFind.table[15938, 68:86]
plot(c(1:20),c(writeFind.table[15938,c(68:86)],0),  ylab = "Log2 Expression Ratio", col = c(FPCAcolors, Xocol) , pch =19,  xaxt="n" , xlab= "", main = "hh_Ab_HarvardCol_000008099 in Females",  cex=2)
axis( side = 1, at = c(1:20), labels = c(colnames(writeFind.table[11:29]),"Xo116f") , las = 2 )
abline(0,0)

### lets look as noise, stdeve of the features in individual analysis
mean(as.numeric(F4PCA[1093, 200:218])) # 1.349038
mean(as.numeric(M4PCA[1093, 200:218])) # 1.191557
mean(as.numeric(F4PCA[4407, 200:218])) # 1.3208
mean(as.numeric(M4PCA[4407, 200:218])) # 1.736459
mean(as.numeric(writeFind.table[15938, 87:105])) # 1.552551
################################
########
######## Making GENE LISTS
################################

#############################################################################
### get gene annotations
annotR <- read.table (file = "fishprint4.03_AnnotationsBroadandInhouseU20130725.txt", header = TRUE, sep = "\t", quote = "")
head(annotR) # not sure why some have quotes, oh well...
annotR[which(annotR[,2]=="ab.gene.s26.40"),]
annotR[which(annotR[,1]=="hh_Ab_HarvardCol_000008099"), ]
annotR[which(annotR[,1]=="TC1685"), ]
#################################################################################
## merge Annotations (BROAD Plus Alternates) with results and write the table 
################# writing the table to play with in Excel
dim(annotR) # [1] 7915    5
annotatedR<- merge(annotR, dataB, by.x = "TCorIDannot", by.y =  "ID", all.y = TRUE)
dim(annotatedR) # [1] 5232  118
write.table(annotatedR,file="20140206_EctodiniMaskedResultsAllAnnotated_wMT104.tab",  row.names = FALSE, sep="\t")
head(annotatedR)
cbind(colnames(annotatedR))
  # [1,] "TCorIDannot"           
  # [2,] "BROADgeneID.x"    ### this includes "none" as well as NA change that       
  # [3,] "BROADgeneTreeInfo1"  ### this is ENSEMBLE number for annotation but includes NONE; ok  
  # [4,] "BROADgeneTreeName"     
  # [5,] "BROADgenePlusInhouse"  
  # [6,] "Block"                 
  # [7,] "Row"                   
  # [8,] "Column"                
  # [9,] "DFCI_TC_2009"          
 # [10,] "DFCI_TC_2011"          
 # [11,] "TCorID"                
 # [12,] "TCorGB"                
 # [13,] "TCorGBorHH"            
 # [14,] "BROADgeneID.y"         

 # [54,] "n005.M.AlvsMt"         
 # [55,] "n005.M.AlvsXf"         
 # [56,] "n005.M.AlvsXo"         
 # [57,] "n005.M.MtvsXf"         
 # [58,] "n005.M.MtvsXo"         
 # [59,] "n005.M.XfvsXo"         
 # [60,] "n005.M.XfXovsAlMt"     
 # [61,] "n005.M.XfAlvsXoMt"     
 # [62,] "n005.M.AlvsAll"        
 # [63,] "n005.M.XfvsAll"        
 # [64,] "n005.M.MtvsAll"        
 # [65,] "n005.M.XovsAll"        
 # [66,] "n005.M.XfXoDiffAlMt"   
 
# [106,] "n005.F.AlvsMt"         
# [107,] "n005.F.AlvsXf"         
# [108,] "n005.F.AlvsXo"         
# [109,] "n005.F.MtvsXf"         
# [110,] "n005.F.MtvsXo"         
# [111,] "n005.F.XfvsXo"         
# [112,] "n005.F.XfXovsAlMt"     
# [113,] "n005.F.XfAlvsXoMt"     
# [114,] "n005.F.AlvsAll"        
# [115,] "n005.F.XfvsAll"        
# [116,] "n005.F.MtvsAll"        
# [117,] "n005.F.XovsAll"        
# [118,] "n005.F.XfXoDiffAlMt"   
head(annotatedR[c(1,2,3)], 50)
annotatedR[which(annotatedR[,2]=="none"),2] <- NA
head(annotatedR[c(1,2,3)], 50) # looks good


##################################################################
# species specific (Using 4 columns
AlUp <- which(annotatedR$n005.M.AlvsAll ==1 & annotatedR$n005.F.AlvsAll ==1)
AlUp2<- cbind(AlUp, rep("Al-up", length(AlUp)))
AlDown <-which(annotatedR$n005.M.AlvsAll ==-1 & annotatedR$n005.F.AlvsAll ==-1)
AlDown2<- cbind(AlDown, rep("Al-down", length(AlDown)))
XfUp <- which(annotatedR$n005.M.XfvsAll ==1 & annotatedR$n005.F.XfvsAll ==1)
XfUp2<- cbind(XfUp, rep("Xf-up", length(XfUp)))
XfDown <- which(annotatedR$n005.M.XfvsAll ==-1 & annotatedR$n005.F.XfvsAll ==-1)
XfDown2<- cbind(XfDown, rep("Xf-down", length(XfDown)))
MtUp <- which(annotatedR$n005.M.MtvsAll ==1 & annotatedR$n005.F.MtvsAll ==1)
MtUp2<- cbind(MtUp, rep("Mt-up", length(MtUp)))
MtDown <- which(annotatedR$n005.M.MtvsAll ==-1 & annotatedR$n005.F.MtvsAll ==-1)
MtDown2<- cbind(MtDown, rep("Mt-down", length(MtDown)))
XoUp <- which(annotatedR$n005.M.XovsAll ==1 & annotatedR$n005.F.XovsAll ==1)
XoUp2<- cbind(XoUp, rep("Xo-up", length(XoUp)))
XoDown <- which(annotatedR$n005.M.XovsAll ==-1 & annotatedR$n005.F.XovsAll ==-1)
XoDown2<- cbind(XoDown, rep("Xo-down", length(XoDown)))

Al3 <- rbind(AlUp2, AlDown2)
Xf3 <- rbind(XfUp2, XfDown2)
Mt3 <- rbind(MtUp2, MtDown2)
Xo3 <- rbind(XoUp2, XoDown2)

temp1 <- merge(Al3,Mt3, by = 1, all = TRUE)
temp2 <- merge (Xf3, Xo3, by = 1, all = TRUE)
temp4 <- merge (temp1, temp2, by=1, all = TRUE)
temp5 <- temp4[ order(temp4[,5], temp4[,2], temp4[,4], temp4[,3] ),]
#### shit, there are a handful that are called down in more than one species
# 2634    <NA>   Mt-up Xf-down Xo-down
# 330    <NA>   Mt-up Xf-down Xo-down
# 481    <NA>   Mt-up Xf-down Xo-down
# 2442   Al-up Mt-down Xf-down    <NA>
# 719    <NA> Mt-down   Xf-up   Xo-up
##### looking at the data they really are, with weak P but they are.
##### when sorted they appear at the top of the list which will call attention
##### to me this says the threshold is not stringent enough to eliminate noise
colnames(temp5) <- c("row","Al.bias","Mt.bias","Xf.bias","Xo.bias")

SuppTableSpecies2 <- cbind(temp5[,c(2,3,4,5)], annotatedR[as.character(temp5[,1]),c(1,2,3,5)])
dim(SuppTableSpecies2) #221 (lost 2 with mt104 included)
sum(!is.na(SuppTableSpecies2[,6])) #135 (down 5 now that mt104 is in)
sum(!is.na(SuppTableSpecies2[,8])) #116 (down 6 now that mt104 is in)

write.table(SuppTableSpecies2,file="SuppTableSpecies2_20140205.tab",  row.names = FALSE, sep="\t")


##################################################################
# mating strategy (will need to merge based on row number somehow)
## make a male reg column a female reg column and a concordant/discordant column
mMo <- which(annotatedR$n005.M.XfAlvsXoMt ==1) #  male monogamy
mMo <- cbind(mMo, rep("mMo", length(mMo)))
mPo <- which(annotatedR$n005.M.XfAlvsXoMt ==-1) #  male polygamy
mPo <- cbind(mPo, rep("mPo", length(mPo)))
mMoPo <- rbind(mMo,mPo)
colnames(mMoPo)[2] <- ("Male.reg")
#
fMo <- which(annotatedR$n005.F.XfAlvsXoMt ==1) #  female monogamy
fMo <- cbind(fMo, rep("fMo", length(fMo)))
fPo <- which(annotatedR$n005.F.XfAlvsXoMt ==-1) #  female polygamy
fPo <- cbind(fPo, rep("fPo", length(fPo)))
fMoPo <- rbind(fMo,fPo)
colnames(fMoPo)[2] <- ("Female.reg")
#
mfReg <- merge(mMoPo, fMoPo, by=1, all = TRUE)
#
ConMo <- which(annotatedR$n005.F.XfAlvsXoMt==1 & annotatedR$n005.M.XfAlvsXoMt==1) # concordant monogamy 71
ConMo <- cbind(ConMo,rep("ConMo", length(ConMo)))
ConPo <- which(annotatedR$n005.F.XfAlvsXoMt==-1 & annotatedR$n005.M.XfAlvsXoMt==-1) # concordantpolygamy  17 
ConPo <- cbind(ConPo,rep("ConPo", length(ConPo)))
Dis <- which(annotatedR$n005.F.XfAlvsXoMt==1 & annotatedR$n005.M.XfAlvsXoMt==-1) # discordant  2 (none in other direction)
Dis <- cbind(Dis, rep("Dis", length(Dis)))
#
ConDis <- rbind(ConMo,ConPo,Dis)
colnames(ConDis)[2] <- "Concordance"
SuppTableMating  <- merge(mfReg, ConDis, by= 1, all = TRUE)
SuppTableMating <- SuppTableMating[ order(SuppTableMating[,4], SuppTableMating[,2], SuppTableMating[,3]), ]
dim(SuppTableMating) # [1] 488   4 #### now 426 with Mt104 included 
#########
####  test<- as.character(SuppTableMating[,1]) ### this is necessary so numbers are really numbers...!!!!!!!!!
SuppTableMatingOut <- cbind(SuppTableMating[,c(2,3,4)], annotatedR[as.character(SuppTableMating[,1]),c(1,2,3,5)])
colnames(SuppTableMatingOut) 
# "Male.reg" "Female.reg" "Concordance" "TCorIDannot" "BROADgeneID.x" "BROADgeneTreeInfo1"  "BROADgenePlusInhouse"
dim(SuppTableMatingOut) # 488  7  ## 426 now with Mt104 in
sum(!is.na(SuppTableMatingOut[,5])) #257 ## 224 now with Mt104 in
sum(!is.na(SuppTableMatingOut[,7])) #227  ##194 now with Mt104 in
write.table(SuppTableMatingOut,file="SuppTableMating_20140205.tab",  row.names = FALSE, sep="\t")
#### could add any data (fold change, or P-value) we want to here  ######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


##################################################################
# Lineage (merged based on row number) (use as.character)
## make an X.lineage column a nonX.lineageg column and a concordant column
mX <- which(annotatedR$n005.M.XfXovsAlMt ==1) #  male X-lineage
mX <- cbind(mX, rep("mX", length(mX)))
mNX <- which(annotatedR$n005.M.XfXovsAlMt ==-1) #  male non-X-lineage
mNX <- cbind(mNX, rep("mNX", length(mNX)))
mXNX <- rbind(mX,mNX)
colnames(mXNX)[2] <- ("Male.lin.reg")
#
fX <- which(annotatedR$n005.F.XfXovsAlMt ==1) #  female X-lineage
fX <- cbind(fX, rep("fX", length(fX)))
fNX <- which(annotatedR$n005.F.XfXovsAlMt ==-1) #  female non-X-lineage
fNX <- cbind(fNX, rep("fNX", length(fNX)))
fXNX <- rbind(fX,fNX)
colnames(fXNX)[2] <- ("Female.lin.reg")
#
mfLReg <- merge(mXNX, fXNX, by=1, all = TRUE)
#
ConX <- which(annotatedR$n005.F.XfXovsAlMt==1 & annotatedR$n005.M.XfXovsAlMt==1) # concordant Xlineage 35
ConX <- cbind(ConX,rep("ConX", length(ConX)))
ConNX <- which(annotatedR$n005.F.XfXovsAlMt==-1 & annotatedR$n005.M.XfXovsAlMt==-1) # concordant non-X-lineage  30 
ConNX <- cbind(ConNX,rep("ConNX", length(ConNX)))  
# there are no lineage discordant genes
ConDisL <- rbind(ConX, ConNX)
colnames(ConDisL)[2] <- "Concordance"
SuppTableLineage  <- merge(mfLReg, ConDisL, by= 1, all = TRUE)
SuppTableLineage <- SuppTableLineage[ order(SuppTableLineage[,4], SuppTableLineage[,2], SuppTableLineage[,3]), ]
dim(SuppTableLineage) #    [1] 515   4  #### now 513 with Mt104 in
#########
####  test<- as.character(SuppTableMating[,1]) ### this is necessary so numbers are really numbers...!!!!!!!!!
SuppTableLineageOut <- cbind(SuppTableLineage[,c(2,3,4)], annotatedR[as.character(SuppTableLineage[,1]),c(1,2,3,5)])
colnames(SuppTableLineageOut) 
# "[1] "Male.lin.reg" "Female.lin.reg" "Concordance" "TCorIDannot" "BROADgeneID.x" "BROADgeneTreeInfo1" "BROADgenePlusInhouse"
dim(SuppTableLineageOut) # 515  7  ### 513 now with Mt104 in
sum(!is.na(SuppTableLineageOut[,5])) #280   #### 276 with Mt104 in
sum(!is.na(SuppTableLineageOut[,7])) #232	### 229 with Mt104 in
write.table(SuppTableLineageOut,file="SuppTableLineage_20140205.tab",  row.names = FALSE, sep="\t")
#### could add any data (fold change, or P-value) we want to here


##################################################################
## within lineage mating strategy

mNxM <- which(annotatedR$n005.M.AlvsMt == 1)  
mNxP <- which(annotatedR$n005.M.AlvsMt == -1)
mNxMP <- cbind(c(mNxM, mNxP), c(rep("mNxM",length(mNxM)), rep("mNxP", length(mNxP))))
         
mXM <- which(annotatedR$n005.M.XfvsXo == 1) 
mXP <- which(annotatedR$n005.M.XfvsXo == -1)        
mXMP <- cbind(c(mXM, mXP), c(rep("mXM",length(mXM)), rep("mXP", length(mXP))))


fNxM <- which(annotatedR$n005.F.AlvsMt == 1)   
fNxP <- which(annotatedR$n005.F.AlvsMt == -1)   
fNxMP <- cbind(c(fNxM, fNxP), c(rep("fNxM",length(fNxM)), rep("fNxP", length(fNxP))))
         
fXM <- which(annotatedR$n005.F.XfvsXo == 1)
fXP <- which(annotatedR$n005.F.XfvsXo == -1)
fXMP <- cbind(c(fXM, fXP), c(rep("fXM",length(fXM)), rep("fXP", length(fXP))))

## just work with these and do the discordant in Excel
mXNXMP <- merge(mNxMP, mXMP, by=1, all = TRUE)
fXNXMP <- merge(fNxMP, fXMP, by=1, all = TRUE)
cntcount <- merge(mXNXMP, fXNXMP, by=1, all = TRUE)
temp<- apply(!is.na(cntcount), 1, sum) # make something to use for sorting so those in more than one are at the top
temp2 <- cbind(temp,cntcount)
temp3 <- temp2[order(-temp2[,1],  temp2[,3], temp2[,4], temp2[,5], temp2[,6]),  ]
SuppWithinLineageMating <- cbind(temp3[,c(3,4,5,6)], annotatedR[as.character(temp3[,2]),c(1,2,3,5)])
write.table(SuppWithinLineageMating,file="SuppTableWithinLineageMating_20140205.tab",  row.names = FALSE, sep="\t")


###20140214 Ectodini_GO_Analysis  -- copying from EctoniniGO_20130808
cbind(colnames(dataB))
dim(dataB)#5232   114
# columns I want for analysis of mating system
# [10,] "BROADgeneID"  USED FOR GO
# [57,] "n005.M.XfAlvsXoMt"   
# [109,] "n005.F.XfAlvsXoMt"
alive <- unique(dataB[,10])
length(alive) # 2019 ## though there are 5397 features analyzed on 2019 hit unique genes (didn't collabse by gene hit collapsed by sequence)
## not all of these have a GO term... do they just get ignored?
write.table(alive,file="BINGO_suvived_20140214.txt", sep ="", quote = FALSE, row.names=FALSE )## file for bingo
#### HOW MANY HAVE ANNOTATIONS
GObroad2013 <- read.table (file = "20130807_arrayASSOCIATIONfileSlimBroad.fb", header = FALSE)
head( GObroad2013 )
GObroad2013IDu <- unique(GObroad2013[,2])
length(GObroad2013IDu) # 3777
length(intersect(alive, GObroad2013IDu)) #1419 total genes are annotated to GO
## male
M.ms<-which(dataB$n005.M.XfAlvsXoMt!=0)  
length(M.ms) # 177  # same as before  ### but it shouldn't be should it? well maybe
M.ms.B <- dataB[M.ms,10]
M.ms.Bn <- M.ms.B[which(!is.na(M.ms.B))]
length(M.ms.Bn) # 177  # but before it was 392 according to my notes; but that was beofre compresionn
M.ms.BnU <- unique(M.ms.Bn)
length(M.ms.BnU) # 93 but was 352???  unique genes with annotations (these are things that coulda been compressed)
## female
F.ms <- which(dataB$n005.F.XfAlvsXoMt!=0) 
length(F.ms)  # 331 but was 475
F.ms.B <- dataB[F.ms,10]
F.ms.Bn <- F.ms.B[which(!is.na(F.ms.B))]
length(F.ms.Bn) # 331 but was 261
#summary(F.ms.Bn)  # ~27 are there 2 time, 2 are there 3 times and 1 is there 4 times
F.ms.BnU <- unique(F.ms.Bn)
length(F.ms.BnU) # 159 was 235  unique genes with annotations (these are things that coulda been compressed)
####################### so numbers are very different than August.... but that was before I filtered faint
## male up in monogamy
Mu.ms<-which(dataB$n005.M.XfAlvsXoMt==1)  
length(Mu.ms) # 127  was 473
Mu.ms.B <- dataB[Mu.ms,10]
Mu.ms.Bn <- Mu.ms.B[which(!is.na(Mu.ms.B))]
length(Mu.ms.Bn) # 131 was  242
#summary(Mu.ms.Bn) # 2 with 4  ~ 16 with 2
Mu.ms.BnU <- unique(Mu.ms.Bn)
length(Mu.ms.BnU) # 74 was 221  unique genes with annotations (these are things that coulda been compressed)
## female up in monogamy
Fu.ms <- which(dataB$n005.F.XfAlvsXoMt==1) 
length(Fu.ms)  # 260 was  373
Fu.ms.B <- dataB[Fu.ms,10]
Fu.ms.Bn <- Fu.ms.B[which(!is.na(Fu.ms.B))]
length(Fu.ms.Bn) # 250 was   205
#summary(Fu.ms.Bn)  # 1 with 3 ~ 17 with 2
Fu.ms.BnU <- unique(Fu.ms.Bn)
length(Fu.ms.BnU) # 127   was 186  unique genes with annotations (these are things that coulda been compressed)
#
## male down in monogamy (up in polygamy)
Md.ms<-which(dataB$n005.M.XfAlvsXoMt==-1)  
length(Md.ms) # 46 was  288  ### why are the numbers so much smaller and why were there NAs before that aren't now
Md.ms.B <- dataB[Md.ms,10]
Md.ms.Bn <- Md.ms.B[which(!is.na(Md.ms.B))]
length(Md.ms.Bn) # 46 was 154
#summary(Md.ms.Bn) # 1 with 5 , 1 with 4, 3 with 3  7 with 2
Md.ms.BnU <- unique(Md.ms.Bn)
length(Md.ms.BnU) # 20 was 134  unique genes with annotations (these are things that coulda been compressed)
## female down in monogamy (up in polygamy))
Fd.ms <- which(dataB$n005.F.XfAlvsXoMt==-1) 
length(Fd.ms)  # 71  was ?012
Fd.ms.B <- dataB[Fd.ms,10]
Fd.ms.Bn <- Fd.ms.B[which(!is.na(Fd.ms.B))]
length(Fd.ms.Bn) # 71 was  56
#summary(Fd.ms.Bn)  # 1 with 4 , 1 with 3, 
Fd.ms.BnU <- unique(Fd.ms.Bn)
length(Fd.ms.BnU) # 33 was 51  unique genes with annotations (these are things that coulda been compressed)
### last time I made a table all of same lenght columns 2019
# length(F.ms.BnU) <- 2019
# length(Fu.ms.BnU) <- 2019
# length(Fd.ms.BnU) <- 2019
# length(M.ms.BnU) <- 2019
# length(Mu.ms.BnU) <- 2019
# length(Md.ms.BnU) <- 2019
# length(alive) <- 2019  # already is
### must be a better wya
temp.F.ms.BnU<- as.character(F.ms.BnU)
temp.Fu.ms.BnU<- as.character(Fu.ms.BnU)
temp.Fd.ms.BnU<- as.character(Fd.ms.BnU)
temp.M.ms.BnU<- as.character(M.ms.BnU)
temp.Mu.ms.BnU<- as.character(Mu.ms.BnU)
temp.Md.ms.BnU<- as.character(Md.ms.BnU)
c(temp.Md.ms.BnU)
### why not keep as vector but put batch and proper header between each for direct BINGO input
outGO.ms <-as.character(c("cluster_F.ms.BnU", temp.F.ms.BnU, "batch", "cluster_Fu.ms.BnU", temp.Fu.ms.BnU, "batch", "cluster_Fd.ms.BnU", temp.Fd.ms.BnU, "batch", "cluster_M.ms.BnU", temp.M.ms.BnU, "batch", "cluster_Mu.ms.BnU", temp.Mu.ms.BnU, "batch", "cluster_Md.ms.BnU", temp.Md.ms.BnU) ) ## 
temp <- as.character (outGO.ms) ## 
write.table(temp, file="BINGO_MatingSystem_20140214.txt", sep ="", quote = FALSE, row.names=FALSE )  ## 

### HEATHER DID ALL OF THIS FOR LINEAGE DIFFERENCE ALSO
## make LINEAG bias lists (should this be up and down separate?)
## male
M.l<-which(dataB$n005.M.XfXovsAlMt!=0)  
length(M.l) # 161
M.l.B <- dataB[M.l,10]
M.l.Bn <- M.l.B[which(!is.na(M.l.B))]
length(M.l.Bn) # 161
# summary(M.l.Bn)
M.l.BnU <- unique(M.l.Bn)
length(M.l.BnU) # 89 # was 318  unique genes with annotations (these are things that coulda been compressed)
## female
F.l <- which(dataB$n005.F.XfXovsAlMt!=0) 
length(F.l)  # 431 was 612
F.l.B <- dataB[F.l,10]
F.l.Bn <- F.l.B[which(!is.na(F.l.B))]
length(F.l.Bn) # 431
#summary(F.l.Bn)  # 
F.l.BnU <- unique(F.l.Bn)
length(F.l.BnU) # 215 was 298  unique genes with annotations (these are things that coulda been compressed)
#

## male up in monogamy
Mu.l<-which(dataB$n005.M.XfXovsAlMt==1)  
length(Mu.l) # 92 was 347
Mu.l.B <- dataB[Mu.l,10]
Mu.l.Bn <- Mu.l.B[which(!is.na(Mu.l.B))]
length(Mu.l.Bn) # 92 was 188
Mu.l.BnU <- unique(Mu.l.Bn)
length(Mu.l.BnU) # 53 was 178  unique genes with annotations (these are things that coulda been compressed)
## female up in monogamy
Fu.l <- which(dataB$n005.F.XfXovsAlMt==1) 
length(Fu.l)  # 274 was 368
Fu.l.B <- dataB[Fu.l,10]
Fu.l.Bn <- Fu.l.B[which(!is.na(Fu.l.B))]
length(Fu.l.Bn) # 274 was 190
Fu.l.BnU <- unique(Fu.l.Bn)
length(Fu.l.BnU) # 137 was 174  unique genes with annotations (these are things that coulda been compressed)
#

## male down in monogamy (up in polygamy)
Md.l<-which(dataB$n005.M.XfXovsAlMt==-1)  
length(Md.l) # 69 was 258
Md.l.B <- dataB[Md.l,10]
Md.l.Bn <- Md.l.B[which(!is.na(Md.l.B))]
length(Md.l.Bn) # 69 was 152
Md.l.BnU <- unique(Md.l.Bn)
length(Md.l.BnU) # 37 was 141  unique genes with annotations (these are things that coulda been compressed)
## female down in monogamy (up in polygamy))
Fd.l <- which(dataB$n005.F.XfXovsAlMt==-1) 
length(Fd.l)  # 157 was 244
Fd.l.B <- dataB[Fd.l,10]
Fd.l.Bn <- Fd.l.B[which(!is.na(Fd.l.B))]
length(Fd.l.Bn) # 157 was 133
Fd.l.BnU <- unique(Fd.l.Bn)
length(Fd.l.BnU) # 80 was 126  unique genes with annotations (these are things that coulda been compressed)
###
(dim(intersect(Fd.l.BnU,Fu.l.BnU))) # NULL None are on both the up and down list
(dim(intersect(Md.l.BnU,Mu.l.BnU))) # NULL
###
temp.F.l.BnU<- as.character(F.l.BnU)
temp.Fu.l.BnU<- as.character(Fu.l.BnU)
temp.Fd.l.BnU<- as.character(Fd.l.BnU)
temp.M.l.BnU<- as.character(M.l.BnU)
temp.Mu.l.BnU<- as.character(Mu.l.BnU)
temp.Md.l.BnU<- as.character(Md.l.BnU)
c(temp.Md.l.BnU)
### why not keep as vector but put batch and proper header between each for direct BINGO input
outGO.l <-as.character(c("cluster_F.l.BnU", temp.F.l.BnU, "batch", "cluster_Fu.l.BnU", temp.Fu.l.BnU, "batch", "cluster_Fd.l.BnU", temp.Fd.l.BnU, "batch", "cluster_M.l.BnU", temp.M.l.BnU, "batch", "cluster_Mu.l.BnU", temp.Mu.l.BnU, "batch", "cluster_Md.l.BnU", temp.Md.l.BnU) ) ## 
temp <- as.character (outGO.l) ## FUCK
write.table(temp, file="BINGO_Lineage_20140214.txt", sep ="", quote = FALSE, row.names=FALSE ) 
##### Workin in BINGO to get stats,
##### lood at DAGs to determine which to keep, leef most and anything with genes added to parent and more or less same sig
allSlim <- read.csv(file = "ectodiniGO2014results.csv", header = TRUE, na.strings = "NA", , row.names = 1)
head(allSlim)
allSlimX <- as.matrix(allSlim)
###### 20170706 Reviewers wantto know what genes are in the GO categories.
## problem here is that GO files are suing AB.GENE in caps and my annotations files are in lower case so merging won't work.
## Also, GO formate from BINGO is a nightmare.

class(allSlimX)
dim(allSlimX)
library (gplots)
breaks = c(0.000001, 0.00001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
col <- colorpanel( 8, low = "red", high = "white")
heatmap.2(allSlimX, scale = c("none"), dendrogram = c("none"),Colv = NA, Rowv= NA, col= col,  colsep = c(1:9), rowsep = c(1:23), sepcolor="black", sepwidth=c(0.01,0.01), trace = c("none"), breaks = breaks, key = FALSE,  cexRow = 0.83, main = "GO (SLIM) Ectodini 20140217", na.color="black")



### COUNTING cichlid slim based on BROAD ID
godb <- read.table(file="20130807_arrayASSOCIATIONfileSlimBroad.fb" , header = FALSE, sep = "\t")
head(godb)
alivecb<- data.frame(alive)
head(alivecb)
survivedGO <- merge(godb, alivecb, by.x= c(2), by.y = c(1))
head (survivedGO)  # 
dim(survivedGO)  # [1] 6304  15  different GO annotations for BROAD IDs that survived on array of a total 20758
length(which(survivedGO[,9]=="P")) # 2483 total Biological Process of 7801
length(which(survivedGO[,9]=="C")) # 1532 total Cellular component of 4455
length(which(survivedGO[,9]=="F")) # 2289 total Molecular Funciton of 8502

length(unique(survivedGO[,5])) # 132 total Slim terms interrogated 
length(unique(survivedGO[which(survivedGO[,9]=="P"),5])) # 62 unique Biological Process terms interogated of 178
length(unique(survivedGO[which(survivedGO[,9]=="C"),5])) # 30 unique Cellular Component terms interogated of 66
length(unique(survivedGO[which(survivedGO[,9]=="F"),5])) # 40 unique Molecular Function terms interogated of 97        

length(unique(survivedGO[,1])) # 1419 of 3777 total unique annotated array genes interogated 
length(unique(survivedGO[which(survivedGO[,9]=="P"),1])) # 1063 2813 unique annotated BROAD genes on array  interogated for Biological Process
length(unique(survivedGO[which(survivedGO[,9]=="C"),1])) # 981 of 2555 unique annotated BROAD genes on array  interogated for  Cellular Component
length(unique(survivedGO[which(survivedGO[,9]=="F"),1])) # 1299 of 3467 unique annotated BROAD genes on array interogated for Molecular Function




####2017 reviewer wants to know more about which genes are in teh go terms for different phenotype
## file formats suck here : inconsisten caps,anduse of both , and | for sep :: so super hack to get it in EXCEL added to supplemental tables
