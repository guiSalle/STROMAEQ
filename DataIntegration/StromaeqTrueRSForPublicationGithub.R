## function to perform pre-filtering on OTU data
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.01 # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
bins <- function(xMin,xMax,yMin,yMax,height,width,minBins) {
  if(width > height) {
    hbins = ((width/height)*minBins)
    vbins = minBins
  } else if (width < height) { 
    vbins = ((height/width)*minBins)
    hbins = minBins
  } else { 
    vbins = hbins = minBins
  }
  binwidths <- c(((xMax-xMin)/hbins),((yMax-yMin)/vbins))
  return(binwidths)
}

####===================================================================================
#### STROMAEQ data analysis
####===================================================================================

###====--------------
# flattenCorrMatrix
# https://github.com/heuselm/mocode/blob/master/R/flattenCorrMatrix.R
###====--------------
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut] ## modified
  )
}

###===== Requested packages
require(ggplot2)
require(Hmisc)
require(lme4)
require(lmerTest)
require(ROCR)
require(car)
require(nlme)
require(lsmeans)
require(corrplot)
require(PerformanceAnalytics)
require(glmnet)
library(plyr)
require(tidyr)
require(reshape2)
require(factoextra)
require(FactoMineR)
require(viridis)

###==========----------- Packages for metabolomic data analysis
library(Hotelling)     # Function needed: clr
library(ade4)          # Functions needed: dist.binary, is.euclid
library(vegan)         # Functions needed: rda, stressplot
library(pls)           # Function needed: cppls
library(RVAideMemoire) # Functions needed: MVA.synt, MVA.plot, MVA.cmv, MVA.test
require(mixOmics)
theme_set(theme_bw())
options(digits=3)
###==========-----------

setwd("~/Documents/INRA/STROMAEQ/Paper/")

###==== Define colors
colsRS = viridis_pal(option = 'D')(8)[c(2,6)]
colsStat = viridis_pal(option = 'D')(8)[c(3,7)]
colDay = viridis_pal(option = 'C')(8)[c(2,4,6)]

#paracol = viridis_pal(option = 'D')(8)[c(2,6)]

###==== Load needed data

##--- Parasite and blood data
sumOPG = read.csv(file="../STROMAEQ_OPG.csv",header=T,sep=",") ## FEC
sumHt = read.csv(file="../STROMAEQ_Hemato.csv",header=T,sep=";") ## Blood parameters
sumBio = read.csv(file="../STROMAEQ_Bioch.csv",header=T,sep=";") ## Serum biochemistry

##--- Age data
age = read.csv(file = '../age_nov2014.csv',header=T,sep=';')
age = age[,c('NUM','age')]
head(age)
colnames(age)[1] = 'ID'

sumOPG = merge(age,sumOPG,by=c('ID'))
sumOPG$age = as.numeric(as.character(sumOPG$age))

## Add group information
sumOPG$gp="R"

sumOPG$gp[sumOPG$ID=="W632"] = "S"
sumOPG$gp[sumOPG$ID=="W646"] = "S"
sumOPG$gp[sumOPG$ID=="W652"] = "S"
sumOPG$gp[sumOPG$ID=="W654"] = "S"
sumOPG$gp[sumOPG$ID=="W657"] = "S"
sumOPG$gp[sumOPG$ID=="W666"] = "S"
sumOPG$gp[sumOPG$ID=="W667"] = "S"
sumOPG$gp[sumOPG$ID=="W668"] = "S"
sumOPG$gp[sumOPG$ID=="W672"] = "S"
sumOPG$gp[sumOPG$ID=="W675"] = "S"

sumOPG$Date = as.Date(sumOPG$Date,"%d/%m/%Y")
sumHt$Date = as.Date(sumHt$Date,"%d/%m/%Y")
sumBio$Date = as.Date(sumBio$Date,"%d/%m/%Y")

##--- Weight data
BW = read.table(file='./data/ADG_stromaeq.txt', header=T)
## Weights were measured a day before/after FEC in July and September
## Modify dates so that they match with other traits
BW$Date = as.character(BW$Date)
BW$Date[BW$Date=='09/07/2015'] = '10/07/2015'
BW$Date[BW$Date=='18/09/2015'] = '17/09/2015'

BW$Date = as.Date(BW$Date,"%d/%m/%Y")

colnames(BW)[1] = 'ID'
BW$ID = factor(BW$ID)

head(BW)
#     ID LOTS       Date Time  BW    sex      Tim Value
# 1 W632    S 2015-06-17    0 242 female ADG_0_24 1.409
# 2 W636    R 2015-06-17    0 329 female ADG_0_24 2.218
# 3 W637    R 2015-06-17    0 310 female ADG_0_24 1.018
# 4 W641    R 2015-06-17    0 307 female ADG_0_24 1.627
# 5 W642    R 2015-06-17    0 334 female ADG_0_24 1.182
# 6 W643    R 2015-06-17    0 227 female ADG_0_24 0.918
dim(BW)
#[1] 100   8

adg = BW[,c('ID','Time','BW')]
adg$ID = factor(adg$ID)
w0=0
for(i in unique(adg$ID)){
  w0 = adg$BW[adg$ID==i & adg$Time==0]
  n=0
  for(t in unique(adg$Time)){
    if(t!=0){
    adg$ADG[adg$ID==i & adg$Time==t] = (adg$BW[adg$ID==i & adg$Time==t] - w0)/(t-n)}else{
      adg$ADG[adg$ID==i & adg$Time==t] = 0}
    n = t
  }
  rm(w0)
}

sumADG = merge(BW[,c('ID','Time','Date')], adg, by=c('ID','Time'))
sumADG = sumADG[,c('ID','Date','ADG')]
head(sumADG)
#     ID       Date  BW  ADG
# 1 W632 2015-06-17 242 0.00
# 2 W632 2015-10-26 318 1.90
# 3 W632 2015-07-10 273 1.29
# 4 W632 2015-07-30 282 2.07
# 5 W632 2015-09-17 305 1.28
# 6 W636 2015-06-17 329 0.00

dim(sumADG)
#[1] 100   4

##--- Format parasite and blood data
##--- Merge by anim
dat = merge(sumOPG,sumHt,by=c("ID","Date"))
dim(dat)
#[1] 100  37
dat = merge(dat,sumBio,by=c("ID","Date"))
dim(dat)
#[1] 100  46
dat = merge(dat,sumADG,by=c("ID","Date"))
dim(dat)
#[1] 100  52
head(dat)

##--- Add transformed FEC
shapiro.test(dat$opgSD)
shapiro.test(log(50+dat$opgSD))
dat$FEC = log(50+dat$opgSD)

dat$RSK = 0
dat$RSK[dat$opgSD > 200] = 1
dat$RSK = factor(dat$RSK)
table(dat$RSK,dat$gp)
#    R  S
# 0 49 38
# 1  1 12

## Add day to recode date
dd = array("D",dim(dat)[1])
dat$Day = factor(paste(dd,match(as.factor(dat$Date),levels(factor(dat$Date))),sep=""))

dat$TG = 'unk'
trueR = dat$ID[dat$Day=='D5' & dat$gp=='R' & dat$RSK==0]
trueS = dat$ID[dat$Day=='D5' & dat$gp=='S' & dat$RSK==1]
dat$TG[dat$ID %in% trueR] ='TR'
dat$TG[dat$ID %in% trueS] ='TS'

table(dat$TG,dat$Day)
#    D1 D2 D3 D4 D5
# TR  9  9  9  9  9
# TS  8  8  8  8  8

a = aggregate(opgSD ~ TG, FUN=mean,data=dat[dat$Day=='D5',])
a$std = aggregate(opgSD ~ TG, FUN=sd,data=dat[dat$Day=='D5',])[,2]
a
#    TG opgSD   std
# 1  TR  55.6  76.8
# 2  TS 418.8 148.7
# 3 unk 300.0 477.0

### Those with conflicting potential
falseR = dat$ID[dat$Day=='D5' & dat$gp=='R' & dat$RSK==1]
falseS = dat$ID[dat$Day=='D5' & dat$gp=='S' & dat$RSK==0]
falseOPG = dat[which((dat$ID %in% falseR|dat$ID %in% falseS) & dat$Day=='D5'),
               c(4,7,50,51)]
#    opgSD gp RSK Day
# 15   850  R   1  D5
# 50    50  S   0  D5
# 60     0  S   0  D5

## Time-match with NMR
dat = dat[dat$Day %in% c('D1','D2','D5'),]
dat$RSK = NULL

###-------- Hemato - non-redundant variables 
par(mfrow=c(1,1))
tmp = dat[!is.na(dat$PAL),3:dim(dat)[2]]

tmp$gp = NULL
tmp$RSK = NULL
tmp$FD = NULL
tmp$oxyures = NULL
tmp$opgPEQ = NULL
#tmp$opgSD = NULL
tmp$pLym = NULL
tmp$pMono = NULL
tmp$pBaso = NULL
tmp$pMono = NULL
tmp$pEos = NULL
tmp$pNeut = NULL
tmp$pOther = NULL
tmp$Other = NULL
tmp$Med = NULL
tmp$mode = NULL
## Redundant variables
tmp$IDP = NULL
tmp$IDRSD = NULL
tmp$IDRCV = NULL
tmp$micPlt = NULL ##redundant with VPM
tmp$macPlt = NULL ##redundant with VPM
tmp$GR = NULL ##-- redundant with Ht (91% correlation)
tmp$micGR = NULL ##redundant with VGM
tmp$macGR = NULL ##redundant with VGM
tmp$Pct = NULL
tmp$TCMH = NULL ##-- computed from CCMH
tmp$Hg = NULL ##-- Redundant with Ht

####====------------------- QC on a limited set of variables
lv = colnames(tmp)
lv = c("ID","Date","TG","age",'Day',lv)

### Keep variables of interest
wk = dat[,which(colnames(dat) %in% lv)]
wk2 = melt(wk,c(1:4,dim(wk)[2]))

##--- Reorder data
wk = wk[order(wk$Date,wk$ID),]
wk = wk[wk$TG!='unk',]
wk$Day = factor(wk$Day)
table(wk$Day)
# D1 D2 D5 
# 17 17 17

##--- Non-continuous data: ID, Date, TG, Day
inf = wk[,-c(1:2,4,dim(wk)[2]-1,dim(wk)[2])]
dim(inf)
#[1] 51  22

####====------------------- 16S biodiversity data
# div = read.table(file='../alpha_diversity_animal_time_modified.txt',sep='\t',header=T)
# 
# bdiv = read.table('../matrix.distance.Unifrac.Weighted.D0.txt')
# colnames(bdiv) = rownames(bdiv)
# rownames(bdiv) = div$ID[match(rownames(bdiv),div$SampleID)]
# colnames(bdiv) = div$ID[match(colnames(bdiv),div$SampleID)]
# bdiv2 = as.dist(bdiv[order(rownames(bdiv)),order(colnames(bdiv))])

##--- OTU
metaotu = read.csv(file = '../Mapfie_StromaEq.txt',header=T,sep='\t')

##--- Genus
genus = read.csv(file = './Genus_count_Clarck.tsv', header = T, sep = '\t',strip.white = T)
genus$Genus = as.character(genus$Genus)
## 2 Clostridium genus:
# first one is OTU 8e015d7462eae7d64a0914a922321782 (f__Clostridiaceae)
# second one is OTU d2815209132419e56bf2d8ec485eda4a (f__Lachnospiraceae)
n=''
for(i in grep('Clostridium',genus$Genus)){genus$Genus[i]=paste0('Clostridium',n);n='.1'}
n=''
for(i in grep('Prevotella',genus$Genus)){genus$Genus[i]=paste0('Prevotella',n);n='.1'}

## 2 Prevotella genus:
#                                otuid         Genus
# 55  73278f64dcf2dfc74dc3ad990aac8762    Prevotella (f__Prevotellaceae)
# 103 c4d76fdb0e0e67e2e1d0fb094da01cf5  [Prevotella] (f__[Paraprevotellaceae])

rownames(genus) = genus$Genus
genus = genus[,-1]
genus = data.frame(t(genus))
genus = genus[match(metaotu$SampleID_Sequencing,rownames(genus)),]

### Retain timepoints matching that from NMR
metaotu = metaotu[metaotu$Time %in% c(0,24,132),]
dim(metaotu)
metaotu = metaotu[metaotu$ID %in% wk$ID,]
dim(metaotu)
#[1] 51 10
genus = genus[metaotu$Time %in% c(0,24,132),]
genus = genus[which(rownames(genus) %in% metaotu$SampleID_Sequencing),]
dim(genus)
#[1] 51 91

####====------------------- METABOLOMICS
options(digits=4)
buck = read.csv(file="./data/P2016-0037-bucket-plasma-final.csv",
                header=TRUE,sep=",")
bucketnames = read.csv(file = './data/P2016-0037-buckets-identifies.csv',header=F,sep=';')
bc = buck[,-(seq(1:3))]
#colnames(bc) = bucketnames$V2[match(gsub('X','',colnames(bc)),bucketnames$V1)]
meta = buck[,1:3]
colnames(meta) = c('Ind','Day','Group')
meta$Day = factor(meta$Day)

## Filter for samples correctly predicted
bc = bc[meta$Ind %in% c(as.character(trueR),as.character(trueS)),]
dim(bc)
# 51 791
meta = meta[meta$Ind %in% c(as.character(trueR),as.character(trueS)),]
dim(meta)
#[1] 51  3
### Apply TRUE R/S Group
meta$TG = factor(wk$TG[match(meta$Ind, wk$ID)])
meta$Group = NULL

### Compile lipid shifts
# lipnames = colnames(bc)[grep('Lipides',colnames(bc))]
# lipshift = bucketnames[grep('Lipides',colnames(bc)),]
# dflip = data.frame(lipnames,lipshift)
# write.table(dflip, file='./lipide_name_shifts.tsv',quote = F,row.names = F)

### Output individual bucket data for TR and TS individuals only
bucketindiv = cbind(meta,bc)
colnames(bucketindiv) = gsub('X','',colnames(bucketindiv))
write.csv(bucketindiv,file = './RawBucketData_TR_TS_byInd_byDay.csv',quote=F)

### Remove weird characters imported during csv reading
colnames(bc) = bucketnames$V2[match(gsub('X','',colnames(bc)),bucketnames$V1)]
peakToKept = which(colnames(bc) != 'Bruit' & colnames(bc) != 'bruit')

### Find noise
as.numeric(as.character(substr(colnames(buck)[which(colnames(bc)=='Bruit')][-c(1:3)],2,6)))

dim(bc)
#[1]  51 791

### Add number to identifcal peak names
bc = data.frame(bc)

### Denoise data
bc.dnoised = bc[,peakToKept] 
dim(bc.dnoised)
#[1]  51 356

### Find unknown
length(grep('U',colnames(bc)))
# 36
as.numeric(as.character(substr(colnames(buck)[grep('U',colnames(bc))][-c(1:3)],2,6)))
# [1] 7.915 7.885 7.875 7.865 7.155 7.145 7.135 7.125 7.115 7.105 7.095 6.845 6.835 4.335 4.325
# [16] 4.075 4.065 4.055 3.415 3.185 1.465 1.455 1.445 1.435 1.425 1.415 1.165 1.155 1.145 0.945
# [31] 0.935 0.925 0.915

###---=== Compile Metabolite id with their respective shifts

metabo_names_analysis = colnames(bc) ## metabolite peak names used in analysis
shifts = bucketnames$V1 ## peak shift
bucketnames$peak = colnames(bc)

## Remove noise
bucketnames = bucketnames[bucketnames$V2!='Bruit' & bucketnames$V2!='bruit',]

## Assign shift interval for each metabolite
rep = 0
k = 0
metabo = NULL
for(i in 1:dim(bucketnames)[1]){
  if(rep==0 | bucketnames$V2[i]!=rep){
    rep = bucketnames$V2[i]
    k = k + 1
    metabo$metabolite[k] = bucketnames$V2[i]
    metabo$down[k] = bucketnames$V1[i]
    metabo$up[k] = bucketnames$V1[i]
  }else{
    metabo$up[k] = bucketnames$V1[i]
  }
}
metabo = data.frame(metabo)
metabo$metabolite = gsub('\\+',' ',metabo$metabolite)

##Expand table to include metabolite peak names for re-annotation of lipid-associated peaks and unknown spectra
bucketnames$V2 = gsub('\\+',' ',bucketnames$V2)
for(i in 1:dim(bucketnames)[1]){
  #met = bucketnames$V2[i]
  peak = bucketnames$V1[i]
  for(j in 1:dim(metabo)[1]){
    if(peak <= metabo$down[j] & peak >= metabo$up[j]){
      bucketnames$metabo[i] = metabo$metabolite[j]
      bucketnames$up[i] = metabo$up[j]
      bucketnames$down[i] = metabo$down[j]
    }
  }
}

nmrindv = t(bc.dnoised)
nmrindv = cbind(bucketnames,nmrindv)
meta2 = meta
meta2$Day = as.character(meta2$Day)
meta2$Day[meta2$Day=='J24']='d24'
meta2$Day[meta2$Day=='RB']='d132'
meta2$Day[meta2$Day=='J0']='d0'
colnames(nmrindv)[7:dim(nmrindv)[2]] = paste(meta2$Ind,meta2$Day,meta2$TG,sep='_')
colnames(nmrindv)[1:6] = c('Shift','Metabolite','Peak name','Metabolite2','Lower','Upper')

write.csv(nmrindv,file = './Table_HNMR_annotated.csv',quote = F,row.names = F)

##Import new annotation and update bc.dnoised colnames - Done 19/02/2021
bucknames2 = read.csv(file = './data/Table_HNMR_reannotated_NM.csv',header = T)
bucknames2$X = NULL
colnames(bucknames2) = gsub('RB','d132',colnames(bucknames2))
head(bucknames2[,1:10])
#   Shift Metabolite Peak.name Metabolite_Reannotated Metabolite_name_R Lower Upper W632_d0_TS
# 1 8.225         U1        U1                Inosine                   8.215 8.225  0.0001421
# 2 8.215         U2        U2                Inosine                   8.215 8.225  0.0001293
# 3 7.995  Histidine Histidine              Histidine                   7.985 7.995  0.0004731
# 4 7.985  Histidine Histidine              Histidine                   7.985 7.995  0.0003869
# 5 7.895         U3        U3                     U3                U3 7.895 7.895  0.0001885
# 6 7.885         U4        U4                     U4                U4 7.885 7.885  0.0000758
#   W632_d24_TS W632_d132_TS
# 1   0.0001199    0.0001301
# 2   0.0001883    0.0001672
# 3   0.0002798    0.0002675
# 4   0.0004407    0.0003281
# 5   0.0002014    0.0001885
# 6   0.0001188    0.0001022

##Replace missing symbols and polish names
bucknames2$Metabolite_Reannotated = gsub('_','alpha',bucknames2$Metabolite_Reannotated)
bucknames2$Metabolite_Reannotated = gsub(' ','-',bucknames2$Metabolite_Reannotated)
bucknames2$Metabolite_Reannotated = stringr::str_to_title(bucknames2$Metabolite_Reannotated)
bucknames2$Metabolite_Reannotated = gsub('Lipides','Lipid',bucknames2$Metabolite_Reannotated)
bucknames2$Metabolite_Reannotated = gsub('Lipides','Lipid',bucknames2$Metabolite_Reannotated)
bucknames2$Metabolite_Reannotated = gsub('Tryptophane','Tryptophan',bucknames2$Metabolite_Reannotated)

##Renumber unknown metabolites
for(i in 1:length(bucknames2$Metabolite_Reannotated[grep('U',bucknames2$Metabolite_Reannotated)])){
  bucknames2$Metabolite_Reannotated[grep('U',bucknames2$Metabolite_Reannotated)][i] = paste0('U',i)
}

### Find unknown metabolites
length(grep('U',bucknames2$Metabolite_Reannotated))
# 29
bucknames2$Shift[grep('U',bucknames2$Metabolite_Reannotated)]
# [1] 7.895 7.885 7.855 7.845 7.835 7.125 7.115 7.105 7.095 7.085 7.075 7.065 6.815 6.805 4.305
# [16] 4.295 4.045 4.035 4.025 3.385 1.435 1.425 1.415 1.405 1.395 1.385 1.135 1.125 1.115

## Assign new metabolite id to the working NMR set
colnames(bc.dnoised) = bucknames2$Metabolite_Reannotated
head(bc.dnoised[,1:10])
#    Inosine Inosine.1 Histidine Histidine.1       U3       U4       U5       U6       U7 1-methylhistidine
# 1 0.000142  1.29e-04  0.000473    0.000387 0.000188 7.58e-05 3.27e-05 5.39e-05 2.96e-05          8.02e-05
# 2 0.000120  1.88e-04  0.000280    0.000441 0.000201 1.19e-04 6.28e-05 9.92e-05 7.81e-05          2.24e-05
# 3 0.000130  1.67e-04  0.000267    0.000328 0.000189 1.02e-04 7.47e-05 1.00e-04 9.33e-05          2.78e-05
# 4 0.000191  6.89e-05  0.000657    0.000282 0.000198 9.96e-05 4.77e-05 5.64e-05 5.86e-05          1.08e-04
# 5 0.000181  1.42e-04  0.000445    0.000397 0.000221 1.24e-04 1.24e-04 1.56e-04 1.48e-04          4.54e-05
# 6 0.000139  1.22e-04  0.000322    0.000274 0.000177 1.09e-04 7.84e-05 1.12e-04 7.89e-05          3.30e-05

##Add numbers to identical metabolite peak names
bc.dnoised = data.frame(bc.dnoised) 
colnames(bc.dnoised) = gsub('X','',colnames(bc.dnoised))

## Use colnames to update bucknames2
bucknames2$Peak.name = colnames(bc.dnoised)

bucknames2$Metabolite_name_R = NULL
head(bucknames2[,1:7])
#   Shift Metabolite   Peak.name Metabolite_Reannotated Lower Upper W632_d0_TS
# 1 8.225         U1     Inosine                Inosine 8.215 8.225  0.0001421
# 2 8.215         U2   Inosine.1                Inosine 8.215 8.225  0.0001293
# 3 7.995  Histidine   Histidine              Histidine 7.985 7.995  0.0004731
# 4 7.985  Histidine Histidine.1              Histidine 7.985 7.995  0.0003869
# 5 7.895         U3          U1                     U1 7.895 7.895  0.0001885
# 6 7.885         U4          U2                     U2 7.885 7.885  0.0000758


####-----====== Working set = Aggregate buckets by metabolite peaks =====------
meta_peaks = unique(paste0(bucknames2$Lower,'-',bucknames2$Upper))
length(meta_peaks)
#[1] 117

length(unique(bucknames2$Metabolite_Reannotated))
#[1] 87

## Melt nmr peak
buck_melt = reshape2::melt(bucknames2,c(1:6))
buck_melt$meta_peak = paste0(bucknames2$Lower,'-',bucknames2$Upper)
head(buck_melt)
#   Shift Metabolite   Peak.name Metabolite_Reannotated Lower Upper   variable     value
# 1 8.225         U1     Inosine                Inosine 8.215 8.225 W632_d0_TS 0.0001421
# 2 8.215         U2   Inosine.1                Inosine 8.215 8.225 W632_d0_TS 0.0001293
# 3 7.995  Histidine   Histidine              Histidine 7.985 7.995 W632_d0_TS 0.0004731
# 4 7.985  Histidine Histidine.1              Histidine 7.985 7.995 W632_d0_TS 0.0003869
# 5 7.895         U3          U1                     U1 7.895 7.895 W632_d0_TS 0.0001885
# 6 7.885         U4          U2                     U2 7.885 7.885 W632_d0_TS 0.0000758
#     meta_peak
# 1 8.215-8.225
# 2 8.215-8.225
# 3 7.985-7.995
# 4 7.985-7.995
# 5 7.895-7.895
# 6 7.885-7.885

dim(buck_melt)
#[1] 18156     9

### Sum intensities for metabolite peaks
sum_meta = aggregate(value ~ meta_peak + variable, FUN = sum, data = buck_melt)
head(sum_meta)
#     meta_peak   variable    value
# 1 0.805-0.875 W632_d0_TS 0.071839
# 2 0.885-0.915 W632_d0_TS 0.024219
# 3 0.935-0.945 W632_d0_TS 0.006224
# 4 0.955-0.955 W632_d0_TS 0.003652
# 5 0.965-0.975 W632_d0_TS 0.006510
# 6 0.985-0.995 W632_d0_TS 0.003242

dim(sum_meta)
#[1] 5967    3

### Merge metabolite data with summed intensities
merged_buck = unique(merge(buck_melt[,-c(1,2,3,5,6,8)],sum_meta,by=c('variable','meta_peak')))
dim(merged_buck)
#[1] 6069    4

head(merged_buck)
#      variable   meta_peak  Metabolite_Reannotated    value
# 1  W632_d0_TS 0.805-0.875 Methyl-Moieties-From-Fa 0.071839
# 9  W632_d0_TS 0.885-0.915                Butyrate 0.024219
# 13 W632_d0_TS 0.935-0.945              Isoleucine 0.006224
# 15 W632_d0_TS 0.955-0.955      Leucine-Isoleucine 0.003652
# 16 W632_d0_TS 0.965-0.975                 Leucine 0.006510
# 18 W632_d0_TS 0.985-0.995                  Valine 0.003242

nmr_peak = merged_buck %>% dcast(meta_peak + Metabolite_Reannotated ~ variable, 
                           value.var = 'value')
dim(nmr_peak)
#[1] 119  54

head(nmr_peak[,1:4])
#     meta_peak  Metabolite_Reannotated W632_d0_TS W632_d24_TS
# 1 0.805-0.875 Methyl-Moieties-From-Fa   0.071839    0.088622
# 2 0.885-0.915                Butyrate   0.024219    0.022367
# 3 0.935-0.945              Isoleucine   0.006224    0.007141
# 4 0.955-0.955      Leucine-Isoleucine   0.003652    0.004128
# 5 0.965-0.975                 Leucine   0.006510    0.006415
# 6 0.985-0.995                  Valine   0.003242    0.004412

## ####------===== Output final supplementary Table 1
nmr_peak = as.matrix(nmr_peak)
rownames(nmr_peak) = nmr_peak[,2]
nmr_peak = data.frame(nmr_peak)
rownames(nmr_peak) = gsub('X','',rownames(nmr_peak))
nmr_peak$Metabolite_Reannotated = rownames(nmr_peak)
colnames(nmr_peak)[1]='Chemical_shift'
write.csv(nmr_peak,file = './Supplementary_Table1.csv',quote=F,row.names=F)

#### working NMR set
bc.dnoised = t(as.matrix(nmr_peak[,-c(1,2)]))
bc.dnoised = data.frame(bc.dnoised)
colnames(bc.dnoised) = gsub('X','',colnames(bc.dnoised))
write.csv(bc.dnoised,file='./bc.dnoised.csv',quote=F,row.names = T)
bc.dnoised = as.matrix(sapply(bc.dnoised,as.numeric))

###=====================================================================================
###= FEC and biochemical data analysis
###=====================================================================================

###============--------- Clinical data trajectory and predictive ability for FEC

# ggplot(dat, aes(x = factor(Date), y = opgSD, color = gp,group = ID)) + 
#   geom_point() + geom_line() +
#   scale_color_manual(values = colsRS) +
#   facet_wrap( ~ ID, nrow=4) +		
#   #  scale_x_discrete(name="Day ", breaks = seq(1:length(unique(tstFec$DATE))), 
#   #                   labels = unique(tstFec$DATE)) +
#   scale_y_continuous(name="OPG", limits = c(0, max(dat$opgSD))) +
#   ggtitle("OPG") +
#   theme(plot.title = element_text(lineheight=.8, face="bold"))
# 
# ggplot(dat, aes(x = factor(Date), y = opgSD, color = gp)) + 
#   geom_boxplot() +
#   scale_color_manual(values = colsRS) +
#   #facet_wrap( ~ ID, nrow=4) +		
#   #  scale_x_discrete(name="Day ", breaks = seq(1:length(unique(tstFec$DATE))), 
#   #                   labels = unique(tstFec$DATE)) +
#   scale_y_continuous(name="OPG", limits = c(0, max(dat$opgSD))) +
#   ggtitle("OPG") +
#   theme(plot.title = element_text(lineheight=.8, face="bold")) +
#   theme(axis.text=element_text(size=12),
#         legend.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         title=element_text(size=10,face="bold"))
# 
# ### Correlation between clinical variables
# cormat = rcorr(as.matrix(inf))
# corrplot::corrplot(cormat$r, type="upper", order="hclust", 
#                    tl.col="black", tl.srt=45)
# 
# cormat2 = rcorr(as.matrix(inf),type='spearman')
# corrplot::corrplot(cormat2$r, type="upper", order="hclust", 
#                    tl.col="black", tl.srt=45)

### Relationship between FEC and clinical variables
# inf %>%
#   gather(-FEC, key = "variable", value = "value") %>% 
#   ggplot(aes(x = value, y = FEC)) +
#   geom_point(size =.8) +
#   facet_wrap(~ variable, scales = "free") 

### Age has an impact on a few features: ADG, BW and Lym <> growth vs. parasite exposure

# wk[,c('age','Lym','ADG')] %>%
#   gather(- age, key = "variable", value = "value") %>% 
#   ggplot(aes(x = value, y = age)) +
#   geom_point(size = 2) +
#   facet_wrap(~ variable, scales = "free") 
# 
# inf %>%
#   gather(-age, key = "variable", value = "value") %>% 
#   ggplot(aes(x = value, y = age)) +
#   geom_point(size =.8) +
#   facet_wrap(~ variable, scales = "free") 

# ## Age does not impact on FEC
# 
# # inf %>%
# #   gather(-opgSD, key = "variable", value = "value") %>% 
# #   ggplot(aes(x = value, y = opgSD)) +
# #   geom_point(size = .8) +
# #   facet_wrap(~ variable, scales = "free") 
# # 
# # ## Plot density
# # inf %>%
# #   gather(-opgSD, key = "variable", value = "value") %>% 
# #   ggplot(aes(x = as.numeric(as.character(value)))) +
# #   geom_density() + xlab('Variable') + 
# #   facet_wrap(~ variable, scales = "free") 

## Plot clinical variable ~ density
# inf %>%
# #  keep(is.numeric) %>%                     # Keep only numeric columns
#   gather() %>%                             # Convert to key-value pairs
#   ggplot(aes(value)) +                     # Plot the values
#   facet_wrap(~ key, scales = "free") +   # In separate panels
#   geom_density()                         # as density

###--- Check normality for every variable
norm = as.data.frame(sapply(inf[,-1],shapiro.test))[-c(3,4),]

norm[,which(norm[1,]<0.9)]
#                Eos    Baso      PAL  Glucose  ABLGLOB      ADG      FEC
# statistic    0.856   0.649    0.791    0.747    0.751    0.893    0.635
# p.value   1.92e-05 8.8e-10 4.42e-07 5.25e-08 6.37e-08 0.000256 5.21e-10

######==============================================================================================================
###########----------- Feature engineering: Replace old variables by transformed ones & rm Basophils (no variation)
######==============================================================================================================

###-- Create new variable for non-FEC variables and Baso (not retained as only 2 values)
for(i in which(norm[1,]<0.9)){
  assign(paste0(colnames(inf)[i+1],"log"),log(inf[,i+1]+1))
}

##-- Remove useless info
inf$Baso = NULL ## no variation

##--Replace variables by their normalized counterparts

inf$lEos = Eoslog
inf$lPAL = PALlog
inf$lGlc = Glucoselog
inf$lABGL = ABLGLOBlog

inf$Eos = NULL
inf$PAL = NULL
inf$Glucose = NULL
inf$ABLGLOB = NULL

###--- Check normality following data transformation and sample removal
as.data.frame(sapply(inf,shapiro.test))[-c(3,4),]
# age     GB   Lym    Mono    Neut   VGM    Ht  CCMH   Plt     VPM   UREE Albumine
# statistic    0.896  0.953  0.98   0.938   0.915 0.982 0.974 0.987 0.982   0.934  0.925    0.952
# p.value   0.000318 0.0407 0.525 0.00981 0.00138 0.613 0.335 0.857 0.612 0.00729 0.0033   0.0371
# Cholesterol    PT Globines      ADG      FEC   lEos    lPAL  lGlc    lABGL
# statistic       0.964 0.983    0.968    0.893    0.635  0.953   0.925 0.857    0.897
# p.value         0.125 0.684    0.175 0.000256 5.21e-10 0.0418 0.00327 2e-05 0.000324

## New working set with appropriate variables
clinical = cbind(wk[,c('ID','TG','Day')],inf)
clinical$TG = factor(clinical$TG)
clinical$ID = factor(clinical$ID)

head(clinical)
#      ID TG Day age   GB Lym Mono Neut  VGM   Ht CCMH Plt VPM UREE Albumine Cholesterol PT
# 1  W632 TS  D1   6 7.21 3.7  0.5  2.6 44.0 34.4 33.9  67 7.2 0.21       24         0.9 63
# 6  W636 TR  D1   6 7.55 3.5  0.5  2.9 45.6 32.2 34.2 136 6.2 0.21       25         0.8 54
# 16 W641 TR  D1   6 6.83 3.5  0.4  2.7 40.4 29.5 34.8 113 6.9 0.26       25         0.7 55
# 21 W642 TR  D1   5 7.03 3.7  0.4  2.5 44.0 37.7 31.2 151 6.4 0.17       24         0.8 58
# 26 W643 TR  D1   5 6.53 4.4  0.8  0.4 39.4 27.9 32.7  46 8.1 0.15       25         0.9 53
# 31 W646 TS  D1   5 7.01 3.8  0.4  2.3 43.7 36.6 33.0  53 7.3 0.24       23         0.5 63
# Globines ADG  FEC  lEos lPAL  lGlc lABGL
# 1        38   0 3.91 0.262 5.26 0.693 0.470
# 6        29   0 3.91 0.405 5.34 0.610 0.588
# 16       30   0 3.91 0.182 5.40 0.495 0.588
# 21       34   0 3.91 0.262 5.16 0.571 0.531
# 26       28   0 3.91 0.531 5.37 0.593 0.642
# 31       40   0 3.91 0.336 5.00 0.560 0.470

#######-------------------- Apply age correction to clinical dataset if needed -----------------------###########
age.vec = clinical$age
mod.stor = NULL

for(i in 6:dim(clinical)[2]){
  y = clinical[,i]
  p = data.frame(summary(lm(y ~ age.vec, data = clinical))$coefficients)
  p$param = colnames(clinical)[i]
  mod.stor = rbind(mod.stor,p[2,])
}
p = mod.stor$Pr...t..
mod.stor$fdr = p.adjust(p, method = 'BH', n = length(p))

mod.stor[mod.stor$fdr < .05,]
# Estimate Std..Error t.value Pr...t.. param      fdr
# age.vec   -0.468     0.0688    -6.8 1.35e-08   Lym 2.56e-07

## Before correction
clinical[clinical$Day=='D1',c('age','Lym')] %>%
  gather(- age, key = "variable", value = "value") %>% 
  ggplot(aes(x = value, y = age)) +
  geom_point(size = 2) +
  facet_wrap(~ variable, scales = "free")

### Apply correction for identified parameters
list.param = mod.stor$param[mod.stor$fdr < .05]
n = 0
for(i in match(list.param, colnames(clinical))){
  n = n + 1
  clinical = cbind(clinical,residuals(lm(clinical[,i] ~ age.vec)))
  colnames(clinical)[dim(clinical)[2]] = paste0(list.param[n],'c')
}

##check new relationship
clinical[,c('age','Lymc')] %>%
  gather(- age, key = "variable", value = "value") %>% 
  ggplot(aes(x = value, y = age)) +
  geom_point(size = 2) +
  facet_wrap(~ variable, scales = "free") 

##-- Remove uncorrected variables
clinical = clinical[,-c(match(list.param, colnames(clinical)))]
dim(clinical)
#[1] 51  24

head(clinical)
#      ID TG Day age   GB Mono Neut  VGM   Ht CCMH Plt VPM UREE Albumine Cholesterol PT Globines
# 1  W632 TS  D1   6 7.21  0.5  2.6 44.0 34.4 33.9  67 7.2 0.21       24         0.9 63       38
# 6  W636 TR  D1   6 7.55  0.5  2.9 45.6 32.2 34.2 136 6.2 0.21       25         0.8 54       29
# 16 W641 TR  D1   6 6.83  0.4  2.7 40.4 29.5 34.8 113 6.9 0.26       25         0.7 55       30
# 21 W642 TR  D1   5 7.03  0.4  2.5 44.0 37.7 31.2 151 6.4 0.17       24         0.8 58       34
# 26 W643 TR  D1   5 6.53  0.8  0.4 39.4 27.9 32.7  46 8.1 0.15       25         0.9 53       28
# 31 W646 TS  D1   5 7.01  0.4  2.3 43.7 36.6 33.0  53 7.3 0.24       23         0.5 63       40
#    ADG  FEC  lEos lPAL  lGlc lABGL   Lymc
# 1    0 3.91 0.262 5.26 0.693 0.470  0.332
# 6    0 3.91 0.405 5.34 0.610 0.588  0.132
# 16   0 3.91 0.182 5.40 0.495 0.588  0.132
# 21   0 3.91 0.262 5.16 0.571 0.531 -0.136
# 26   0 3.91 0.531 5.37 0.593 0.642  0.564
# 31   0 3.91 0.336 5.00 0.560 0.470 -0.036

##-- Final check on correlation between features
cormat4 = rcorr(as.matrix(clinical[,-seq(1:3)]),type='pearson')
corrplot::corrplot(cormat4$r, type="upper", order="hclust", 
                   tl.col="black", tl.srt=45)

summary(cormat4$r[cormat4$r!=1])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  -0.799  -0.033   0.146   0.161   0.326   0.838  

dim(clinical)
# dim(clinical)
# [1] 51  24

###----- Rename clinical variables - switch to English and annotations used in paper
colnames(clinical)
# [1] "ID"          "TG"          "Day"         "age"         "GB"          "Mono"        "Neut"        "VGM"         "Ht"         
# [10] "CCMH"        "Plt"         "VPM"         "UREE"        "Albumine"    "Cholesterol" "PT"          "Globines"    "ADG"        
# [19] "FEC"         "lEos"        "lPAL"        "lGlc"        "lABGL"       "Lymc"

colnames(clinical) = c("ID","TG","Day","Age",
                       "Leukocytes","Monocytes","Neutrophil",
                       "MGV","Haematocrit",         
                       "MCHC","Platelet","MPV",
                       "Urea",
                       "Albumin","Cholesterol", "TP","Globins",
                       "ADG","FEC","Eosinophils","ALP","Glucose","AGR","Lymphocytes")

##----==== Output supplementary Table 4 ====----
write.csv(clinical,
          file = './Supplementary_Table4.csv', quote = F, row.names = F)

#######-------------------- DATA ready for analysis -----------------------###########

###---==== Differences in clinical parameters between R and S individuals at D0
mod0clinic = NULL
nclinical = length(seq(5,dim(clinical)[2]-1,1))
d0 = clinical[clinical$Day=='D1',-c(1:3)]
TG = clinical$TG[clinical$Day=='D1']
d0$FEC = NULL
d0$RSK = NULL

for(i in seq(2,dim(d0)[2]-1,1)){
  test = t.test(d0[,i] ~ TG)
  tt = test$p.value
  tt = data.frame(tt)
  tt$var = colnames(d0)[i]
  mod0clinic = rbind(mod0clinic,tt)
}
na.omit(mod0clinic[mod0clinic$tt<.05,])
# [1] tt  var
# <0 rows> (or 0-length row.names)

###---==== Run mixed model on every clinical data (longitudinal trend)
modclinic = NULL
nclinical = length(seq(6,dim(clinical)[2]-1,1))

## model with predicted groups
for(i in seq(6,dim(clinical)[2]-1,1)){
  y = clinical[,i] 
  mod = lme(y ~ TG*Day,
            random =~ 1|ID,
            data = clinical) 
  ano = as.data.frame(anova.lme(mod,type="marginal", adjustSigma = F))
  ano$param = colnames(clinical)[i]
  ano$fact = rownames(ano)
  modclinic = rbind(modclinic,ano)
  rm(ano,mod)
}
rownames(modclinic) = seq(1,dim(modclinic)[1],1)
p = modclinic$`p-value`[modclinic$fact=='TG']
modclinic$fdr[modclinic$fact=='TG'] = p.adjust(p, method = 'BH', n = length(p))

p = modclinic$`p-value`[modclinic$fact=='TG:Day']
modclinic$fdr[modclinic$fact=='TG:Day'] = p.adjust(p, method = 'BH', n = length(p))

p = modclinic$`p-value`[modclinic$fact=='Day']
modclinic$fdr[modclinic$fact=='Day'] = p.adjust(p, method = 'BH', n = length(p))

## Significant differences
modclinic[modclinic$fact == 'TG' & modclinic$fdr<0.05,] 
# [1] numDF   denDF   F-value p-value param   fact    fdr    
# <0 rows> (or 0-length row.names)

modclinic[modclinic$fact == 'TG:Day' & modclinic$fdr < 0.05,]
# numDF denDF F-value  p-value param   fact      fdr
# 56     2    30    18.2 6.63e-06   FEC TG:Day 0.000119

modclinic[modclinic$fact == 'TG:Day' & modclinic$`p-value` < 0.05,]
# numDF denDF F-value  p-value param   fact      fdr
# 56     2    30   18.21 6.63e-06                 FEC TG:Day 0.000119
# 64     2    30    3.49 4.34e-02 AlkalinePhosphatase TG:Day 0.390876

dim(modclinic[modclinic$fact == 'Day' & modclinic$fdr<0.05,])
#[1] 14  7
modclinic[modclinic$fact == 'Day' & modclinic$fdr<0.05,]
#    numDF denDF F-value  p-value               param fact      fdr
# 3      2    30   18.16 6.81e-06           Monocytes  Day 1.75e-05
# 7      2    30    3.89 3.14e-02          Neutrophil  Day 4.34e-02
# 11     2    30   28.24 1.27e-07  MeanGlobularVolume  Day 5.71e-07
# 15     2    30   14.92 3.17e-05         Haematocrit  Day 7.14e-05
# 23     2    30   10.01 4.66e-04            Platelet  Day 9.32e-04
# 31     2    30   38.01 5.97e-09                Urea  Day 5.37e-08
# 35     2    30    7.89 1.77e-03             Albumin  Day 3.18e-03
# 39     2    30    4.80 1.55e-02         Cholesterol  Day 2.32e-02
# 43     2    30   18.29 6.41e-06        TotalProtein  Day 1.75e-05
# 47     2    30    3.63 3.87e-02             Globins  Day 4.98e-02
# 51     2    30  176.83 0.00e+00                 ADG  Day 0.00e+00
# 55     2    30    5.64 8.36e-03                 FEC  Day 1.37e-02
# 59     2    30   34.36 1.74e-08         Eosinophils  Day 1.04e-07
# 63     2    30   18.37 6.19e-06 AlkalinePhosphatase  Day 1.75e-05
# > 

## Check model details for FECc, lEos, lPAL
y = clinical[,c('FEC')]
summary(lme(y ~ TG*Day ,
            random =~ 1|ID,
            data = clinical))
# Linear mixed-effects model fit by REML
# Data: clinical 
# AIC  BIC logLik
# 85.2 99.7  -34.6
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept) Residual
# StdDev:       0.277    0.388
# 
# Fixed effects: y ~ TG * Day 
#             Value Std.Error DF t-value p-value
# (Intercept)  3.91     0.159 30   24.62  0.0000
# TGTS         0.17     0.232 15    0.75  0.4659
# DayD2        0.00     0.183 30    0.00  1.0000
# DayD5        0.53     0.183 30    2.91  0.0068
# TGTS:DayD2   0.22     0.267 30    0.84  0.4077
# TGTS:DayD5   1.49     0.267 30    5.60  0.0000
# Correlation: 
#            (Intr) TGTS   DayD2  DayD5  TGTS:DD2
# TGTS       -0.686                              
# DayD2      -0.576  0.395                       
# DayD5      -0.576  0.395  0.500                
# TGTS:DayD2  0.395 -0.576 -0.686 -0.343         
# TGTS:DayD5  0.395 -0.576 -0.343 -0.686  0.500  
# 
# Standardized Within-Group Residuals:
#   Min      Q1     Med      Q3     Max 
# -1.0947 -0.5779 -0.0836  0.2758  2.6841 
# 
# Number of Observations: 51
# Number of Groups: 17 

aggregate(opgSD ~ TG*Day, data= wk[wk$Day =='D5',],FUN = median)
#   TG Day opgSD
# 1  TR  D5    0
# 2  TS  D5   375

length(which(wk$opgSD[wk$Day=='D5' & wk$TG=='TR']>0))
# [1] 4

y = clinical[,c('ALP')]
# summary(lme(y ~ TG*Day ,
#             random =~ 1|ID,
#             data = clinical))
# Linear mixed-effects model fit by REML
# Data: clinical 
# AIC  BIC logLik
# 10 24.4      3
# 
# Random effects:
#   Formula: ~1 | ID
# (Intercept) Residual
# StdDev:       0.139    0.162
# 
# Fixed effects: y ~ TG * Day 
#             Value Std.Error DF t-value p-value
# (Intercept)  5.30    0.0710 30    74.6   0.000
# TGTS        -0.04    0.1035 15    -0.4   0.730
# DayD2        0.44    0.0762 30     5.8   0.000
# DayD5        0.10    0.0762 30     1.3   0.203
# TGTS:DayD2  -0.29    0.1111 30    -2.6   0.014
# TGTS:DayD5  -0.11    0.1111 30    -0.9   0.350
# Correlation: 
#   (Intr) TGTS   DayD2  DayD5  TGTS:DD2
# TGTS       -0.686                              
# DayD2      -0.537  0.368                       
# DayD5      -0.537  0.368  0.500                
# TGTS:DayD2  0.368 -0.537 -0.686 -0.343         
# TGTS:DayD5  0.368 -0.537 -0.343 -0.686  0.500  
# 
# Standardized Within-Group Residuals:
#   Min      Q1     Med      Q3     Max 
# -1.9122 -0.5184 -0.0189  0.4038  1.9145 
# 
# Number of Observations: 51
# Number of Groups: 17 

######============================================================================================
#######----------------------------- Are TR and TS statuses stable through time
######============================================================================================

dir='~/Documents/INRA/HorseDrugResistance/LongTermTST/'

###--- Check whether FEC data after the experiment match predicted status
strom = c(trueR,trueS)

### This file is not disclosed unless an MTA has been agreed upon 
fec = read.csv(file = paste0(dir,'FEC.csv'),header=T,sep=';')
dim(fec)
#[1] 6009   10
fec = fec[!is.na(fec$RESULTATS_OPG),]
dim(fec)
#[1] 5928   10
fec$DATE_PRELEVEMENT = as.Date(fec$DATE_PRELEVEMENT,'%d/%m/%Y')
## Individual FEC only
fec = fec[fec$MODE_DE_PRELEVEMENT=='individuel',]
## Strongylid only
fec = fec[!(fec$VERS %in% c('Anoplocephale_opg','Anoplocephales_vus','oxyures',
                            'OXYURES','oxyures prEsence(1) absence(0)',
                            'P_equorum_opg','P_equorum_vus','PASCARIS')),]
dim(fec)
#[1] 4596   10

fec = fec[fec$RESULTATS_OPG <10000,] ## 1 outlier removed
dim(fec)
#[1] 4595   10

### Put TR / TS status
fec$Group = 'Herd'
fec$Group[fec$NUM %in% trueR]='TR'
fec$Group[fec$NUM %in% trueS]='TS'
fec$DATE_PRELEVEMENT = as.Date(fec$DATE_PRELEVEMENT,'%d/%m/%Y')
fec$Time='Before the experiment \n (2010 - 2015)'
fec$Time[fec$DATE_PRELEVEMENT>'2015-06-17' & fec$DATE_PRELEVEMENT<'2015-10-26']='During the experiment \n (June 2015 - October 2015)'
fec$Time[fec$DATE_PRELEVEMENT>'2015-10-26']='Following the experiment \n (2016 - 2020)'

table(fec$Group,fec$Time)
# Before the experiment \n (2010 - 2015) During the experiment \n (June 2015 - October 2015)
# Herd                                   2020                                                  74
# TR                                      155                                                  26
# TS                                      110                                                  24
# 
# Following the experiment \n (2016 - 2020)
# Herd                                      1881
# TR                                         183
# TS                                         177

### Plot and add+1 to FEC values for log-transformed axis representation

###----- Do the group differ before or after the experiment
#### Need to correct for age, season and treatment

##-- Age
age = read.csv(file = paste0(dir,'herd.csv'),header = T,sep =';')

##-- Treatment
trt = read.csv(file = paste0(dir,'drench.csv'),header = T,sep =';')
trt = trt[,c('NumAnimal','DateTrait','NomProduit')]
trt$drench = 'IVM'
trt$drench[trt$NomProduit=='PANACUR PATE']='FBZ'
trt$drench[trt$NomProduit=='STRONGID']='PYR'
trt$drench[grep('EQUEST',trt$NomProduit)]='MOX'

colnames(trt) = c('NUM','DrenchDate','Product','Molecule')
trt$DrenchDate = as.Date(trt$DrenchDate,'%d/%m/%Y')

trtold = read.csv(file = paste0(dir,'drench2.csv'),header = T,sep =';')
trtold$HEURE = NULL
trtold$Molecule = 'IVM'
trtold$Molecule[grep('EQUEST',trtold$NomProduit)]='MOX'
trtold$Molecule[grep('STRONGID',trtold$NomProduit)]='PYR'
trtold$Molecule[grep('PANACUR',trtold$NomProduit)]='FBZ'
trtold$Molecule[grep('CESTOCUR',trtold$NomProduit)]='PZQ'
trtold$Molecule[trtold$NomProduit=='YOMESANE']='NCL'
trtold$Molecule[trtold$NomProduit=='VERMEQUINE']='OXB'
trtold$Molecule[trtold$NomProduit=='MEFLOSYL']='Other'
trtold$DateTrait=as.Date(trtold$DateTrait,'%m/%d/%Y')
colnames(trtold) = c('NUM','DrenchDate','Product','Molecule')
trt = rbind(trtold[trtold$Molecule!='Other',],trt)
dim(trt)
# [1] 9579    4
table(trt$Molecule)
#  FBZ  IVM  MOX  NCL  OXB  PYR  PZQ 
# 1592 5912  566   58   57  912  482

### Retrieve individuals age
fec$birth = age$NAISSANCE[match(fec$NUM,age$NUM)]
fec$birth = as.Date(fec$birth,'%d/%m/%Y')
fec$age = (fec$DATE_PRELEVEMENT - fec$birth)/365

### Retrieve latest anthelmintic drug received
fetch_trt = function(x,d){
  tmp = trt[trt$NUM==x,]
  ## Retain anthelmintic drugs administered before weighing
  tmp = tmp[which((d-tmp$DrenchDate)>=0),]
  ## fetch FEC closest to weight
  i = which.min(abs(d-tmp$DrenchDate))
  ff = c(tmp$Molecule[i],as.character(tmp$DrenchDate[i]))
  return(ff)
}

for(i in 1:nrow(fec)){
  ff = fetch_trt(fec$NUM[i],fec$DATE_PRELEVEMENT[i])  
  fec$AH[i] = ff[1]
  fec$date_trt[i] = ff[2]
  rm(ff)
}

fec$AH = unlist(fec$AH)
fec = transform(fec, date_trt = as.Date(as.character(date_trt), "%Y-%m-%d"))
fec$difftrt = fec$DATE_PRELEVEMENT-fec$date_trt
summary(as.numeric(fec$difftrt))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       0      65     124     203     249    3090     244 

## Remove cases with missing info
fec = na.omit(fec)
dim(fec)
#[1] 2174   17
fec$age = as.numeric(fec$age)
fec$difftrt = as.numeric(fec$difftrt)

#      Before the experiment \n (2010 - 2015) During the experiment \n (June 2015 - October 2015)
# Herd                                    457                                                  47
# TR                                      100                                                  15
# TS                                       92                                                  21
# 
#       Following the experiment \n (2016 - 2020)
# Herd                                      1176
# TR                                         108
# TS                                         158

fecToPlot = fec
fecToPlot$Group = factor(fecToPlot$Group, levels = c('TR', 'Herd', 'TS'))

pdf(file='./Figure1.pdf')
ggplot(fecToPlot, aes(x = Group, y = RESULTATS_OPG, fill = Group)) +
  geom_point(aes(x = Group, y = RESULTATS_OPG,col = Group),
             position = position_jitterdodge(),pch=21, alpha = .4, size=2) +
  geom_boxplot(alpha = .7,outlier.shape=NA) + 
  facet_wrap(~ Time) + 
  scale_fill_manual(values = c(colsRS[1], 'grey',colsRS[2])) +
  scale_colour_manual(values = c(colsRS[1], 'grey',colsRS[2])) +
  ylab('Faecal Egg Count (eggs/g)') +
  xlab('') + theme_classic() +
  theme(legend.position = 'bottom', text = element_text(size = 12),strip.background = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 12))
dev.off()

#### Plot relationship between measured weight and FEC
fec$yr = format(fec$DATE_PRELEVEMENT,'%Y')

table(fec$yr)
# 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 
# 21   58  101  186  199  281  407  165  140  610    6 

fec = fec[fec$yr!=2020,] ## too few observations
fec$Time = factor(fec$Time)
fec$Tim = factor(match(fec$Time,levels(fec$Time)))

### Were TR, Herd and TS different before the experiment
m1 = lme(log(50 + RESULTATS_OPG) ~ age + difftrt + AH + yr + Group,
         random =~ 1|NUM,
         data = fec[fec$Tim==1,])
summary(m1)
# Linear mixed-effects model fit by REML
# Data: fec[fec$Tim == 1, ] 
# AIC  BIC logLik
# 2080 2146  -1025
# 
# Random effects:
#   Formula: ~1 | NUM
# (Intercept) Residual
# StdDev:      0.7317    1.061
# 
# Fixed effects: log(50 + RESULTATS_OPG) ~ age + difftrt + AH + yr + Group 
#              Value Std.Error  DF t-value p-value
# (Intercept)  5.555    1.1489 564   4.835  0.0000
# age         -0.160    0.0301 564  -5.318  0.0000
# difftrt      0.002    0.0003 564   6.856  0.0000
# AHIVM        0.685    1.1147 564   0.615  0.5390
# AHMOX       -0.322    1.1118 564  -0.290  0.7723
# AHPYR        0.650    1.1145 564   0.583  0.5601
# yr2011      -0.256    0.2753 564  -0.929  0.3531
# yr2012      -0.056    0.2760 564  -0.202  0.8401
# yr2013      -0.467    0.2719 564  -1.719  0.0862
# yr2014      -0.607    0.2775 564  -2.188  0.0291
# yr2015      -0.771    0.2960 564  -2.605  0.0094
# GroupTR     -0.862    0.3696  72  -2.333  0.0224
# GroupTS      0.826    0.3239  72   2.550  0.0129

### Were TR, Herd and TS different during the experiment
m2 = lme(log(50+RESULTATS_OPG) ~ age + difftrt + AH + Group, ## no year effect as all measures were done in 2015
         random =~ 1|NUM,
         data = fec[fec$Tim==2,])
summary(m2)
# Linear mixed-effects model fit by REML
# Data: fec[fec$Tim == 2, ] 
# AIC   BIC logLik
# 194.4 217.6  -87.2
# 
# Random effects:
#   Formula: ~1 | NUM
# (Intercept) Residual
# StdDev:     0.09923   0.6254
# 
# Fixed effects: log(50 + RESULTATS_OPG) ~ age + difftrt + AH + Group 
#               Value Std.Error DF t-value p-value
# (Intercept)  2.5103    0.7888 44   3.182  0.0027
# age         -0.0514    0.0239 44  -2.149  0.0372
# difftrt      0.0182    0.0038 44   4.728  0.0000
# AHIVM        0.4207    0.7709 44   0.546  0.5880
# AHMOX       -0.5908    0.6846 44  -0.863  0.3928
# AHPYR        0.9388    0.6789 44   1.383  0.1737
# GroupTR      0.1391    0.2130 31   0.653  0.5185
# GroupTS      0.6779    0.1868 31   3.629  0.0010

### Were TR, Herd and TS different before the experiment ?
m3 = lme(log(50+RESULTATS_OPG) ~ age + difftrt + AH + yr + Group,
         random =~ 1|NUM,
         data = fec[fec$Tim==3,])
summary(m3)
# Linear mixed-effects model fit by REML
# Data: fec[fec$Tim == 3, ] 
# AIC  BIC logLik
# 3891 3965  -1932
# 
# Random effects:
#   Formula: ~1 | NUM
# (Intercept) Residual
# StdDev:      0.8952    0.816
# 
# Fixed effects: log(50 + RESULTATS_OPG) ~ age + difftrt + AH + yr + Group 
#              Value Std.Error   DF t-value p-value
# (Intercept)  3.955    0.2575 1288  15.360  0.0000
# age         -0.114    0.0202 1288  -5.663  0.0000
# difftrt      0.002    0.0002 1288   8.876  0.0000
# AHIVM        1.764    0.2371 1288   7.441  0.0000
# AHMOX        1.105    0.2553 1288   4.329  0.0000
# AHPYR        1.619    0.2426 1288   6.673  0.0000
# yr2016       0.426    0.0903 1288   4.720  0.0000
# yr2017      -0.788    0.1367 1288  -5.768  0.0000
# yr2018      -0.027    0.1379 1288  -0.194  0.8463
# yr2019       0.106    0.1166 1288   0.913  0.3614
# GroupTR     -0.989    0.4261  136  -2.320  0.0218
# GroupTS      0.820    0.3598  136   2.278  0.0243

fec_time = aggregate(RESULTATS_OPG ~ Group, data = fec[fec$Tim==1,], FUN = summary)
fec_time$Time = 1
fec_time2 = aggregate(RESULTATS_OPG ~ Group, data = fec[fec$Tim==2,], FUN = summary)
fec_time2$Time=2
fec_time3 = aggregate(RESULTATS_OPG ~ Group, data = fec[fec$Tim==3,], FUN = summary)
fec_time3$Time=3

rbind(fec_time, fec_time2, fec_time3)
#   Group RESULTATS_OPG.Min. RESULTATS_OPG.1st Qu. RESULTATS_OPG.Median RESULTATS_OPG.Mean
# 1  Herd              0.000                 0.000              100.000            526.319
# 2    TR              0.000                 0.000                0.000            158.000
# 3    TS              0.000               100.000              575.000            733.696
# 4  Herd              0.000                 0.000                0.000             53.191
# 5    TR              0.000                 0.000                0.000              3.333
# 6    TS              0.000                 0.000                0.000            128.571
# 7  Herd              0.000                 0.000               50.000            320.595
# 8    TR              0.000                 0.000                0.000             43.056
# 9    TS              0.000               150.000              650.000            756.614
#   RESULTATS_OPG.3rd Qu. RESULTATS_OPG.Max. Time
# 1               650.000           6850.000    1
# 2                50.000           2200.000    1
# 3              1150.000           3200.000    1
# 4                 0.000            750.000    2
# 5                 0.000             50.000    2
# 6               200.000            700.000    2
# 7               350.000           5400.000    3
# 8                 0.000            900.000    3
# 9              1200.000           2800.000    3

#######=====================================================================================================
#######----------------------------- Can we distinguish R from S at D1 based on blood parameters - oPLS 
#######=====================================================================================================
sample = clinical$ID[clinical$Day=='D1']

clin = clinical[clinical$Day == 'D1',]
clin = clin[order(clin$ID,clin$Day),]

clin$FEC = NULL ## not to be considered as FEC negative individuals
clin$ADG = NULL ## 0 on day 1 so variance is 0

### 1st select top clinical phenotypes
## Y : keep the only variables showing differences between R and S at D5; top 10 of PLS-DA
Y = factor(clin$TG)
X = as.matrix(data.frame(clin[,-c(1:3)]))
rownames(X) = seq(1:dim(X)[1])
dim(X)
#[1] 17 19
length(Y)
#[1] 17

### Poor prediction and not significant
RS.d1.opls = ropls::opls(X, Y,predI = 1)
# PLS-DA
# 17 samples x 19 variables and 1 response
# standard scaling of predictors and response(s)
#       R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
# Total     0.29    0.304  0.0049 0.443   1   0  0.5 0.25

#######=====================================================================================================
#######----------------------------- Can we distinguish ponies to be treated at D132 based on blood parameters
#######=====================================================================================================

sample = clinical$ID[clinical$Day=='D5'] ##day 132 (5th sampling point)
sample = sample[sample !='W661'] ## remove outlier sample (detected in NMR)

clin = clinical[clinical$Day %in% c('D5'),]
clin = clin[-c(which(clin$ID=='W661' & clin$Day=='D5')),] ## remove outlier point
clin = clin[order(clin$ID,clin$Day),]

clin$FEC = NULL ## we want to predict w/o FEC

Y = factor(clin$TG)
X = as.matrix(clin[,-c(1:4)])
rownames(X) = seq(1:dim(X)[1])

dim(X)
#[1] 16 19

### Poorer prediction of the risk status // can't run prediction 
RS.d132.opls = ropls::opls(X, Y,predI = 1)
# PLS-DA
# 16 samples x 19 variables and 1 response
# standard scaling of predictors and response(s)
#      R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
#Total    0.161    0.437  -0.161 0.401   1   0  0.8 0.8

###=====================================================================================
###================------------ METABOLOMICS data analysis 
###=====================================================================================

###--------============== 1. Age correction for NMR data?
### Check this at D1 only to remove age-biased metabolites
age.vec2 = clinical[clinical$Day %in% c('D1'),c('ID','Day','Age')]
age.vec2 = age.vec2$Age[order(age.vec2$ID,age.vec2$Day)]
age.vec2

bc1 = bc.dnoised[meta$Day=='J0',]

mod.stor2 = NULL

for(i in seq(1:dim(bc1)[2])){
  y = bc1[,i]
  p = data.frame(summary(lm(y ~ age.vec2))$coefficients)
  p$param = colnames(bc1)[i]
  mod.stor2 = rbind(mod.stor2,p[2,])
}
p = mod.stor2$Pr...t..
mod.stor2$fdr = p.adjust(p, method = 'BH', n = length(p))

mod.stor2[mod.stor2$fdr < .05,]
# [1] Estimate   Std..Error t.value    Pr...t..   param      fdr
# <0 rows> (or 0-length row.names)

###--------============== 2. First pass on metabolomics data highlights major shifts due to pasture season
### No differences between TR and TS throughout

sample = as.data.frame(meta$Ind)

PCA.denoised <- mixOmics::pca(bc.dnoised,
                              center = T, scale = T,ncomp = 5, multilevel = sample)
fviz_pca_ind(PCA.denoised,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level=0.95,
             legend.title = "Groups")

fviz_pca_ind(PCA.denoised,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = meta$Day, # color by groups
             #palette = colsD[1:3],
             addEllipses = TRUE, # Concentration ellipses
             ellipse.level=0.95,
             legend.title = "Groups")

### Check Normality
summary(sapply(sapply(lapply(data.frame(bc.dnoised),FUN = shapiro.test),function(x) x[1]),function(x) x[1]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.818   0.952   0.970   0.962   0.981   0.991   

#######----------------------------- ASCA R vs. S throughout study - metaboanalyst
detach("package:mixOmics", unload=TRUE)
require(MetaboAnalystR)

nmr = t(bc.dnoised)
Sample = meta$Ind
Phenotype = as.factor(meta$TG)
Time = meta$Day

nmrSet = rbind(Sample, Phenotype, Time, nmr)
nmrSet = data.frame(nmrSet)
write.table(nmrSet,file = './data/nmr_metabo.csv',sep=",",  quote = F,col.names=FALSE)

mSet<-InitDataObjects("pktable", "ts", FALSE)
mSet<-SetDesignType(mSet, "time")
mSet<-Read.TextData(mSet, "./data/nmr_metabo.csv", "colts", "disc")
mSet<-SanityCheckData(mSet)
# [1] "Successfully passed sanity check!"                                                                                
# [2] "The data is time-series data."                                                                                    
# [3] "2 groups were detected in samples for factor Phenotype"                                                           
# [4] "3 groups were detected in samples for factor Time"                                                                
# [5] "Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed."                             
# [6] "<font color=\"orange\">Other special characters or punctuations (if any) will be stripped off.</font>"            
# [7] "All data values are numeric."                                                                                     
# [8] "A total of 0 (0%) missing values were detected."                                                                  
# [9] "<u>By default, missing values will be replaced by 1/5 of min positive values of their corresponding variables</u>"
# [10] "Click the <b>Skip</b> button if you accept the default practice;"                                                 
# [11] "Or click the <b>Missing value imputation</b> to use other methods."    

# Perform data processing - Minimum Value Replacing
mSet<-ReplaceMin(mSet);

# Perform data processing - Variable Filtering and Normalization
mSet<-PreparePrenormData(mSet);
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", "S10T0", ratio=FALSE, ratioNum=20);

mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA);
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA);

## Perform ASCA, specifying model components
## implementation of ASCA in MetaboAnalyst was based on the
## algorithm described by AK Smildle, et al. 2005 (PMID: 15890747)

# a, specify the number of components for facA
# b, specify the number of components for facB
# x, specify the number of components for interaction AB
# res, specify the number of model residuals type is string, 
# indicating the type of analysis "abc" separately 
# "aab" facA joins with AB 
# "bab" facB joins with AB

mSet<-Perform.ASCA(mSet, a=1, b=1, x=2, res=2)

## Create scree plots of each model
mSet<-PlotModelScree(mSet, "asca_scree_0_", format = "png", dpi = 72, width=NA)
## Plot ASCA model A = Group
mSet<-PlotASCAModel(mSet, "asca_fa_0_", format = "png", dpi = 72, width=NA, "a",FALSE)
## Plot ASCA model B = Time
mSet<-PlotASCAModel(mSet, "asca_fb_0_", format = "png", dpi = 72, width=NA, "b",FALSE)

# Plot ASCA Interaction
mSet<-PlotInteraction(mSet, "asca_fab_0_", "png", 72,FALSE, width=NA)
# Perform model validation, 1000 permutations
mSet<-Perform.ASCA.permute(mSet, 1000)
# Plot model validation
mSet<-PlotASCA.Permutation(mSet, "asca_perm_0_", "png", 72, width=NA)

# Calculate significant features, specifying the alpha threshold (spe.thresh) and leverage threshold (lev.thresh)
mSet<-CalculateImpVarCutoff(mSet, spe.thresh = 0.05, lev.thresh = 0.9)

# Plots of significant features for each model
mSet<-PlotAscaImpVar(mSet, "asca_impa_0_", "png", 72, width=NA, "a")
mSet<-PlotAscaImpVar(mSet, "asca_impb_0_", "png", 72, width=NA, "b")
mSet<-PlotAscaImpVar(mSet, "asca_impab_0_", "png", 72, width=NA, "ab")

### Results
mSet$analSet$asca$perm.p ## no group effect, significant day effect, no group x day interaction
#[1] "p = 0.344 (344/1000)" "p < 0.001 (1/1000)"   "p = 0.544 (544/1000)"

### Identify outliers to not consider these further
# higher Leverage means more importance) 
# SPE (or squared prediction error, 
# an evaluation of the goodness of fit of the model to a particular variable
# higher SPE means less fit 

## Group effect
outlier.a = data.frame(mSet$analSet$asca$out.list$Model.a)
outlier.a
#                            Leverage       SPE
# Methylene.Moieties.From.Fa  0.29141 1.643e-07
# Methyl.Moieties.From.Fa     0.26704 1.506e-07
# Glucose.2                   0.26459 1.492e-07
# Choline                     0.06224 3.509e-08
# Glucose.1                   0.05633 3.176e-08

## Time effect
outlier.b = data.frame(mSet$analSet$asca$out.list$Model.b)
outlier.b
#                            Leverage       SPE
# Methylene.Moieties.From.Fa 0.334726 8.130e-04
# Methyl.Moieties.From.Fa    0.079385 6.442e-04
# Choline                    0.051071 1.661e-04
# Lactate.Threonine.Lipid    0.000426 1.316e-04
# Acetate                    0.003164 1.236e-04
# Glucose.2                  0.280021 9.458e-05

outlier.ab = data.frame(mSet$analSet$asca$out.list$Model.ab)
outlier.ab
#                  Leverage       SPE
# Glucose.2       0.0173803 7.050e-34
# Glucose.Proline 0.0006901 6.533e-34
# Glucose.1       0.0013252 5.465e-34
# Acetamide       0.0091719 3.793e-34
# Choline         0.2011916 1.651e-34

outlier.features = unique(c(rownames(outlier.a),rownames(outlier.b),rownames(outlier.ab)))
outlier.features
# [1] "Methylene.Moieties.From.Fa" "Methyl.Moieties.From.Fa"    "Glucose.2"                 
# [4] "Choline"                    "Glucose.1"                  "Lactate.Threonine.Lipid"   
# [7] "Acetate"                    "Glucose.Proline"            "Acetamide"

## Plot significant metabolites associated with Time
## w/o outliers

sig.b = data.frame(mSet$analSet$asca$sig.list$Model.b)
sig.b = sig.b[!(seq(1:dim(sig.b)[1]) %in% (which(rownames(sig.b) %in% outlier.b))),]
sig.b = sig.b[order(rownames(sig.b)),]
sig.b
#                           Leverage       SPE
# Alkalene.Moieties.From.Fa  0.02087 1.524e-05
# Glucose.1                  0.07826 2.220e-05
# Glucose.Proline            0.03685 1.745e-05
# Leucine                    0.01488 5.159e-06
# Valine.1                   0.01893 5.638e-07

summary(sig.b$Leverage)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.0149  0.0189  0.0209  0.0340  0.0369  0.0783

summary(sig.b$SPE)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 5.64e-07 1.16e-05 1.63e-05 1.39e-05 1.86e-05 2.22e-05 

## Match NMR data with sig.b for plotting
nmr.dat = data.frame(t(nmrSet))

nmr.b = melt(nmr.dat[,c(1:3,match(rownames(sig.b),colnames(nmr.dat)))],1:3)
nmr.b$value = as.numeric(as.character(nmr.b$value))
nmr.b$Shift = nmr_peak$meta_peak[match(nmr.b$variable,colnames(bc.dnoised))]
nmr.b$variable2 = paste0(nmr.b$variable,' - ',nmr.b$Shift,' ppm')

## Figure showing differential metabolites through time
pdf(file = './supplementary_Figure1.pdf',width=14,height=8)
ggplot(nmr.b,aes(x = Time, y = value)) +
  geom_boxplot()+
  facet_wrap(~variable2,scales ='free') + 
  geom_point(aes(col=Phenotype),position = position_jitterdodge())+
  scale_color_manual(values = colsRS)+ theme_classic() +
  scale_x_discrete(labels = c('Day 0','Day 24','Day 132')) +
  xlab('') + ylab('Peak intensity')+coord_flip()+
  theme(legend.position = 'none', 
        strip.background = element_blank(),text = element_text(size = 12)) 
invisible(dev.off())

### Update bc.dnoised2 to bc.dnoised.cln
bc.dnoised.cln = bc.dnoised[,-c(match(outlier.features,colnames(bc.dnoised)))]
dim(bc.dnoised.cln)
#[1]  51 110

####----------------------- PLS-DA on day effect
Y = meta$Day
X = bc.dnoised.cln
rownames(X) = seq(1:dim(X)[1])

dim(X)
#[1]  51 110

### Model
Day.opls = ropls::opls(X, Y)
# PLS-DA
# 51 samples x 110 variables and 1 response
# standard scaling of predictors and response(s)
# R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
# Total    0.683    0.886   0.836 0.166   3   0 0.05 0.05

###=====================================================================================
###= OPLS analyses: R vs. S at day0 
###=====================================================================================

####---- Can we pick metabolites that differ betw. R and S ponies at day 0 ?
Y = meta$TG[meta$Day=='J0']
X = bc.dnoised.cln[meta$Day=='J0',]
rownames(X) = seq(1:dim(X)[1])

dim(X)
#[1] 17 110

### Poorer prediction of the risk status 
RS.Day0.opls = ropls::opls(X, Y,predI = 1)
# PLS-DA
# 17 samples x 110 variables and 1 response
# standard scaling of predictors and response(s)
# R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
# Total    0.292    0.366 -0.0286 0.423   1   0 0.75 0.35

### Metabolites with most significant differences between R and S before pasture turnout
v = NULL
for(i in 1:ncol(X)){
  v[i] = t.test(X[,i] ~ Y)$p.value
}
fdr = p.adjust(v, method = 'BH', n = length(v))
min(fdr)
#[1] 0.7189

## Nominal p-values
colnames(X[,which(v<.05)])
#[1] "Lysine.2"           "Lysine.3"           "Dimethyl.Sulfone.1"

df.nmr.0 = melt(data.frame(stat = Y, X[,which(v<.05)]),1)
levels(df.nmr.0$variable)
#[1] "Lysine.2"           "Lysine.3"           "Dimethyl.Sulfone.1"

###--- Is there any relationship between differential metabolites at day 0 and final FEC ?

### Correlation with final FEC
FEC = clinical[clinical$Day=='D5','FEC']
df.nmrs.05 = data.frame(FEC = FEC, day = rep(132,17), bc.dnoised.cln[meta$Day=='RB',which(v<.05)])
df.nmrs.05$day = factor(df.nmrs.05$day)
df.nmrs.05 = df.nmrs.05[df.nmrs.05$day==132,]
df.nmrs.05$day = NULL
cormat = rcorr(as.matrix(df.nmrs.05),type = 'spearman')$r

corrplot::corrplot(cormat,type = 'upper',#insig = 'blank',sig.level = .05,
                   col = viridis_pal(option = 'D')(100))

rcorr(as.matrix(df.nmrs.05),type = 'spearman')$r[1,]
#     FEC           Lysine.2           Lysine.3 Dimethyl.Sulfone.1 
# 1.00000           -0.08871           -0.07997            0.34110 
#not significant
rcorr(as.matrix(df.nmrs.05),type = 'spearman')$P[1,][rcorr(as.matrix(df.nmrs.05),type = 'spearman')$P[1,]<.05]
# <NA>           
# NA             

### Plot
dfplot.corFEC=melt(df.nmrs.05,1)
df.nmr.0$variable2 = paste0(df.nmr.0$variable,' \n ',
                             nmr_peak$meta_peak[match(df.nmr.0$variable,colnames(bc.dnoised))], ' ppm')
dfplot.corFEC$variable2 = paste0(dfplot.corFEC$variable,' \n ',
                              nmr_peak$meta_peak[match(dfplot.corFEC$variable,colnames(bc.dnoised))], ' ppm')

FigS2A = ggplot(df.nmr.0,aes(x = stat, y = value,
                            group = stat)) + 
  facet_wrap(~ variable2, scales = 'free',ncol = 3) +
  geom_boxplot(aes(fill = stat, group = stat,alpha = .5),outlier.color = NA) + 
  geom_point(aes(col = stat,group = stat), position = position_jitterdodge()) +
  scale_color_manual(values = c(viridis_pal(option="D")(6)[c(2,5)])) +
  scale_fill_manual(values = c(viridis_pal(option="D")(6)[c(2,5)])) +
  theme_classic() +
  xlab('Pony resistance status') + ylab('Metabolite peak intensity') +
  ggtitle('a') +
  theme(legend.position = 'none',
        strip.background = element_blank(),text = element_text(size = 10))

FigS2B = ggplot(dfplot.corFEC,aes(y = value, x = FEC)) +
  geom_smooth(method = 'lm',col = 'black', lty = 2,lwd = .5,alpha =.2,
              fill = c('skyblue')) +
  geom_point(size = 3, alpha = 0.5) +
  facet_wrap(~ variable2, scales = "free",ncol=3) + theme_classic() +
  xlab('Log-transformed Faecal Egg Count') + #coord_flip() +
  ylab('Metabolite peak intensity') + ggtitle('b')+
  theme(legend.position = 'none',
        strip.background = element_blank(),
       text = element_text(size = 10))
 
pdf(file = './Supplementary_Figure2.pdf') #,width = 14, height = 8)
multiplot(FigS2A,FigS2B, cols = 1)
dev.off()

###--- Metabolites that differ between R and S at day 132
Y1 = factor(meta$TG[meta$Day=='RB'])
X1 = bc.dnoised.cln[which(meta$Day=='RB'),]
rownames(X1) = seq(1:dim(X1)[1])

dim(X1)
#[1] 17 110

### Poorer prediction of the risk status 
RS.Day5.opls = ropls::opls(X1, Y1,predI = 1)
# PLS-DA
# 17 samples x 110 variables and 1 response
# standard scaling of predictors and response(s)
#       R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y pQ2
# Total     0.42    0.239  0.0282 0.463   1   0 0.75 0.3

### Metabolites with most significant differences between R and S at pasture turnout
v1=NULL
for(i in 1:ncol(X1)){
  v1[i] = t.test(X1[,i] ~ Y1)$p.value
}
fdr1 = p.adjust(v1, method = 'BH', n = length(v1))
min(fdr1)
#[1] 0.361

## Nominal p-values
colnames(X1[,which(v1<.05)])
#[1] "Phenylalanine.2" "U13"             "Phenylalanine.3"  

df.nmr.5 = melt(data.frame(stat = Y1, X1[,which(v1<.05)]),1)
levels(df.nmr.5$variable)
#[1] "Phenylalanine.2" "U13"             "Phenylalanine.3"  

df.nmr.5$variable2 = paste0(df.nmr.5$variable,' \n ',
                            nmr_peak$Chemical_shift[match(df.nmr.5$variable,colnames(bc.dnoised))], ' ppm')

### Correlation with final FEC
FEC = clinical[clinical$Day=='D5','FEC']
df.nmrs.55 = data.frame(FEC = FEC, X1[,which(v1<.05)])

cormat = rcorr(as.matrix(df.nmrs.55),type = 'spearman')$r
corrplot::corrplot(cormat,type = 'upper',#insig = 'blank',sig.level = .05,
                   col = viridis_pal(option = 'D')(100))

rcorr(as.matrix(df.nmrs.55),type = 'spearman')$r[1,]
#    FEC Phenylalanine.2             U13 Phenylalanine.3 
# 1.0000          0.5560          0.3911          0.6060 

rcorr(as.matrix(df.nmrs.55),
      type = 'spearman')$P[1,][rcorr(as.matrix(df.nmrs.55),
                                     type = 'spearman')$P[1,]<.05][-1]
# Phenylalanine.2 Phenylalanine.3 
#        0.020471        0.009925 

corFEC2 = names(rcorr(as.matrix(df.nmrs.55),type = 'spearman')$P[1,])
corFEC2 = corFEC2[-1]
dfplot.corFEC2 = melt(df.nmrs.55[,c('FEC',corFEC2)],1)

dfplot.corFEC2$variable2 = paste0(dfplot.corFEC2$variable, ' \n ',
                                  nmr_peak$Chemical_shift[match(dfplot.corFEC2$variable,colnames(bc.dnoised))],
                                  ' ppm')

Fig3A = ggplot(df.nmr.5[df.nmr.5$variable!='U13',],aes(x = stat, y = value,
                             group = stat)) + 
  facet_wrap(~ variable2, scales = 'free',ncol = 3) +
  geom_boxplot(aes(fill = stat, group = stat,alpha = .5),outlier.color = NA) + 
  geom_point(aes(col = stat,group = stat), position = position_jitterdodge()) +
  scale_color_manual(values = c(viridis_pal(option="D")(6)[c(2,5)])) +
  scale_fill_manual(values = c(viridis_pal(option="D")(6)[c(2,5)])) +
  theme_classic() +
  xlab('Pony resistance status') + ylab('Metabolite peak intensity') +
  ggtitle('A') +
  theme(legend.position = 'none',
        strip.background = element_blank(),text = element_text(size = 10))

Fig3B = ggplot(dfplot.corFEC2[df.nmr.5$variable!='U13',],aes(y = value, x = FEC)) +
  geom_smooth(method = 'lm',col = 'black', lty = 2,lwd = .5,alpha =.2,
              fill = c('skyblue')) + 
  geom_point(size = 3, alpha = 0.5) + 
  facet_wrap(~ variable2, scales = "free", ncol = 3) + theme_classic() + 
  xlab('Log-transformed Faecal Egg Count') + #coord_flip() +
  ylab('Metabolite peak intensity') + ggtitle('B')+
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 10)) 

#####---------------=========== Run validation on Peachey et al. 2019 Sci rep dataset
metavalid = read.csv(file = './PeacheyData/metadata.v6.csv', header = T, sep =',')
nmrvalid = read.csv(file = './PeacheyData/Metabolomicdata_correctedforDM.txt', sep='\t')
colnames(nmrvalid)[2:dim(nmrvalid)[2]] = substr(colnames(nmrvalid)[2:dim(nmrvalid)[2]], 2, nchar(colnames(nmrvalid)[2:dim(nmrvalid)[2]]))
rownames(nmrvalid) = nmrvalid$OTU.ID
nmrvalid = t(nmrvalid[,-1])

dim(nmrvalid)
#[1] 58 28

## Keep samples of interest
metavalid = metavalid[(metavalid$Infection==0 & metavalid$Timepoint=='a')|
                        (metavalid$Infection==1 & metavalid$Timepoint=='a'),]
nmrvalid = nmrvalid[!(is.na(match(rownames(nmrvalid),metavalid$Subject))),]
dim(nmrvalid)
#[1] 23 28

## Try to match subject across datasets
sapply(paste0("_",meta$ID),function(x) grep(x,metavalid$Subject))
metavalid$Subject = as.character(metavalid$Subject)

## issues: 2 ML but no HL
## inspection of infection x timepoint suggests that HL 0 a = 47_10_ML1 / ML1 = 3_108_ML1
metavalid$Subject[metavalid$Subject=='47_10_ML1']='47_10_HL1'
ord = sapply(paste0("_",meta$ID),function(x) grep(x,metavalid$Subject))
ord

## Reorder NMR data to match 16S dataset
metavalid = metavalid[ord,]
nmrvalid = nmrvalid[ord,]
nmrvalid = data.frame(nmrvalid)

nmvalid = colnames(nmrvalid)[as.integer(na.omit(match(unique(unlist(stringr::str_split(colnames(X1[,which(v1<.05)]),'\\.'))),colnames(nmrvalid))))]
nmvalid
#[1] "Phenylalanine"

## Test whether Phe is differential
wilcox.test(nmrvalid[,'Phenylalanine'] ~ factor(metavalid$Infection)) ##
# Wilcoxon rank sum test
# 
# data:  nmrvalid[, "Phenylalanine"] by factor(metavalid$Infection)
# W = 27, p-value = 0.02
# alternative hypothesis: true location shift is not equal to 0

vn = cbind(metavalid[,c('Subject','Infection')],
           nmrvalid[,nmvalid])
colnames(vn)[3] = nmvalid
vn$Infection = factor(vn$Infection)
v_nmr_melt = reshape2::melt(vn,c('Subject','Infection'))
v_nmr_melt$stat='Low'
v_nmr_melt$stat[v_nmr_melt$Infection==1]='High'
v_nmr_melt$stat=factor(v_nmr_melt$stat,levels=c('Low','High'))

Fig3C = ggplot(v_nmr_melt,aes(x = stat, y = value,fill=Infection)) +
  theme_classic() + 
  geom_boxplot(alpha=.8) + scale_fill_manual(values = colsStat) + xlab('') +
  theme(legend.position = 'none') + ylab('Phenylalanine \n peak intensity') +
  ggtitle("C") + scale_y_continuous(limits = c(0.01, 0.45), breaks = seq(0.01, 0.45,0.05)) +
  theme(legend.position = 'none',
        text = element_text(size = 12))

pdf(file='./Figure3.pdf')
multiplot(Fig3A,Fig3C, Fig3B,cols=2)
dev.off()


###=====================================================================================
###===== COVARIATION data analysis - DIABLO between R and S before pasture
###=====================================================================================

###-- Reorder every dataset as NMR profiles (by ind Day)
metatot = metaotu[metaotu$Time == 0,] ## Metadata
metatot = metatot[order(metatot$ID,metatot$Time),]

clin = clinical[clinical$ID %in% meta$Ind & clinical$Day %in% c('D1'),]
clin = clin[order(clin$ID,clin$Day),]
stat = clin$TG
clin = clin[,-c(1:4)]
clin$FEC = NULL ## remove FEC data 
clin$ADG = NULL ## no gain

#####--------- Prep ASVs data

## Filter Genus data
# keep day0
genus0 = genus[rownames(genus) %in% metatot$SampleID_Sequencing,] 
result.filter = low.count.removal(genus0, percent = 0.05)
genus.filter = result.filter$data.filter
length(result.filter$keep.otu) #
##[1] 44 ASVs
genuscov0 = genus.filter

# function is applied to each row (i.e. each sample)
genuscov = mixOmics::logratio.transfo(result.filter$data.filter, 
                                      logratio = 'CLR',offset = 1)
dim(genuscov)
#[1] 17 44

summary(unlist(sapply(apply(genuscov,2, shapiro.test),function(x) x[1])))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.694   0.851   0.903   0.890   0.946   0.994 

### NMR data
nmrcov = bc.dnoised.cln[meta$Day=='J0',]

### rownames
rownames(clin) = seq(1:dim(clin)[1])
rownames(genuscov) = seq(1:dim(genuscov)[1])
rownames(nmrcov) = seq(1:dim(nmrcov)[1])
colnames(nmrcov) = gsub("X",'',colnames(nmrcov))

## Rename [Prevotella]
colnames(genuscov)=gsub('X.Prevotella.','Prevotella.1',colnames(genuscov))

### X matrix
data = list(genus = genuscov,
            nmr = nmrcov,
            clin = clin)

# check dimension
lapply(data, dim)
# $genus
# [1] 17 44
# 
# $nmr
# [1]  17 110
# 
# $clin
# [1] 17 18

#####------------------- D0 - DIABLO
detach('package:MetaboAnalystR',unload = T)
detach('package:ade4',unload = T)
detach('package:amap',unload = T)

require(mixOmics)
### Get insight of correlations between datasets
tot.pls1 <- pls(genuscov, clin, ncomp = 1, mode = "regression")
cor(tot.pls1$variates$X, tot.pls1$variates$Y)
#       comp1
# comp1 0.844

tot.pls2 <- pls(genuscov, nmrcov, ncomp = 1, mode = "regression")
cor(tot.pls2$variates$X, tot.pls2$variates$Y)
#       comp1
# comp1 0.8703

tot.pls3 <- pls(nmrcov, clin, ncomp = 1, mode = "regression")
cor(tot.pls3$variates$X, tot.pls3$variates$Y)
#         comp 1
# comp 1   0.8013

### Parameter choice
# Indicate what parameter should be connected
# A compromise needs to be achieved between maximising the correlation between data sets (design value between 0.5 and 1) 
# and maximising the discrimination with the outcome Y (design value between 0 and 0.5)
design1 = matrix(c(0, 0.8703, 0.844,
                   0.8703, 0, 0.8013,
                   0.844, 0.8013, 0), 
                 byrow = T,
                 ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
design1 
#        genus    nmr   clin
# genus 0.0000 0.8703 0.8440
# nmr   0.8703 0.0000 0.8013
# clin  0.8440 0.8013 0.0000

### Use block.plsda to infer relationship between K's
nc = 2 ## K - 1

sgccda.res1 = block.plsda(X = data, Y = stat, design = design1, 
                          near.zero.var = T,
                          ncomp = nc)

perf.1 = perf(sgccda.res1, validation = 'Mfold', folds = 5, 
              nrepeat = 10)

perf.1$error.rate
# $genus
#       max.dist centroids.dist mahalanobis.dist
# comp1   0.4647         0.4941           0.4941
# comp2   0.3353         0.3353           0.3176
# 
# $nmr
#       max.dist centroids.dist mahalanobis.dist
# comp1   0.6235         0.6059           0.6059
# comp2   0.4588         0.4647           0.4529
# 
# $clin
#       max.dist centroids.dist mahalanobis.dist
# comp1   0.5353         0.6000           0.6000
# comp2   0.4412         0.4235           0.4412

#### Optimize number of features with updated design
ncomp = perf.1$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp 
#[1] 2

## Plot matrix relationship on component 1
plotIndiv(sgccda.res1,pch = stat)
plotDiablo(sgccda.res1, ncomp = 1) ## R vs. S
#plotDiablo(sgccda.res1, ncomp = 2)

## Update design
design2 = matrix(c(0, 0.87, 0.68,
                   0.87, 0, 0.83,
                   0.68, 0.83, 0), 
                 byrow = T,
                 ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
design2
#       genus  nmr clin
# genus  0.00 0.87 0.68
# nmr    0.87 0.00 0.83
# clin   0.68 0.83 0.00

## Optimization for the number of features to be retained in each block
set.seed(123) # for reproducibility, only when the `cpus' argument is not used
test.keepX = list (genus = seq(10,20,5),
                   nmr = seq(10,80,10),
                   clin = seq(5,10,1))

tune.TCGA = tune.block.splsda(X = data, Y = stat, ncomp = 2, 
                              test.keepX = test.keepX, 
                              near.zero.var = T,design = design2,
                              validation = 'Mfold', folds = 5, nrepeat = 10,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX
list.keepX
# $genus
# [1] 20 10
# 
# $nmr
# [1] 10 30
# 
# $clin
# [1] 6 6

### Final model with appropriate design and number of features
sgccda.rs0 = block.splsda(X = data, Y = stat, ncomp = 2, 
                          near.zero.var = T,
                          keepX = list.keepX, design = design2)
sgccda.rs0$design
#       genus  nmr clin Y
# genus  0.00 0.87 0.68 1
# nmr    0.87 0.00 0.83 1
# clin   0.68 0.83 0.00 1
# Y      1.00 1.00 1.00 0

### Selected variables // status 
gn0 = factor(sort(selectVar(sgccda.rs0, block = 'genus', comp = 1)$genus$name))
gn0
# [1] Adlercreutzia     Anaerofustis      BF311             CF231             Clostridium.48   
# [6] Coprococcus       Desulfovibrio     Dietzia           Fibrobacter       Gordonia         
# [11] Mycobacterium     Paludibacter      Prevotella        Pseudomonas       Ruminococcus     
# [16] Saccharopolyspora Solibacillus      Sphingomonas      Treponema         YRC22  

nm0 = factor(sort(selectVar(sgccda.rs0, block = 'nmr', comp = 1)$nmr$name))
nm0
# [1] Lysine                     Lysine.3                   Threonine.1               
# [4] U13                        U18                        U21                       
# [7] U4                         U5                         X1.Methylhistidine.1      
# [10] X3.Hydroxybutyrate.Proline

## Shifts for methylhistidine, lysine
nmr_peak$Chemical_shift[grep('1.Methylhistidine',nmr_peak$Metabolite_Reannotated)][2]
#[1] "7.785-7.785"
nmr_peak$Chemical_shift[grep('Lysine',nmr_peak$Metabolite_Reannotated)][1]
#[1] "1.445-1.465"
nmr_peak$Chemical_shift[grep('Lysine',nmr_peak$Metabolite_Reannotated)][4]
#[1] "1.845-1.915"
nmr_peak$Chemical_shift[grep('3.Hydroxybutyrate.Proline',nmr_peak$Metabolite_Reannotated)][1]
#[1] "4.125-4.135"

cl0 = factor(sort(selectVar(sgccda.rs0, block = 'clin', comp = 1)$clin$name))
cl0
#[1] Eosinophils Glucose     Leukocytes  Lymphocytes MGV         Neutrophil 

plotDiablo(sgccda.rs0, ncomp = 1, col = colsRS) ## R vs. S
plotIndiv(sgccda.rs0, legend = TRUE, title = 'DIABLO',pch = stat)
plotVar(sgccda.rs0, var.names = TRUE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(.8,.8,.8))

### Loadings on component 1, associated with R / S differences
pdf(file = './Supplementary_Figure6.pdf', width = 14, height = 8)
plotLoadings(sgccda.rs0, comp = 1, contrib = 'max', #method = 'median',
             size.name = 1.1,legend = F,legend.color = colsRS,
             subtitle = c('Component 1 \n Bacteria genera',
                          'Component 1 \n Metabolites',
                          'Component 1 \n Clinical features'))
dev.off()

#### Circos plot 
pdf(file = './Supplementary_Figure5.pdf', width = 14, height = 14)
circp = circosPlot(sgccda.rs0, comp = 1, 
           cutoff = 0.5, line = TRUE, size.variables = 1.7,
           color.blocks = viridis_pal(option='D')(10)[c(3,7,10)],
           color.Y = colsRS,legend = F,
           color.cor = c("chocolate3","grey20"), size.labels = 0.0000001)
dev.off()

#### Correlation between the parameters of interest
gn0 = as.character(gn0)
nm0 = as.character(nm0)
cl0 = as.character(cl0)

feat_mat0 = circp[c(gn0,nm0,cl0),c(gn0,nm0,cl0)]

pdf(file ='./Supplementary_Figure7.pdf')
print(pheatmap::pheatmap(feat_mat0,fontsize = 8,
                   cluster_rows = T,cluster_cols = T))
dev.off()

### Strongest association for Immune cells on component 1
factor(names(feat_mat0[,'Lymphocytes'][abs(feat_mat0[,'Lymphocytes'])  > 0.5]))
# [1] Anaerofustis  Coprococcus   Dietzia       Fibrobacter   Gordonia      Mycobacterium
# [7] Ruminococcus  Solibacillus  Treponema     U4            U5            Leukocytes   
# [13] Lymphocytes 

### Strongest association for Immune cells on component 1
factor(names(feat_mat0[,'Neutrophil'][abs(feat_mat0[,'Neutrophil'])  > 0.5]))
#[1] Coprococcus  Ruminococcus Treponema    U4           Leukocytes  

feat_mat0[c('Coprococcus','Ruminococcus','Anaerofustis','Neutrophil','Lymphocytes'),
          c('Coprococcus','Ruminococcus','Anaerofustis','Neutrophil','Lymphocytes')]
#              Coprococcus Ruminococcus Anaerofustis Neutrophil Lymphocytes
# Coprococcus       0.8303       0.8489       0.6288     0.5030      0.6697
# Ruminococcus      0.8489       0.8679       0.6429     0.5143      0.6847
# Anaerofustis      0.6288       0.6429       0.4763     0.3810      0.5072
# Neutrophil        0.5030       0.5143       0.3810     0.3047      0.4057
# Lymphocytes       0.6697       0.6847       0.5072     0.4057      0.5402

summary(feat_mat0[c('Coprococcus','Ruminococcus','Anaerofustis','Lymphocytes'),
                  c('Coprococcus','Ruminococcus','Anaerofustis','Lymphocytes')])
#   Coprococcus     Ruminococcus    Anaerofustis    Lymphocytes   
# Min.   :0.629   Min.   :0.643   Min.   :0.476   Min.   :0.507  
# 1st Qu.:0.659   1st Qu.:0.674   1st Qu.:0.499   1st Qu.:0.532  
# Median :0.750   Median :0.767   Median :0.568   Median :0.605  
# Mean   :0.744   Mean   :0.761   Mean   :0.564   Mean   :0.601  
# 3rd Qu.:0.835   3rd Qu.:0.854   3rd Qu.:0.632   3rd Qu.:0.673  
# Max.   :0.849   Max.   :0.868   Max.   :0.643   Max.   :0.685 

summary(feat_mat0[c('Coprococcus','Ruminococcus','Anaerofustis','Neutrophil'),
                  c('Coprococcus','Ruminococcus','Anaerofustis','Neutrophil')])
#   Coprococcus     Ruminococcus    Anaerofustis     Neutrophil   
# Min.   :0.503   Min.   :0.514   Min.   :0.381   Min.   :0.305  
# 1st Qu.:0.597   1st Qu.:0.611   1st Qu.:0.453   1st Qu.:0.362  
# Median :0.730   Median :0.746   Median :0.553   Median :0.442  
# Mean   :0.703   Mean   :0.719   Mean   :0.532   Mean   :0.426  
# 3rd Qu.:0.835   3rd Qu.:0.854   3rd Qu.:0.632   3rd Qu.:0.506  
# Max.   :0.849   Max.   :0.868   Max.   :0.643   Max.   :0.514

### Bacteria vs. TR status
sort(sgccda.rs0$loadings$genus[,1][sgccda.rs0$loadings$genus[,1]>0],decreasing=T)
# Dietzia     Mycobacterium      Solibacillus          Gordonia Saccharopolyspora             BF311     Desulfovibrio 
# 0.42273           0.30776           0.26584           0.23855           0.07966           0.03405           0.02793 

sort(sgccda.rs0$loadings$nmr[,1][sgccda.rs0$loadings$nmr[,1]>0],decreasing=T)
#     U21                    Lysine               Threonine.1 
# 0.29058                   0.28117                   0.13189 
# Lysine.3                       U18                       U13 
#  0.07124                   0.04893                   0.04635 
# 3.Hydroxybutyrate.Proline 
#                   0.03365

sort(sgccda.rs0$loadings$clin[,1][sgccda.rs0$loadings$clin[,1]<0],decreasing=T)

### Bacteria vs. TS status
sort(sgccda.rs0$loadings$genus[,1][sgccda.rs0$loadings$genus[,1]<0])
# Treponema   Ruminococcus    Coprococcus    Fibrobacter   Paludibacter   Anaerofustis Clostridium.48  Adlercreutzia          CF231 
# -0.370528      -0.355255      -0.316883      -0.308908      -0.221248      -0.157632      -0.145704      -0.137138      -0.121571 
# Pseudomonas     Prevotella          YRC22   Sphingomonas 
# -0.046613      -0.019928      -0.010305      -0.009905 

### Strongest association for Immune cells on component 1
factor(names(feat_mat0[,'Ruminococcus'][abs(feat_mat0[,'Ruminococcus'])  > 0.5]))
#[1] Coprococcus  Ruminococcus Treponema    U4           Leukocytes  


### Most discriminant at day 0
nm.disc0 = names(sort(abs(sgccda.rs0$loadings$nmr[,1])[abs(sgccda.rs0$loadings$nmr[,1])>0.25],
           decreasing = T))
cl.disc0 = names(sort(abs(sgccda.rs0$loadings$clin[,1])[abs(sgccda.rs0$loadings$clin[,1])>0.25],
           decreasing = T))
gn.disc0 = names(sort(abs(sgccda.rs0$loadings$genus[,1])[abs(sgccda.rs0$loadings$genus[,1])>0.25],
           decreasing = T))



###=====================================================================================
###= COVARIATION data analysis - DIABLO for TR vs. TS under infection (day 132)
###=====================================================================================

###-- Reorder every dataset as NMR profiles (by ind Day)
metatot = metaotu[metaotu$Time %in% c(132),] ## Metadata
metatot = metatot[order(metatot$ID,metatot$Time),]

clin = clinical[clinical$ID %in% meta$Ind & clinical$Day %in% c('D5'),]
stat = clin$TG ## Define status we are interested in
clin = clin[order(clin$ID,clin$Day),-c(1:4)]
clin$FEC = NULL ## remove FEC data as it is used to build status variable

#####--------- Prep ASVs data

## Filter Genus data
# keep day5
genus5 = genus[rownames(genus) %in% metatot$SampleID_Sequencing,] 
result.filter = low.count.removal(genus5, percent = 0.05)
genus.filter = result.filter$data.filter
length(result.filter$keep.otu) #
##[1] 31 ASVs

# function is applied to each row (i.e. each sample)
genuscov = mixOmics::logratio.transfo(result.filter$data.filter,
                                      logratio = 'CLR',offset =1) #t(apply(result.filter$data.filter, 1, clr))
dim(genuscov)
#[1] 16 31

summary(unlist(sapply(apply(genuscov,2, shapiro.test),function(x) x[1])))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.711   0.842   0.929   0.897   0.950   0.980  

### NMR data
nmrcov = bc.dnoised.cln[meta$Day=='RB',]

### rownames
rownames(clin) = seq(1:dim(clin)[1])
rownames(genuscov) = seq(1:dim(genuscov)[1])
rownames(nmrcov) = seq(1:dim(nmrcov)[1])

### X matrix
data = list(genus = genuscov,
            nmr = nmrcov,
            clin = clin)

# check dimension
lapply(data, dim)
# $genus
# [1] 17 31
# 
# $nmr
# [1]  17 110
# 
# $clin
# [1] 17 19

#####------------------- D5 - DIABLO
### Get insight of correlations between datasets
tot.pls1 <- pls(genuscov, clin, ncomp = 1, mode = "regression")
cor(tot.pls1$variates$X, tot.pls1$variates$Y)
#       comp1
# comp1 0.8

tot.pls2 <- pls(genuscov, nmrcov, ncomp = 1, mode = "regression")
cor(tot.pls2$variates$X, tot.pls2$variates$Y)
#       comp1
# comp1 0.7754

tot.pls3 <- pls(nmrcov, clin, ncomp = 1, mode = "regression")
cor(tot.pls3$variates$X, tot.pls3$variates$Y)
#         comp 1
# comp 1   0.6567

### Parameter choice
# Indicate what parameter should be connected
# A compromise needs to be achieved between maximising the correlation between data sets (design value between 0.5 and 1) 
# and maximising the discrimination with the outcome Y (design value between 0 and 0.5)
design1 = matrix(c(0, 0.7754, 0.8,
                   0.7754, 0, 0.6567,
                   0.8, 0.6567, 0), 
                 byrow = T,
                 ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
design1 
#       genus   nmr  clin
# genus 0.0000 0.7754 0.8000
# nmr   0.7754 0.0000 0.6567
# clin  0.8000 0.6567 0.0000

### Use block.plsda to infer relationship between K's
nc = 2 ## K - 1

sgccda.res1 = block.plsda(X = data, Y = stat, design = design1, 
                          near.zero.var = T,
                          ncomp = nc)

perf.1 = perf(sgccda.res1, validation = 'Mfold', folds = 5, 
              nrepeat = 10)
perf.1$error.rate
# $genus
# max.dist centroids.dist mahalanobis.dist
# comp1   0.6529         0.6294           0.6294
# comp2   0.6118         0.6000           0.5882
# 
# $nmr
# max.dist centroids.dist mahalanobis.dist
# comp1   0.4118         0.4294           0.4294
# comp2   0.5529         0.5000           0.5529
# 
# $clin
# max.dist centroids.dist mahalanobis.dist
# comp1   0.5294         0.5118           0.5118
# comp2   0.5941         0.5765           0.5824

#### Optimize number of features with updated design
ncomp = perf.1$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp 
#[1] 1

plotIndiv(sgccda.res1, legend = TRUE, title = 'DIABLO',pch = stat)

## Plot matrix relationship on component 1
plotDiablo(sgccda.res1, ncomp = 1) ## R vs. S
#plotDiablo(sgccda.res1, ncomp = 2)

## Update design using estimated correlations
design2 = matrix(c(0, 0.73, 0.75,
                   0.73, 0, 0.58,
                   0.75, 0.58, 0), 
                 byrow = T,
                 ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
design2
#       genus  nmr clin
# genus  0.00 0.73 0.75
# nmr    0.73 0.00 0.58
# clin   0.75 0.58 0.00

## Optimization for the number of features to be retained in each block
set.seed(123) # for reproducibility, only when the `cpus' argument is not used
test.keepX = list (genus = seq(10,20,5),
                   nmr = seq(10,80,10),
                   clin = c(1))

tune.TCGA = tune.block.splsda(X = data, Y = stat, ncomp = 2, 
                              test.keepX = test.keepX, 
                              near.zero.var = T,design = design2,
                              validation = 'Mfold', folds = 5, nrepeat = 10,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX
list.keepX
# $genus
# [1] 15 15
# 
# $nmr
# [1] 10 10
# 
# $clin
# [1] 9 5

list.keepX = list(genus = c(15,15),nmr = c(10,10),clin = c(9,5))

### Final model with appropriate design and number of features
sgccda.res = block.splsda(X = data, Y = stat, ncomp = 2, 
                          near.zero.var = T,
                          keepX = list.keepX, design = design2)
sgccda.res$design
#       genus  nmr clin Y
# genus  0.00 0.73 0.76 1
# nmr    0.73 0.00 0.57 1
# clin   0.76 0.57 0.00 1
# Y      1.00 1.00 1.00 0

plotIndiv(sgccda.res, legend = TRUE, title = 'DIABLO',pch = stat)

### Selected variables // status - component 1
gn5 = factor(sort(selectVar(sgccda.res, block = 'genus', comp = 1)$genus$name))
gn5
# [1] "Anaeroplasma"  "Anaerovibrio"  "BF311"         "Blautia"       "Clostridium.1"
# [6] "Coprococcus"   "Fibrobacter"   "Mogibacterium" "Paludibacter"  "Prevotella"   
# [11] "Prevotella.1"  "Ruminococcus"  "Sphingomonas"  "Sutterella"    "Treponema" 

nm5 = factor(sort(selectVar(sgccda.res, block = 'nmr', comp = 1)$nmr$name))
nm5
# [1] Alanine                  Choline.Creatinine       Lysine                  
# [4] Methanol                 Phenylalanine.2          Phenylalanine.3         
# [7] Phenylalanine.Tryptophan Serine                   U14                     
# [10] U25  

cl5 = factor(sort(selectVar(sgccda.res, block = 'clin', comp = 1)$clin$name))
cl5
#[1] ADG         ALP         Globins     Haematocrit Leukocytes  MGV         Monocytes   Neutrophil  TP  

plotDiablo(sgccda.res, ncomp = 1, col = colsRS) ## R vs. S
plotIndiv(sgccda.res, legend = TRUE, title = 'DIABLO',pch = stat)
plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = F, 
        pch = c(16, 17, 15), cex = c(.8,.8,.8))

### Loadings on component 1, associated with R / S differences
pdf(file = './Supplementary_Figure4.pdf', width = 14, height = 8)
plotLoadings(sgccda.res, comp = 1, contrib = 'max', #method = 'median',
             size.name = 1.1,legend = F,legend.color = colsRS,
             subtitle = c('Component 1 \n Bacteria genera',
                          'Component 1 \n Metabolites',
                          'Component 1 \n Clinical features'))
dev.off()

#### Circos plot 
pdf(file = './Figure4.pdf',width = 14, height = 14)
circp = circosPlot(sgccda.res, comp = 1, 
                   cutoff = 0.45, line = TRUE, size.variables = 1.6,
                   color.blocks = viridis_pal(option = 'D')(10)[c(5,8,10)],
                   color.Y = colsRS, legend = F,
                   color.cor = c("chocolate3","grey20"), size.labels = .0000001)
dev.off()

### Heatmap for features associated with component 1
gn5 = as.character(gn5)
nm5 = as.character(nm5)
cl5 = as.character(cl5)

feat_mat = circp[c(gn5,nm5,cl5),c(gn5,nm5,cl5)]
colnames(feat_mat) = gsub('X.','',colnames(feat_mat))
rownames(feat_mat) = gsub('X.','',rownames(feat_mat))

pdf(file='./Supplementary_Figure3.pdf')
print(pheatmap::pheatmap(feat_mat,
                         cluster_rows = T,cluster_cols = T,
                         na_col = 'white',
                         fontsize = 9))
dev.off()

#### Most important contributor to TR vs. TS
sort(abs(sgccda.res$loadings$clin[,1])[abs(sgccda.res$loadings$clin[,1])>0.25],
     decreasing = T)
# Haematocrit  Neutrophil   Monocytes         ADG 
#      0.7355      0.4285      0.3186      0.2991 

clin.disc = names(sort(abs(sgccda.res$loadings$clin[,1])[abs(sgccda.res$loadings$clin[,1])>0.25],
                       decreasing = T))
clin.disc
#[1] "Haematocrit" "Neutrophil"  "Monocytes"   "ADG"

## Differences between TR and TS ponies for these
clin132 = clinical[clinical$Day=='D5',]

cbind(aggregate(Haematocrit ~ TG, data = clin132, FUN = mean),
      aggregate(Haematocrit ~ TG, data = clin132, FUN = sd)[,2])
#   TG Haematocrit aggregate(Haematocrit ~ TG, data = clin132, FUN = sd)[, 2]
# 1 TR       39.74                                                      4.145
# 2 TS       41.45                                                      1.745

cbind(aggregate(Monocytes ~ TG, data = clin132, FUN = mean),
      aggregate(Monocytes ~ TG, data = clin132, FUN = sd)[,2])
# TG Monocytes aggregate(Monocytes ~ TG, data = clin132, FUN = sd)[, 2]
# 1 TR     1.178                                                   0.2386
# 2 TS     1.038                                                   0.2446

cbind(aggregate(Neutrophil ~ TG, data = clin132, FUN = mean),
      aggregate(Neutrophil ~ TG, data = clin132, FUN = sd)[,2])
#   TG Neutrophil aggregate(Neutrophil ~ TG, data = clin132, FUN = sd)[, 2]
# 1 TR       2.80                                                    0.5123
# 2 TS       3.05                                                    0.8177

cbind(aggregate(ALP ~ TG, data = clin132, FUN = mean),
      aggregate(ALP ~ TG, data = clin132, FUN = sd)[,2])
#   TG   ALP aggregate(ALP ~ TG, data = clin132, FUN = sd)[, 2]
# 1 TR 5.394                                             0.1510
# 2 TS 5.252                                             0.1247

sort(abs(sgccda.res$loadings$genus[,1])[abs(sgccda.res$loadings$genus[,1])>0.25],decreasing = T)
# Fibrobacter    Treponema Anaerovibrio   Prevotella Anaeroplasma 
# 0.6218       0.4459       0.4018       0.2929       0.2712 

gen.disc = names(sort(abs(sgccda.res$loadings$genus[,1])[abs(sgccda.res$loadings$genus[,1])>0.25],decreasing = T))
gen.disc
#[1] "Fibrobacter"  "Treponema"    "Anaerovibrio" "Prevotella"   "Anaeroplasma"
sort(abs(sgccda.res$loadings$nmr[,1])[abs(sgccda.res$loadings$nmr[,1])>0.25],decreasing = T)
# Phenylalanine.3                 Methanol Phenylalanine.Tryptophan 
#          0.5668                   0.5461                   0.4025 
#    U14          Phenylalanine.2 
# 0.3528                   0.2928 

nm.disc = names(sort(abs(sgccda.res$loadings$nmr[,1])[abs(sgccda.res$loadings$nmr[,1])>0.25],decreasing = T))
nm.disc
# [1] "Phenylalanine.3"          "Methanol"                 "Phenylalanine.Tryptophan"
# [4] "U14"                      "Phenylalanine.2"

### Strongest association for ADG
factor(names(feat_mat[,'ADG'][(feat_mat[,'ADG']) >= 0.45]))
#[1] Anaerovibrio Fibrobacter  Prevotella   Treponema

factor(names(feat_mat[,'ADG'][(feat_mat[,'ADG']) <= -0.45]))
# [1] Lysine                   Phenylalanine.2          Phenylalanine.3         
# [4] Phenylalanine.Tryptophan Serine                   Haematocrit   

feat_mat[c('ADG','Anaerovibrio','Fibrobacter','Prevotella','Treponema'),
          c('ADG','Anaerovibrio','Fibrobacter','Prevotella','Treponema')]
#                 ADG Anaerovibrio Fibrobacter Prevotella Treponema
# ADG          0.4043       0.4973      0.5351     0.5102    0.5106
# Anaerovibrio 0.4973       0.6117      0.6582     0.6276    0.6281
# Fibrobacter  0.5351       0.6582      0.7082     0.6753    0.6759
# Prevotella   0.5102       0.6276      0.6753     0.6439    0.6444
# Treponema    0.5106       0.6281      0.6759     0.6444    0.6450

### Strongest association for Monocytes
factor(names(feat_mat[,'Monocytes'][abs(feat_mat[,'Monocytes']) >= 0.2]))
# [1] Anaerovibrio             Fibrobacter              Prevotella              
# [4] Treponema                Phenylalanine.2          Phenylalanine.3         
# [7] Phenylalanine.Tryptophan Haematocrit  

factor(names(feat_mat[,'Neutrophil'][(feat_mat[,'Neutrophil']) >= 0.45]))
# [1] Phenylalanine.2          Phenylalanine.3          Phenylalanine.Tryptophan
# [4] Haematocrit   

factor(names(feat_mat[,'Neutrophil'][(feat_mat[,'Neutrophil']) <= -0.45]))
#[1] Anaerovibrio Fibrobacter  Prevotella   Treponema  

### Strongest association for Phenylalanine
factor(names(feat_mat[,'Phenylalanine.2'][abs(feat_mat[,'Phenylalanine.2']) >= 0.2]))
# [1] Anaerovibrio             BF311                    Fibrobacter             
# [4] Mogibacterium            Paludibacter             Prevotella              
# [7] Sphingomonas             Treponema                Alanine                 
# [10] Choline.Creatinine       Lysine                   Methanol                
# [13] Phenylalanine.2          Phenylalanine.3          Phenylalanine.Tryptophan
# [16] Serine                   U14                      U25                     
# [19] ADG                      ALP                      Haematocrit             
# [22] Leukocytes               Neutrophil

#### Model performance
set.seed(123)# for reproducibility, only when the `cpus' argument is not used
perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 3, nrepeat =5, 
                   dist = 'centroids.dist')
#perf.diablo  # lists the different outputs

# Performance with Majority vote
perf.diablo$MajorityVote.error.rate
# $centroids.dist
#              comp1  comp2
# TR          0.5111 0.6889
# TS          0.6000 0.7500
# Overall.ER  0.5529 0.7176
# Overall.BER 0.5556 0.7194

# Performance with Weighted prediction
perf.diablo$WeightedVote.error.rate
# $centroids.dist
#              comp1  comp2
# TR          0.5111 0.6889
# TS          0.6000 0.7500
# Overall.ER  0.5529 0.7176
# Overall.BER 0.5556 0.7194

###=====================================================================================
###= Validation of features of interest at day 132 under infection
###=====================================================================================

### Use of blood data recorded in 2020 on more heavily infected individuals
#####------------------- Validation for neutrophil count 
#####------------------- on an independent set of high-sedders 

valid_blood = read.csv(file ='./data/BloodValid.csv',header=T,sep=';')
dim(valid_blood)
#[1] 25 13
valid_blood = valid_blood[valid_blood$FEC >0,] ## discard one non-infected pony
dim(valid_blood)
#[1] 24 13

valid_blood = valid_blood[,c('Monocyte','Neutrophil', 'Haematocrit','FEC')]
valid_blood$TG = 'Validation set'
valid_blood = valid_blood[,c('TG','Monocyte','Neutrophil', 'Haematocrit','FEC')]
valid_blood$Year = 2020

rcorr(as.matrix(valid_blood[,c('Monocyte','Neutrophil', 'Haematocrit','FEC')]),
      type = 'spearman')
#             Monocyte Neutrophil Haematocrit   FEC
# Monocyte        1.00       0.27       -0.25 -0.30
# Neutrophil      0.27       1.00       -0.10 -0.63
# Haematocrit    -0.25      -0.10        1.00 -0.07
# FEC            -0.30      -0.63       -0.07  1.00
# 
# n= 24 
# 
# 
# P
#                Monocyte Neutrophil Haematocrit FEC   
# Monocyte             0.2062     0.2458      0.1493
# Neutrophil  0.2062              0.6426      0.0009
# Haematocrit 0.2458   0.6426                 0.7453
# FEC         0.1493   0.0009     0.7453

### Testing set
blood2015 = clin132[,c('TG','Monocytes','Neutrophil', 'Haematocrit','FEC')]
colnames(blood2015)[2]='Monocyte'
blood2015$Year = 2015

bloodTot = rbind(blood2015,valid_blood)

p2015.1 = ggplot(blood2015,aes(x = TG,y = Haematocrit, fill = TG)) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  geom_point(aes(fill = TG),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colsRS) +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab('Haematocrit (%)') + xlab('') +
  scale_y_continuous(limits = c(32,45), breaks = seq(32,45,2)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 10))

p2015.2 = ggplot(blood2015,aes(x = TG,y = Monocyte, fill = TG)) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  geom_point(aes(fill = TG),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colsRS) +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab('Monocyte count \n (M cells/mm3)') + xlab('') +
  scale_y_continuous(limits = c(0.6,1.6), breaks = seq(0.6,1.6,0.2)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 10))

p2015.3 = ggplot(blood2015,aes(x = TG,y = Neutrophil, fill = TG)) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  geom_point(aes(fill = TG),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colsRS) +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab('Neutrophil count \n (M cells/mm3)') + xlab('') +
  scale_y_continuous(limits = c(1.5,5), breaks = seq(1.5,5,0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 10))

p2020.1 = ggplot(valid_blood,aes(x = FEC,y = Haematocrit)) +
  geom_point(colour="black",size=2, alpha = .8) +
  geom_smooth(method = 'lm',alpha = .2, lty = 3, lwd = .8, col = 'black') +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab('Haematocrit (%)') + xlab('Faecal Egg Count (eggs/g)') +
  scale_y_continuous(limits = c(22,40), breaks = seq(22,40,4)) +
  scale_x_log10() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 10))

p2020.2 = ggplot(valid_blood,aes(x = FEC,y = Monocyte)) +
  geom_point(colour="black",size=2, alpha = .8) +
  geom_smooth(method = 'lm',alpha = .2, lty = 3, lwd = .8, col = 'black') +
  theme_classic() +
  scale_x_log10() +
  theme(legend.position = 'none') +
  ylab('Monocyte count \n (M cells/mm3)') + xlab('Faecal Egg Count (eggs/g)') +
  scale_y_continuous(limits = c(0.4,1.2), breaks = seq(0.4,1.2,0.2)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 10))

p2020.3 = ggplot(valid_blood,aes(x = FEC,y = Neutrophil)) +
  geom_point(colour="black",size=2, alpha = .8) +
  geom_smooth(method = 'lm',alpha = .2, lty = 3, lwd = .8, col = 'black') +
  theme_classic() +
  scale_x_log10() +
  theme(legend.position = 'none') +
  ylab('Neutrophil count \n (M cells/mm3)') + xlab('Faecal Egg Count (eggs/g)') +
  scale_y_continuous(limits = c(1,4.5), breaks = seq(1,4.5,0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 10))


#####---------------=========== Run validation on Peachey et al. 2019 Sci rep dataset
#####---------------=========== 16S data
meta_valid_otu = read.csv(file='./OTU_Peachey_minFreq1_minSamp2/sample-metadata_v6.csv', 
                          sep=',',header=T)
## engineering to retrive true subject id
meta_valid_otu$ID =  sapply(stringr::str_split(meta_valid_otu$subject,' '), function(x) x[1])
meta_valid_otu$ID = substr(meta_valid_otu$ID, 1, nchar(meta_valid_otu$ID)-1)
LOW = meta_valid_otu$ID[meta_valid_otu$Infection==0 & meta_valid_otu$Timepoint=='a']
length(LOW)
#[1] 14
HIGH = meta_valid_otu$ID[meta_valid_otu$Infection==1 & meta_valid_otu$Timepoint=='a']
length(HIGH)
#[1] 9

metap = meta_valid_otu[meta_valid_otu$ID %in% c(HIGH,LOW) & 
                         meta_valid_otu$Infection!=2 &
                         meta_valid_otu$Timepoint=='a',]
dim(metap)
#[1] 23  8

valid_otu = read.csv('./Genus_count_Peachey.tsv', sep = '\t', header= T)
colnames(valid_otu)=gsub('X','',colnames(valid_otu))
valid_otu$Genus= as.character(valid_otu$Genus)
valid_otu$Genus[22]='Clostridium.1'
valid_otu$Genus = factor(valid_otu$Genus)
rownames(valid_otu) = valid_otu$Genus
valid_otu = valid_otu[,-1]
valid_otu = data.frame(t(valid_otu))
valid_otu = valid_otu[match(metap$Sample.id,rownames(valid_otu)),]
colnames(valid_otu)=gsub('X.','',colnames(valid_otu))
valid_otu = mixOmics::logratio.transfo(as.matrix(valid_otu),logratio = 'CLR',offset =1) 
dim(valid_otu)
#[1] 23 61

#### most discriminant bacterial genera under infection
valid_otu = valid_otu[,colnames(valid_otu) %in% gen.disc]
dim(valid_otu)
#[1] 23 14

vo = cbind(metap[,c('subject','Infection')],valid_otu)
v_otu_melt = reshape2::melt(vo,c('subject','Infection'))

####------- Statistical test
mod = NULL
for(i in unique(v_otu_melt$variable)){
  tt = wilcox.test(value ~ Infection, data= v_otu_melt[v_otu_melt$variable ==i,])$p.value
  g = i
  tmp = data.frame(genus = i, test = tt)
  mod = rbind(mod,tmp)
  rm(tmp,g,tt)
}
rm(p)
p = mod$test
mod$fdr = p.adjust(mod$test, method = 'BH', n = length(p))

mod[mod$fdr<0.05,]
# [1] genus test  fdr  
# <0 rows> (or 0-length row.names)

mod[mod$test<0.05,]
#           genus    test    fdr
# 2    Prevotella 0.04559 0.2279

peachey_sig = v_otu_melt[v_otu_melt$variable %in% mod$genus[mod$test<.05],]
peachey_sig$stat = 'High'
peachey_sig$stat[peachey_sig$Infection==0] = 'Low'
peachey_sig$stat = factor(peachey_sig$stat,levels = c('Low', 'High'))
peachey_sig$variable = as.character(peachey_sig$variable)
peachey_sig = peachey_sig[order(peachey_sig$variable),]
peachey_sig$variable = factor(peachey_sig$variable)

aggregate(value ~ Infection + variable, data = peachey_sig, FUN=mean)
#   Infection   variable value
# 1         0 Prevotella 5.245
# 2         1 Prevotella 4.757

## Our data
tmp_otu = cbind(data.frame(stat),data.frame(meta$TG[meta$Day=='RB']),
                genuscov[,gen.disc])
colnames(tmp_otu)[2] ='TG'
tmp_otu = reshape2::melt(tmp_otu,c('stat','TG'))
tmp_otu$variable = as.character(tmp_otu$variable)
tmp_otu$variable = factor(tmp_otu$variable, levels = unique(sort(tmp_otu$variable)))

dt = tmp_otu[tmp_otu$variable %in% mod$genus[mod$test<.05],]
aggregate(value ~ stat + variable, data = dt, FUN=mean)
#   stat   variable value
# 1   TR Prevotella 3.269
# 2   TS Prevotella 2.983

## Bacteria count
p1 = ggplot(dt,aes(x = stat, y = value, fill = stat)) +
  geom_boxplot(alpha = .4) + 
  scale_fill_manual(values = colsRS) + ylab('Prevotella abundance \n (CLR-transformed)') +
  ggtitle("") + xlab('') + 
  theme_classic() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size=10)) 

p2 = ggplot(peachey_sig,aes(x = stat, y = value, fill = stat)) +
  geom_boxplot(alpha = .4) + 
  scale_fill_manual(values = colsStat) +ylab('Prevotella abundance \n (CLR-transformed)') +
  ggtitle("") + xlab('') +
  theme_classic() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size=10)) 

pdf(file = './Figure5.pdf')
multiplot(p2015.1,p2015.2,p2015.3,p1,
          p2020.1,p2020.2,p2020.3,p2,
          cols = 2)
dev.off()

# multiplot(p2015.2,p2015.3,
#           p2020.2,p2020.3,
#           cols = 2)

###=====================================================================================
###= Core features discriminating TR and TS across conditions
###=====================================================================================

#####------------------- Clinical parameters found differential between TR and TS ponies
clc = as.character(na.omit(clin.disc[match(clin.disc,cl.disc0)]))
clc
#[1] "Neutrophil"

clin = clinical[clinical$ID %in% meta$Ind & clinical$Day %in% c('D1','D5'),]
clin = clin[order(clin$ID,clin$Day),c('ID','Day','TG',clc)]
clin$day[clin$Day=='D1']='day0'
clin$day[clin$Day=='D5']='day132'
clin$Day=NULL
clin = reshape2::melt(clin,c(1,2,dim(clin)[2]))
#clin$variable = 'Neutrophil'
colnames(clin)=c('ID','stat','day','variable','value')
clin = clin[,c('ID','stat','variable','value','day')]

#####------------------- Metabolite found differential between TR and TS ponies
nmc = as.character(na.omit(nm.disc[match(nm.disc,nm.disc0)]))
nmc
#character(0)

nmc = as.character(na.omit(nm0[match(nm5,nm0)]))
nmc
#[1] "Lysine"

##NMR at day 0
nmr0 = bc.dnoised.cln[meta$Day=='J0',nmc] ## day0
nmr0 = data.frame(nmr0)
nmr0$stat = stat
nmr0$day = 'day0'
nmr0$ID = meta$Ind[meta$Day=='J0']
colnames(nmr0)[1]='value'
##NMR at day 132
nmr5 = bc.dnoised.cln[meta$Day=='RB',nmc] ## day0
nmr5 = data.frame(nmr5)
nmr5$stat = stat
nmr5$day = 'day132'
nmr5$ID = meta$Ind[meta$Day=='RB']
colnames(nmr5)[1]='value'
nmrtot = rbind(nmr0,nmr5)
nmrtot$variable = paste0(nmr_peak$Metabolite_Reannotated[match(nmc,nmr_peak$Metabolite_Reannotated)],
                         ' \n ',
                         nmr_peak$Chemical_shift[match(nmc,nmr_peak$Metabolite_Reannotated)],
                         ' ppm')
nmrtot = nmrtot[,c('ID','stat','variable','value','day')]

#####------------------- Genera found differential between TR and TS ponies
#### Common parameters between day0 and day5
gnc = as.character(na.omit(gn0[match(gen.disc,gn0)]))
gnc
#[1] "Coprococcus"  "Fibrobacter"  "Ruminococcus" "Treponema" 

#### Is any of these parameters associated with FEC at housing ?
### Genera data at day0
samp = metaotu$SampleID_Sequencing[metaotu$Time %in% c(0)] ## Metadata
bact = genus[rownames(genus) %in% samp,]
bact = data.frame(bact[,gnc])
bact = mixOmics::logratio.transfo(bact,'CLR',offset=1)
class(bact) <- "matrix"
bact = data.frame(bact)
bact$stat = stat

bactplot.rs = melt(bact,c(dim(bact)[2]))
bactplot.rs$day = 'day0'
bactplot.rs$ID = metaotu$ID[metaotu$Time %in% c(0)]
rm(bact)

# Retrieve genera data at day 132
samp = metaotu$SampleID_Sequencing[metaotu$Time %in% c(132)] ## Metadata
bact = genus[rownames(genus) %in% samp,]
bact = data.frame(bact[,gnc])
bact = mixOmics::logratio.transfo(bact,'CLR',offset=1)
class(bact) <- "matrix"
bact = data.frame(bact)
bact$stat = stat

bactplot.TG = melt(bact,c(dim(bact)[2]))
bactplot.TG$day = 'day132'
bactplot.TG$ID = metaotu$ID[metaotu$Time %in% c(132)]

bactot = rbind(bactplot.rs,bactplot.TG)
bactot = bactot[,c('ID','stat','variable','value','day')]
colstot = c(colsRS,colsStat)

## Altogether
tot = rbind(clin,nmrtot,bactot)
bac = factor(unique(bactot$variable))

pS8.1 = ggplot(tot[tot$variable=='Leukocytes',],aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colstot) +
  theme_classic() +
  theme(legend.position = 'none') + ggtitle('Leukocytes') +
  ylab('Leukocyte count \n (M cells/mm3)') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 11),
        axis.title = element_text(size = 10))
pS8.2 = ggplot(tot[tot$variable=='Neutrophil',],aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colstot) +
  theme_classic() +
  theme(legend.position = 'none') + ggtitle('Neutrophil') +
  ylab('Neutrophil count \n (M cells/mm3)') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 11),
        axis.title = element_text(size = 10))
pS8.3 = ggplot(tot[tot$variable=='Lysine \n 1.445-1.465 ppm',],aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colstot) +
  theme_classic() + ggtitle('Lysine') +
  theme(legend.position = 'none') +
  ylab('Metabolite peak intensity') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 11),
        axis.title = element_text(size = 10))
pS8.4 = ggplot(tot[tot$variable=='Coprococcus',],aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colstot) +
  theme_classic() + ggtitle('Coprococcus') +
  theme(legend.position = 'none') +
  ylab('Bacterial abundance \n (CLR-transformed)') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 11),
        axis.title = element_text(size = 10))
pS8.5 = ggplot(tot[tot$variable=='Fibrobacter',],aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colstot) +
  theme_classic() +
  theme(legend.position = 'none') + ggtitle('Fibrobacter') +
  ylab('Bacterial abundance \n (CLR-transformed)') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 11),
        axis.title = element_text(size = 10))
pS8.6 = ggplot(tot[tot$variable=='Ruminococcus',],aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colstot) +
  theme_classic() +
  theme(legend.position = 'none') +ggtitle('Ruminococcus') +
  ylab('Bacterial abundance \n (CLR-transformed)') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 11),
        axis.title = element_text(size = 10))
pS8.7 = ggplot(tot[tot$variable=='Treponema',],aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = colstot) +
  theme_classic() +
  theme(legend.position = 'none') + ggtitle('Treponema') +
  ylab('Bacterial abundance \n (CLR-transformed)') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 11),
        axis.title = element_text(size = 10))

pdf(file='./supplementary_Figure8.pdf',width = 14, height = 8)
multiplot(pS8.4,pS8.1,pS8.5,pS8.2,
          pS8.6,pS8.3,pS8.7,cols=4)
dev.off()

### Identify significant variations
tmp = NULL
for(i in unique(tot$variable)){
  s = summary(lme(value ~ stat*day,
                  random =~ 1|ID,
                  data=tot[tot$variable==i,]))
  tmp = rbind(tmp,coefficients(s)[2:4,])
}

tmp = data.frame(tmp)
tmp$Bact = rep(unique(tot$variable),each =3)

p = tmp$p.value
tmp$fdr = p.adjust(p, method = 'BH', n = length(p))
tmp$Contrast = rep(c('TS vs. TR', 'Day 132 vs day 0', 
                     'Differential temporal trend in TS'),
                   n=length(unique(tot$variable)))
colnames(tmp) = c('Estimate', 'Standard Error','DF', 
                  't-value', 'P-value', 'Variable', 'FDR', 'Contrast')
rownames(tmp) = seq(1:dim(tmp)[1])
tmp = tmp[,c('Variable','Contrast','Estimate', 
             'Standard Error','DF', 't-value', 'P-value','FDR')]
tmp[tmp$`P-value`<.05,]
#                     Variable                          Contrast   Estimate Standard Error DF
# 2                 Leukocytes                  Day 132 vs day 0  1.8411111      4.067e-01 15
# 5                 Neutrophil                  Day 132 vs day 0  0.6666667      2.742e-01 15
# 8  Lysine \n 1.445-1.465 ppm                  Day 132 vs day 0 -0.0002453      6.546e-05 15
# 9  Lysine \n 1.445-1.465 ppm Differential temporal trend in TS  0.0002514      9.542e-05 15
# 15               Fibrobacter Differential temporal trend in TS -0.5527175      1.777e-01 15
# 20                 Treponema                  Day 132 vs day 0 -0.2238013      8.345e-02 15
#    t-value   P-value      FDR
# 2    4.527 0.0004012 0.008426
# 5    2.431 0.0280428 0.098150
# 8   -3.748 0.0019403 0.020373
# 9    2.634 0.0187674 0.078823
# 15  -3.110 0.0071671 0.050170
# 20  -2.682 0.0170700 0.078823

write.csv(tmp, file = './Supplementary_Table3.csv', quote = F, row.names = F)

# require(dplyr)
# dt = tot132 %>% dcast(ID + FEC ~ variable, value.var = 'value')
# dt$ID = NULL
# rcorr(as.matrix(dt,type='spearson'))$r['FEC',]
# # FEC                Leukocytes                       MGV                 Monocytes                Neutrophil 
# # 1.00000                   0.23222                   0.05693                  -0.11903                   0.09030 
# # Lysine \n 1.445-1.465 ppm                     BF311               Coprococcus               Fibrobacter              Paludibacter 
# # 0.32418                   0.30717                   0.03886                  -0.65401                   0.18032 
# # Prevotella              Ruminococcus              Sphingomonas                 Treponema 
# # -0.33715                   0.22252                  -0.11855                  -0.37317
# 
# rcorr(as.matrix(dt,type='spearson'))$P['FEC',]
# # FEC                Leukocytes                       MGV                 Monocytes                Neutrophil 
# # NA                    0.3698                    0.8282                    0.6491                    0.7304 
# # Lysine \n 1.445-1.465 ppm                     BF311               Coprococcus               Fibrobacter              Paludibacter 
# #                    0.2043                    0.2304                    0.8823                    0.0044                    0.4886 
# # Prevotella              Ruminococcus              Sphingomonas                 Treponema 
# #     0.1857                    0.3907                    0.6504                    0.1401

#####---------------=========== Lysine in Peachey et al. 2019 Sci rep dataset

## Fetch corresponding FEC data for Fibrobacter
tot132 = tot[tot$day=='day132' & tot$variable=='Fibrobacter',]
subFEC = clinical[clinical$Day=='D5',c('ID','FEC')]
tot132$FEC = subFEC$FEC[match(tot132$ID,subFEC$ID)]

pFibroFEC = ggplot(tot132,aes(x = exp(FEC)+50,y = value)) +
  geom_point(aes(fill = stat),colour = "black",pch=21, size=3) +
  geom_smooth(method = 'lm',col='black',lwd = .5, lty = 3, 
              fill = 'grey',alpha = .1) +
  scale_fill_manual(values = colstot[c(3,4)]) +
  theme_classic() + 
  ggtitle('D') +
  theme(legend.position = 'none') +
  ylab('Fibrobacter abundance \n (CLR-transformed)') + xlab('Faecal Egg Count (eggs/g)') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12))

pS8.5 = ggplot(tot[tot$variable=='Fibrobacter',],
               aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = rep(colstot[1:2],2)) +
  theme_classic() +
  theme(legend.position = 'none') + ggtitle('C') +
  ylab('Fibrobacter abundance \n (CLR-transformed)') + 
  xlab('Time point') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12))

nmvalid = colnames(nmrvalid)[as.integer(na.omit(match(unique(unlist(stringr::str_split(nmc,'\\.'))),colnames(nmrvalid))))]
nmvalid
#[1] "Lysine"

wilcox.test(nmrvalid[,'Lysine'] ~ factor(metavalid$Infection)) ##
# Wilcoxon rank sum exact test
# 
# data:  nmrvalid[, "Lysine"] by factor(metavalid$Infection)
# W = 16, p-value = 0.002
# alternative hypothesis: true location shift is not equal to 0

####---- NMR for Lysine 
pS8.3 = ggplot(tot[tot$variable=='Lysine \n 1.445-1.465 ppm',],
               aes(x = day,y = value, fill = paste0(day,stat))) +
  geom_boxplot(alpha = .4,outlier.shape=NA) +
  #facet_wrap(~ variable, scales = 'free') +
  geom_point(aes(fill = paste0(day,stat)),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  scale_fill_manual(values = rep(colstot[1:2],2)) +
  theme_classic() + ggtitle('A') +
  theme(legend.position = 'none') +
  ylab('Lysine peak intensity \n 1.445-1.465 ppm') + xlab('') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 11))

####---- Lysine in Peachey's et al.
vn = cbind(metavalid[,c('Subject','Infection')],
           nmrvalid[,nmvalid])
colnames(vn)[3] = nmvalid
vn$Infection = factor(vn$Infection)
v_nmr_melt = reshape2::melt(vn,c('Subject','Infection'))
v_nmr_melt$stat='Low'
v_nmr_melt$stat[v_nmr_melt$Infection==1]='High'
v_nmr_melt$stat=factor(v_nmr_melt$stat,levels=c('Low','High'))

pLysinePeachey = ggplot(v_nmr_melt,aes(x = stat, y = value,fill=stat)) +
  geom_boxplot(alpha = .6, outlier.shape=NA) +
  geom_point(aes(fill = stat),
             position = position_jitterdodge(),colour="black",pch=21, size=2) +
  xlab('') + theme_classic() +
  scale_fill_manual(values = colsStat) + 
  theme(legend.position = 'none') + ylab('Lysine peak intensity') +
  ggtitle("B") + 
  theme(legend.position = 'none',
        text = element_text(size = 12))

pdf(file = './Figure6.pdf')
multiplot(pS8.3, pS8.5, 
          pLysinePeachey,pFibroFEC, 
          cols = 2)
dev.off()


#####---------------============== session info =============--------------
sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] mixOmics_6.12.2            MASS_7.3-53                RVAideMemoire_0.9-78      
# [4] pls_2.7-3                  vegan_2.5-6                permute_0.9-5             
# [7] Hotelling_1.0-5            corpcor_1.6.9              viridis_0.5.1             
# [10] viridisLite_0.3.0          FactoMineR_2.3             factoextra_1.0.7          
# [13] reshape2_1.4.4             tidyr_1.1.2                plyr_1.8.6                
# [16] glmnet_4.0-2               PerformanceAnalytics_2.0.4 xts_0.12.1                
# [19] zoo_1.8-8                  corrplot_0.84              lsmeans_2.30-0            
# [22] emmeans_1.5.2-1            nlme_3.1-150               car_3.0-10                
# [25] carData_3.0-4              ROCR_1.0-11                lmerTest_3.1-3            
# [28] lme4_1.1-25                Matrix_1.2-18              Hmisc_4.4-1               
# [31] Formula_1.2-4              survival_3.2-7             lattice_0.20-41           
# [34] ggplot2_3.3.3             
# 
# loaded via a namespace (and not attached):
# [1] estimability_1.3      ModelMetrics_1.2.2.2  coda_0.19-4           bit64_4.0.5          
# [5] knitr_1.30            multcomp_1.4-15       data.table_1.13.2     rpart_4.1-15         
# [9] doParallel_1.0.16     generics_0.1.0        metap_1.4             BiocGenerics_0.34.0  
# [13] preprocessCore_1.50.0 TH.data_1.0-10        RSQLite_2.2.1         RApiSerialize_0.1.0  
# [17] bit_4.0.4             mutoss_0.1-12         xml2_1.3.2            lubridate_1.7.9      
# [21] gower_0.2.2           xfun_0.19             hms_0.5.3             scrime_1.3.5         
# [25] progress_1.2.2        caTools_1.18.0        readxl_1.3.1          igraph_1.2.6         
# [29] DBI_1.1.0             tmvnsim_1.0-2         htmlwidgets_1.5.2     reshape_0.8.8        
# [33] rARPACK_0.11-0        stats4_4.0.2          purrr_0.3.4           ellipsis_0.3.1       
# [37] RSpectra_0.16-0       ctc_1.62.0            dplyr_1.0.2           ggpubr_0.4.0         
# [41] backports_1.2.0       gbRd_0.4-11           RcppParallel_5.0.2    vctrs_0.3.6          
# [45] Biobase_2.48.0        Cairo_1.5-12.2        abind_1.4-5           caret_6.0-86         
# [49] withr_2.4.0           checkmate_2.0.0       prettyunits_1.1.1     ropls_1.20.0         
# [53] mnormt_2.0.2          cluster_2.1.0         lazyeval_0.2.2        crayon_1.3.4         
# [57] ellipse_0.4.2         edgeR_3.30.3          recipes_0.1.15        pkgconfig_2.0.3      
# [61] labeling_0.4.2        ProtGenerics_1.20.0   nnet_7.3-14           rlang_0.4.10         
# [65] lifecycle_0.2.0       sandwich_3.0-0        affyio_1.58.0         mathjaxr_1.0-1       
# [69] cellranger_1.1.0      matrixStats_0.57.0    graph_1.66.0          boot_1.3-25          
# [73] base64enc_0.1-3       pheatmap_1.0.12       stringfish_0.14.2     png_0.1-7            
# [77] mzR_2.22.0            rsm_2.10.2            bitops_1.0-6          KernSmooth_2.23-18   
# [81] pROC_1.17.0.1         blob_1.2.1            shape_1.4.5           stringr_1.4.0        
# [85] jpeg_0.1-8.1          rstatix_0.6.0         S4Vectors_0.26.1      ggsignif_0.6.0       
# [89] scales_1.1.1          leaps_3.1             memoise_1.1.0         magrittr_2.0.1       
# [93] gplots_3.1.0          gdata_2.18.0          zlibbioc_1.34.0       compiler_4.0.2       
# [97] RColorBrewer_1.1-2    plotrix_3.8-1         pcaMethods_1.80.0     affy_1.66.0          
# [101] htmlTable_2.1.0       mgcv_1.8-33           tidyselect_1.1.0      vsn_3.56.0           
# [105] stringi_1.5.3         forcats_0.5.0         yaml_2.2.1            locfit_1.5-9.4       
# [109] latticeExtra_0.6-29   MALDIquant_1.19.3     ggrepel_0.8.2         fastmatch_1.1-0      
# [113] tools_4.0.2           parallel_4.0.2        rio_0.5.16            rstudioapi_0.13      
# [117] qs_0.23.3             foreach_1.5.1         foreign_0.8-80        crmn_0.0.21          
# [121] gridExtra_2.3         prodlim_2019.11.13    scatterplot3d_0.3-41  farver_2.0.3         
# [125] mzID_1.26.0           digest_0.6.27         BiocManager_1.30.10   lava_1.6.8.1         
# [129] quadprog_1.5-8        ppcor_1.1             Rcpp_1.0.6            siggenes_1.62.0      
# [133] broom_0.7.2           ncdf4_1.17            httr_1.4.2            MSnbase_2.14.2       
# [137] Rdpack_2.1            colorspace_2.0-0      XML_3.99-0.5          IRanges_2.22.2       
# [141] splines_4.0.2         RBGL_1.64.0           statmod_1.4.35        sn_1.6-2             
# [145] multtest_2.44.0       plotly_4.9.2.1        xtable_1.8-4          jsonlite_1.7.2       
# [149] nloptr_1.2.2.2        timeDate_3043.102     glasso_1.11           flashClust_1.01-2    
# [153] ipred_0.9-9           R6_2.5.0              TFisher_0.2.0         pillar_1.4.7         
# [157] htmltools_0.5.0       glue_1.4.2            minqa_1.2.4           BiocParallel_1.22.0  
# [161] class_7.3-17          codetools_0.2-18      fgsea_1.14.0          mvtnorm_1.1-1        
# [165] tibble_3.0.5          numDeriv_2016.8-1.1   huge_1.3.4.1          curl_4.3             
# [169] gtools_3.8.2          zip_2.1.1             openxlsx_4.2.3        limma_3.44.3         
# [173] munsell_0.5.0         iterators_1.0.13      impute_1.62.0         haven_2.3.1          
# [177] gtable_0.3.0          rbibutils_1.4 