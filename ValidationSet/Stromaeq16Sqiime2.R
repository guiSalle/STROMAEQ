library("knitr")
library("BiocStyle")
require('ggplot2')
require('gridExtra')
require('dada2')
require('phyloseq')
require('plyr')
library(data.table)
library(scales)
library(ggridges)
library(dplyr)
library(reshape2)

### Uses code lines from:
# https://f1000research.com/articles/5-1492/v1
# http://joey711.github.io/phyloseq-extensions/DESeq2.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html#import-data-with-phyloseq-convert-to-deseq2
# http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html

set.seed(100)
setwd("~/Documents/INRA/STROMAEQ/Paper/")
theme_set(theme_bw())


######============================================================================
#######----------------------------- OTU analysis from Clarke et al. 2018
######============================================================================

## Working on QIIME2 abundance tables
otu_table = read.table("./OTU_Clarke_minFreq5_minSamp2/dada2_output_mxe1_exported/feature-table.tsv",sep="\t",header=T)
dim(otu_table)
#[1] 6936  101

otu_matrix = read.csv("./OTU_Clarke_minFreq5_minSamp2/dada2_output_mxe1_exported/biom-taxonomy.tsv",sep="\t",header=T)
otu_matrix = as.matrix(otu_matrix)
#row.names(otu_matrix)=otu_matrix[,1]

metadata = read.table("./OTU_Clarke_minFreq5_minSamp2/Mapfie_StromaEq.txt", sep="\t",header=T)
rownames(metadata) = metadata$SampleID_Sequencing

metadata$Time. = factor(metadata$Time.)
metadata$Time = factor(metadata$Time)
metadata$Treatment = factor(metadata$Treatment)
metadata$ID = factor(metadata$ID)

## Data struct
table(metadata$Treatment,metadata$Time)
#    0 24 43 92 132
# R 10 10 10 10  10
# S 10 10 10 10  10

merg = merge(otu_table,otu_matrix,by=c('OTUID'))
dim(merg)
#[1] 6208  103

## Split taxonomy info
merg2 = tidyr::separate(merg,taxonomy,c('Kingdom','Phylum','Class','Order','Family','Genus','Species'),';')
dim(merg2)
#[1] 6208  109

## New tables saved to csv // ALLLL names have to be similar in OTU, TAX ANNDDD phy tree !!!
otu = as.matrix(merg2[,2:(dim(metadata)[1]+1)])
row.names(otu) = merg2$OTUID 
#colnames(otu)=gsub('[.]','-',colnames(otu))

tax = as.matrix(merg2[,(dim(metadata)[1]+2):dim(merg2)[2]])
row.names(tax) = merg2$OTUID 

#####--------------------------------------------- Phyloseq data
OTU = otu_table(as.matrix(otu), taxa_are_rows = T)
TAX = tax_table(as.matrix(tax))
meta = sample_data(metadata)
phy_tree = read_tree("./OTU_Clarke_minFreq5_minSamp2/tree_out/tree_unrooted.1/tree.nwk")
clark = phyloseq(OTU, TAX , meta, phy_tree)
clark
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6208 taxa and 100 samples ]
# sample_data() Sample Data:       [ 100 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 6208 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6208 tips and 6191 internal nodes ]

### Export Genus count table
## Export
gen_count = as(otu_table(clark), 'matrix')
OTUdf = as.data.frame(gen_count)
head(OTUdf)
OTUdf$otuid = rownames(OTUdf)

## Tax
tx = as(tax_table(clark)[,6], 'matrix')
tx = as.data.frame(tx)
tx$otuid = rownames(tx)
#tx$Genus = gsub('g__','',tx$Genus)

#Final
OTUdfe = merge(tx,OTUdf,by='otuid')
OTUdfe$otuid = NULL
head(OTUdfe)
##Rm unassigned Genus
#OTUdfe = OTUdfe[OTUdfe$Genus!=' ',]
dim(OTUdfe)
#[1]  6208 101

write.table(OTUdfe, quote=F,row.names = F,file ='./supplementary_Table2.tsv')


#####===================================================================================================
#### PRE-PROCESSING & DATA EXPLORATION
#####===================================================================================================

#####--------------------------------------------- Seq Depth
sdt = data.table(as(sample_data(clark), "data.frame"),
                 TotalReads = sample_sums(clark), keep.rownames = TRUE)
#setnames(sdt, "rn", "sample.id")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth") + 
  annotation_logticks(scaled = T,sides="b") +
  scale_x_log10(breaks = c(30000,50000,100000,200000), 
                labels = c(30000,50000,100000,200000)) 

pSeqDepth + facet_grid( Treatment ~ Time)

aggregate(TotalReads ~ Treatment + Time, data= sdt, FUN=mean)
#    Treatment Time TotalReads
# 1          R    0    29817.9
# 2          S    0    29719.1
# 3          R   24    28312.2
# 4          S   24    32621.5
# 5          R   43    47445.1
# 6          S   43    24791.6
# 7          R   92    26758.1
# 8          S   92    28097.3
# 9          R  132    20583.6
# 10         S  132    27461.8

## Pups have significantly more reads than others
summary.aov(lm(TotalReads ~ Treatment + Time,data=sdt))
# Df    Sum Sq   Mean Sq F value Pr(>F)
# Treatment    1 1.046e+08 104562895   0.493  0.484
# Time         4 1.582e+09 395430976   1.864  0.123
# Residuals   94 1.994e+10 212091042               

### Need to normalize sequencing depth across samples
summary(sdt$TotalReads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4219   19263   25900   29561   36372   76771 

#####--------------------------------------------- Taxa filtering based on prevalence and total counts
### Total counts
tdt = data.table(tax_table(clark),
                 TotalCounts = taxa_sums(clark),
                 OTU = taxa_names(clark))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts")

# How many unknown Families ? 
length(which(is.na(tdt$Family)))
#[1] 218
# How many unknown Genus ? 
length(which(is.na(tdt$Genus)))
#[1] 529

# How many singletons?
tdt[(TotalCounts <= 0), .N]
#[1] 0
tdt[(TotalCounts <= 1), .N]
#[1] 0

# Min count = 5 
tdt[(TotalCounts < 5), .N]
# [1] 0

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
#- zoom-in
pCumSum + xlim(0, 20000)

###---------------------- Taxa prevalence

## Quicker (from evomics)
source('~/Documents/scripts/taxa_summary.R',local=TRUE)
mdt = fast_melt(clark)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]

ggplot(prevdt, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Histogram of Taxa Prevalence")

## Singletons ?
prevdt[(Prevalence <= 1), .N]
#[1] 0
# How many doubletons?
prevdt[(Prevalence <= 2), .N]
# 1247

prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pPrevCumSum

### Prevalence vs. Total count
ggplot(prevdt, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

## Add Phylum information
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))

# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt2 <- addPhylum[prevdt]
showPhyla = prevdt2[, sum(TotalCounts), by = Phylum][order(-V1)][1:6]$Phylum ## only 6 phyla observed ??
setkey(prevdt2, Phylum)
ggplot(prevdt2[showPhyla], 
       mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_point(size = 4, alpha = 0.75) + 
  geom_vline(xintercept = 4) + geom_hline(yintercept = 5) +
  scale_y_log10()

### Define prevalence threshold as 5% of remaining samples
prevalenceThreshold = 0.05 * nsamples(clark)
prevalenceThreshold
#[1] 5

### Filter entries with unidentified Phylum.
table(addPhylum$Phylum)
# p__   p__Actinobacteria    p__Bacteroidetes    p__Elusimicrobia    p__Fibrobacteres       p__Firmicutes     p__Fusobacteria 
#  26                  85                2065                   9                  63                3234                   2 
# p__Planctomycetes   p__Proteobacteria     p__Spirochaetes      p__Tenericutes  p__Verrucomicrobia 
#                 4                 244                 343                 112                   2 
keepPhyla = table(addPhylum$Phylum)[(table(addPhylum$Phylum) >= 1)]
ff1 = subset_taxa(clark, Phylum %in% names(keepPhyla))
ff1
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6189 taxa and 100 samples ]
# sample_data() Sample Data:       [ 100 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 6189 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6189 tips and 6164 internal nodes ]

## Transform to relative abundance
ff1ra = transform_sample_counts(ff1, function(x){x / sum(x)})

sdt2 = data.table(as(sample_data(ff1), "data.frame"),
                  TotalReads = sample_sums(ff1), keep.rownames = TRUE)


#####--------------------------------------------- Abundances 
### plot rarefaction curve before Genus agglomeration
rc = vegan::rarecurve(t(otu_table(ff1)), step=50, cex=0.5,label=T)

rare <- lapply(rc, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

rare <- purrr::map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")
rare$sample <- levels(sample_data(ff1)$SampleID_Sequencing)[as.integer(rare$sample)]
rare$Treatment = sample_data(ff1)$Treatment[match(rare$sample,sample_data(ff1)$SampleID_Sequencing)]

ggplot(data = rare,aes(x = raw.read, y = OTU, Day=sample,color=Treatment))+
  geom_line() +
  scale_x_continuous(labels =  scales::scientific_format()) +
  ggtitle('Rarefaction curves')

#####===================================================================================================
#### Pre-processing & main observations
#####===================================================================================================

#####--------------------------------------------- Alpha diversity (not to be estimated on filtered data)
## Alpha div estimators account for variation in library sizes
## Shannon influenced by low abundance species
## InvSimpson influenced by major species

pAlpha2 = plot_richness(ff1,
                        shape = "Treatment",
                        color = "Treatment",
                        x = 'Time',
                        measures = c("Shannon"),#c("Observed", "Chao1", "Shannon", "InvSimpson"),
                        title = "Alpha Diversity - Filtering minor OTUs")
pAlpha2 + geom_boxplot() #size = 3)

## By sample
pAlpha3 = plot_richness(ff1,
                        shape = "Time",
                        color = "Treatment",
                        x = 'SampleID_Sequencing',
                        measures = c("Shannon"),
                        title = "Alpha Diversity - Filtering minor OTUs")
pAlpha3 + geom_point(size = 3)

###------------------------------------ Linear modeling of alpha diversity

####-------------- Unagglom. dataset

#### Shannon
ff1_alpha_div <- estimate_richness(ff1, split = TRUE, measure = "Shannon")
ff1_alpha_div$SampleID_Sequencing <- rownames(ff1_alpha_div) %>%  as.factor()
ff1_samp <- sample_data(ff1) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ff1_alpha_div, by = "SampleID_Sequencing") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")

# reorder's facet from lowest to highest diversity
diversity_means <- ff1_samp %>%
  group_by(Treatment,Time) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
# ff1_samp$Stage <- factor(ff1_samp$Stage,
#                          diversity_means$Stage)

alpha_div_model <- lm( alpha_diversity ~ Treatment*Time , data = ff1_samp)
summary.aov(alpha_div_model)
# Treatment       1   0.00   0.000    0.00   0.983    
# Time            4   2.22   0.556    7.11 5.1e-05 ***
# Treatment:Time  4   1.07   0.268    3.43   0.012 *  
# Residuals      90   7.04   0.078  
summary(alpha_div_model)
# Call:
#   lm(formula = alpha_diversity ~ Treatment * Time, data = ff1_samp)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.71374 -0.13406  0.03065  0.16688  0.66211 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         5.85808    0.08843  66.244   <2e-16 ***
# TreatmentS         -0.14895    0.12506  -1.191   0.2368    
# Time24              0.16094    0.12506   1.287   0.2014    
# Time43              0.27980    0.12506   2.237   0.0277 *  
# Time92              0.26995    0.12506   2.159   0.0335 *  
# Time132            -0.20225    0.12506  -1.617   0.1093    
# TreatmentS:Time24   0.29029    0.17686   1.641   0.1042    
# TreatmentS:Time43  -0.16100    0.17686  -0.910   0.3651    
# TreatmentS:Time92   0.18631    0.17686   1.053   0.2950    
# TreatmentS:Time132  0.42312    0.17686   2.392   0.0188 *  

# predictions and confidence intervals.
new_data <- expand.grid(Time = levels(ff1_samp$Time), Treatment = levels(ff1_samp$Treatment))
a = predict(alpha_div_model, newdata = new_data, interval="confidence")
a = data.frame(a)
new_data$pred = a$fit
new_data$lwr= a$lwr
new_data$upr= a$upr

# X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
#                   new_data[-ncol(new_data)])
# pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
# new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2

# fitted values, with error bars
shannon_plot=ggplot(ff1_samp %>% left_join(new_data)) +
  facet_wrap(~ Treatment) +
  geom_errorbar(aes(x = Time, ymin = lwr,
                    ymax = upr),
                col = "#858585", size = .5,width=.2) +
  geom_point(aes(x = Time, y = alpha_diversity,
                 col = Treatment), size = 4) +
  scale_y_continuous(limits = c(5, 7), breaks = seq(5, 7, .5)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Stage", y = "Shannon Diversity", color = "Stage") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90),legend.position = 'bottom') 
shannon_plot

#####--------------------------------------------- Beta diversity

#####------------------ ALL OTUs: PCoA 
### UniFrac accounts for phylogenetic d
### Sequencing depth can greatly affects ordination
### Read counts should be normalized: either with DESeq2 or applying log transf.

## No transf.
DistBC = distance(ff1ra, method = "bray")
DistwUF = distance(ff1ra, method = "wUniFrac")
DistJ = distance(ff1ra, method = "jaccard")
DistuUF = distance(ff1ra, method = "UniFrac")

ordBCp = ordinate(ff1ra, method = "PCoA", distance = DistBC)
ordwUFp = ordinate(ff1ra, method = "PCoA", distance = DistwUF)
ordJp = ordinate(ff1ra, method = "PCoA", distance = DistJ)
orduUFp = ordinate(ff1ra, method = "PCoA", distance = DistuUF)
pBC = plot_ordination(ff1ra, ordBCp, type="samples", color="Time", 
                      shape="Treatment", label="SampleID_Sequencing", title="PCoA - Bray Curtis - No norm.") +
  geom_point(size=2)
pUF = plot_ordination(ff1ra, orduUFp, type="samples", color="Time", 
                      shape="Treatment", label="SampleID_Sequencing", title="PCoA - Unifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=CATEGORY),alpha=.2)
pJ = plot_ordination(ff1ra, ordJp, type="samples", color="Time", 
                     shape="Treatment", label="SampleID_Sequencing", title="PCoA - Jaccard - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=CATEGORY),alpha=.2)
pwUF = plot_ordination(ff1ra, ordwUFp, type="samples", color="Time", 
                       shape="Treatment", label="SampleID_Sequencing", title="PCoA - wUnifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=CATEGORY),alpha=.2)
gridExtra::grid.arrange(pBC,pUF,pJ,pwUF,ncol=2)

ordBCn = ordinate(ff1ra, method = "NMDS", distance = DistBC)
ordwUFn = ordinate(ff1ra, method = "NMDS", distance = DistwUF)
ordJn = ordinate(ff1ra, method = "NMDS", distance = DistJ)
orduUFn = ordinate(ff1ra, method = "NMDS", distance = DistuUF)
nBC = plot_ordination(ff1ra, ordBCn, type="samples", color="Time", 
                      shape="Treatment", label="SampleID_Sequencing", title="NMDS - Bray Curtis - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
nUF = plot_ordination(ff1ra, orduUFn, type="samples", color="Time", 
                      shape="Treatment", label="SampleID_Sequencing", title="NMDS - Unifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
nJ = plot_ordination(ff1ra, ordJn, type="samples", color="Time", 
                     shape="Treatment", label="SampleID_Sequencing", title="NMDS - Jaccard - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
nwUF = plot_ordination(ff1ra, ordwUFn, type="samples", color="Time", 
                       shape="Treatment", label="SampleID_Sequencing", title="NMDS - wUnifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
gridExtra::grid.arrange(nBC,nUF,nJ,nwUF,ncol=2)

rm(DistJ,DistwUF,DistuUF,ordBCp,ordwUFp,ordJp,orduUFp,ordBCn,ordwUFn,ordJn,orduUFn)

# ## 1st log transf.
# pslog <- transform_sample_counts(ff1ra, function(x) log(1 + x))
# DistBC = distance(pslog, method = "bray")
# DistwUF = distance(pslog, method = "wUniFrac")
# DistJ = distance(pslog, method = "jaccard")
# DistuUF = distance(pslog, method = "UniFrac")
# 
# ordBCp = ordinate(pslog, method = "PCoA", distance = DistBC)
# ordwUFp = ordinate(pslog, method = "PCoA", distance = DistwUF)
# ordJp = ordinate(pslog, method = "PCoA", distance = DistJ)
# orduUFp = ordinate(pslog, method = "PCoA", distance = DistuUF)
# 
# pBC = plot_ordination(pslog, ordBCp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Bray Curtis - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pUF = plot_ordination(pslog, orduUFp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Unifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pJ = plot_ordination(pslog, ordJp, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="PCoA - Jaccard - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pwUF = plot_ordination(pslog, ordwUFp, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="PCoA - wUnifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(pBC,pUF,pJ,pwUF,ncol=2)
# 
# ordBCn = ordinate(pslog, method = "NMDS", distance = DistBC)
# ordwUFn = ordinate(pslog, method = "NMDS", distance = DistwUF)
# ordJn = ordinate(pslog, method = "NMDS", distance = DistJ)
# orduUFn = ordinate(pslog, method = "NMDS", distance = DistuUF)
# nBC = plot_ordination(pslog, ordBCn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Bray Curtis - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nUF = plot_ordination(pslog, orduUFn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Unifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nJ = plot_ordination(pslog, ordJn, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="NMDS - Jaccard - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nwUF = plot_ordination(pslog, ordwUFn, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="NMDS - wUnifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(nBC,nUF,nJ,nwUF,ncol=2)
# 
# 
# ##----- 2nd rarefy
# rm(DistBC,DistJ,DistwUF,DistuUF,ordBCp,ordwUFp,ordJp,orduUFp,ordBCn,ordwUFn,ordJn,orduUFn)
# frr = rarefy_even_depth(ff1)
# frr
# DistBC = distance(frr, method = "bray")
# DistwUF = distance(frr, method = "wUniFrac")
# DistJ = distance(frr, method = "jaccard")
# DistuUF = distance(frr, method = "UniFrac")
# 
# ordBCp = ordinate(frr, method = "PCoA", distance = DistBC)
# ordwUFp = ordinate(frr, method = "PCoA", distance = DistwUF)
# ordJp = ordinate(frr, method = "PCoA", distance = DistJ)
# orduUFp = ordinate(frr, method = "PCoA", distance = DistuUF)
# 
# pBC = plot_ordination(frr, ordBCp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Bray Curtis - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pUF = plot_ordination(frr, orduUFp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Unifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pJ = plot_ordination(frr, ordJp, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="PCoA - Jaccard - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pwUF = plot_ordination(frr, ordwUFp, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="PCoA - wUnifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(pBC,pUF,pJ,pwUF,ncol=2)
# 
# ordBCn = ordinate(frr, method = "NMDS", distance = DistBC)
# ordwUFn = ordinate(frr, method = "NMDS", distance = DistwUF)
# ordJn = ordinate(frr, method = "NMDS", distance = DistJ)
# orduUFn = ordinate(frr, method = "NMDS", distance = DistuUF)
# 
# nBC = plot_ordination(frr, ordBCn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Bray Curtis - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nUF = plot_ordination(frr, orduUFn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Unifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nJ = plot_ordination(frr, ordJn, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="NMDS - Jaccard - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nwUF = plot_ordination(frr, ordwUFn, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="NMDS - wUnifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(nBC,nUF,nJ,nwUF,ncol=2)

#####----------------- beta div test
#### On means
### PERMANOVA Bray-Curtis
pslog <- transform_sample_counts(ff1ra, function(x) log(1 + x))
vegan::adonis(DistBC  ~ sample_data(pslog)$Time*sample_data(pslog)$Treatment)

# Call:
#   vegan::adonis(formula = DistBC ~ sample_data(pslog)$Time * sample_data(pslog)$Treatment) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#                                                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# sample_data(pslog)$Time                             4    7.0654 1.76634 10.4440 0.30516  0.001 ***
# sample_data(ff1)$Treatment                          1    0.2341 0.23411  1.3843 0.01011  0.096 .  
# sample_data(pslog)$Time:sample_data(ff1)$Treatment  4    0.6321 0.15803  0.9344 0.02730  0.642    
# Residuals                                          90   15.2213 0.16913         0.65742           
# Total                                              99   23.1528                 1.00000       



##### Aggregate counts by Genus
glom = tax_glom(ff1,taxrank = 'Genus')
## Export
gen_count = as(otu_table(glom), 'matrix')
OTUdf = as.data.frame(gen_count)
head(OTUdf)
OTUdf$otuid = rownames(OTUdf)

## Tax
tx = as(tax_table(ff1)[,6], 'matrix')
tx = as.data.frame(tx)
tx$otuid = rownames(tx)
tx$Genus = gsub('g__','',tx$Genus)

#Final
OTUdfe = merge(tx,OTUdf,by='otuid')
OTUdfe$otuid = NULL
head(OTUdfe)
##Rm unassigned Genus
OTUdfe = OTUdfe[OTUdfe$Genus!=' ',]
dim(OTUdfe)
#[1]  91 101

write.table(OTUdfe, file='./Genus_count_Clarck.tsv', sep='\t',quote=F, row.names = F)

#####------------- Differential families between control and infected
require(structSSI)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###-- ASV level
ff_dds <- phyloseq_to_deseq2(ff1, ~ Treatment + Time + Treatment:Time)

ffMeans = apply(DESeq2::counts(ff_dds), 1, gm_mean)
ffdds = DESeq2::estimateSizeFactors(ff_dds, geoMeans = ffMeans)
#ffdds$Trt <- relevel(ffdds$Trt, ref = "TM")

dds.all =  DESeq2::DESeq(ffdds, fitType="local")
resultsNames(dds.all)
# [1] "Intercept"          "Treatment_S_vs_R"   "Time_24_vs_0"       "Time_43_vs_0"       "Time_92_vs_0"       "Time_132_vs_0"     
# [7] "TreatmentS.Time24"  "TreatmentS.Time43"  "TreatmentS.Time92"  "TreatmentS.Time132"

res.g.primo1 = DESeq2::results(dds.all,name=c('Treatment_S_vs_R'))
alpha = .05
### Results
res.g.primo1 = res.g.primo1[order(res.g.primo1$padj, na.last = NA), ]
sigtab.g.primo1 = res.g.primo1[(res.g.primo1$padj < alpha), ]
sigtab.g.primo1
sigtab.g.primo1 = cbind(as(sigtab.g.primo1, "data.frame"), as(tax_table(ff1)[rownames(sigtab.g.primo1), ], "matrix"))
posigtab.g.primo1 = sigtab.g.primo1[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus", "Species")]

dim(posigtab.g.primo1)
#[1] 26  10

res.g.primo2 = DESeq2::results(dds.all,name=c('TreatmentS.Time92'))
alpha = .05
### Results
res.g.primo2 = res.g.primo2[order(res.g.primo2$padj, na.last = NA), ]
sigtab.g.primo2 = res.g.primo2[(res.g.primo2$padj < alpha), ]
sigtab.g.primo2
sigtab.g.primo2 = cbind(as(sigtab.g.primo2, "data.frame"), as(tax_table(ff1)[rownames(sigtab.g.primo2), ], "matrix"))
posigtab.g.primo2 = sigtab.g.primo2[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus")]

dim(posigtab.g.primo2)
#[1] 25  10
posigtab.g.primo2
#                                  baseMean log2FoldChange lfcSE     padj            Phylum             Class               Order
# 0bec925243d96e8b5625481184f3d5ba    29.16          17.08 1.526 7.56e-26     p__Firmicutes     c__Clostridia    o__Clostridiales
# a774e0c0eb6179153d3774dbca8a8caa    27.47          20.82 2.812 1.05e-10  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 321c8002ae82566c7d55e3adbeab1c43    23.23         -21.32 3.522 7.01e-07     p__Firmicutes     c__Clostridia    o__Clostridiales
# 12800122dabefdcdcae811c867b191cd    22.15         -20.80 3.455 7.01e-07  p__Fibrobacteres  c__Fibrobacteria  o__Fibrobacterales
# 7f73519772de6ef6be106642ffdad700    22.45          17.74 3.433 7.70e-05     p__Firmicutes     c__Clostridia    o__Clostridiales
# f90d0e0ec7ef5bff26b6348cd18fb24a    23.99          18.42 3.678 1.48e-04  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 21205ae97656cc633be7f9f6ce0df895   105.32         -23.07 4.730 2.47e-04  p__Fibrobacteres  c__Fibrobacteria  o__Fibrobacterales
# bda6a27efc639988db56f2e6c4441722    29.08          -5.56 1.152 2.84e-04  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# e0a14426bca3bb8af68555e7566cefd6    49.52          -4.12 0.982 4.96e-03  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 525eb1c79b690d4bc910e2675bedcd01    94.61          14.63 3.719 1.36e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# 42ff27a663b10f2bd6126d3f5a826540    10.81          27.50 7.107 1.61e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# 10044652a1468c62d16762aa7cd84484    18.05           6.32 1.673 2.11e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# b09c33bf7b6445099b6bbc92a4a24978    28.39          17.13 4.609 2.51e-02  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 495029d7f209a698f6d132e6597fb07a    10.03         -19.53 5.415 3.58e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# e4480bd9faa243d4c52e81cb90afc115    16.02         -21.01 5.869 3.69e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# 2b60fa16f41ecdb593e9f68fa98fb50f     7.45         -20.08 5.676 4.06e-02  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
#                                                    Family            Genus
# 0bec925243d96e8b5625481184f3d5ba       f__Ruminococcaceae  g__Ruminococcus
# a774e0c0eb6179153d3774dbca8a8caa                      f__              g__
# 321c8002ae82566c7d55e3adbeab1c43       f__Ruminococcaceae              g__
# 12800122dabefdcdcae811c867b191cd      f__Fibrobacteraceae   g__Fibrobacter
# 7f73519772de6ef6be106642ffdad700       f__Ruminococcaceae              g__
# f90d0e0ec7ef5bff26b6348cd18fb24a    f__Porphyromonadaceae  g__Paludibacter
# 21205ae97656cc633be7f9f6ce0df895      f__Fibrobacteraceae   g__Fibrobacter
# bda6a27efc639988db56f2e6c4441722  f__[Paraprevotellaceae]         g__YRC22
# e0a14426bca3bb8af68555e7566cefd6  f__[Paraprevotellaceae]              g__
# 525eb1c79b690d4bc910e2675bedcd01       f__Ruminococcaceae              g__
# 42ff27a663b10f2bd6126d3f5a826540       f__Ruminococcaceae              g__
# 10044652a1468c62d16762aa7cd84484       f__Ruminococcaceae              g__
# b09c33bf7b6445099b6bbc92a4a24978                      f__              g__
# 495029d7f209a698f6d132e6597fb07a       f__Ruminococcaceae  g__Oscillospira
# e4480bd9faa243d4c52e81cb90afc115       f__Lachnospiraceae              g__
# 2b60fa16f41ecdb593e9f68fa98fb50f                  f__BS11              g__

res.g.primo2 = DESeq2::results(dds.all,name=c('TreatmentS.Time132'))
alpha = .05
### Results
res.g.primo2 = res.g.primo2[order(res.g.primo2$padj, na.last = NA), ]
sigtab.g.primo2 = res.g.primo2[(res.g.primo2$padj < alpha), ]
sigtab.g.primo2
sigtab.g.primo2 = cbind(as(sigtab.g.primo2, "data.frame"), as(tax_table(ff1)[rownames(sigtab.g.primo2), ], "matrix"))
posigtab.g.primo2 = sigtab.g.primo2[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus")]

dim(posigtab.g.primo2)
#[1] 25  10
posigtab.g.primo2
#                                  baseMean log2FoldChange lfcSE     padj            Phylum             Class               Order
# 0bec925243d96e8b5625481184f3d5ba    29.16          17.08 1.526 7.56e-26     p__Firmicutes     c__Clostridia    o__Clostridiales
# a774e0c0eb6179153d3774dbca8a8caa    27.47          20.82 2.812 1.05e-10  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 321c8002ae82566c7d55e3adbeab1c43    23.23         -21.32 3.522 7.01e-07     p__Firmicutes     c__Clostridia    o__Clostridiales
# 12800122dabefdcdcae811c867b191cd    22.15         -20.80 3.455 7.01e-07  p__Fibrobacteres  c__Fibrobacteria  o__Fibrobacterales
# 7f73519772de6ef6be106642ffdad700    22.45          17.74 3.433 7.70e-05     p__Firmicutes     c__Clostridia    o__Clostridiales
# f90d0e0ec7ef5bff26b6348cd18fb24a    23.99          18.42 3.678 1.48e-04  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 21205ae97656cc633be7f9f6ce0df895   105.32         -23.07 4.730 2.47e-04  p__Fibrobacteres  c__Fibrobacteria  o__Fibrobacterales
# bda6a27efc639988db56f2e6c4441722    29.08          -5.56 1.152 2.84e-04  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# e0a14426bca3bb8af68555e7566cefd6    49.52          -4.12 0.982 4.96e-03  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 525eb1c79b690d4bc910e2675bedcd01    94.61          14.63 3.719 1.36e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# 42ff27a663b10f2bd6126d3f5a826540    10.81          27.50 7.107 1.61e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# 10044652a1468c62d16762aa7cd84484    18.05           6.32 1.673 2.11e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# b09c33bf7b6445099b6bbc92a4a24978    28.39          17.13 4.609 2.51e-02  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
# 495029d7f209a698f6d132e6597fb07a    10.03         -19.53 5.415 3.58e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# e4480bd9faa243d4c52e81cb90afc115    16.02         -21.01 5.869 3.69e-02     p__Firmicutes     c__Clostridia    o__Clostridiales
# 2b60fa16f41ecdb593e9f68fa98fb50f     7.45         -20.08 5.676 4.06e-02  p__Bacteroidetes    c__Bacteroidia    o__Bacteroidales
#                                                    Family            Genus
# 0bec925243d96e8b5625481184f3d5ba       f__Ruminococcaceae  g__Ruminococcus
# a774e0c0eb6179153d3774dbca8a8caa                      f__              g__
# 321c8002ae82566c7d55e3adbeab1c43       f__Ruminococcaceae              g__
# 12800122dabefdcdcae811c867b191cd      f__Fibrobacteraceae   g__Fibrobacter
# 7f73519772de6ef6be106642ffdad700       f__Ruminococcaceae              g__
# f90d0e0ec7ef5bff26b6348cd18fb24a    f__Porphyromonadaceae  g__Paludibacter
# 21205ae97656cc633be7f9f6ce0df895      f__Fibrobacteraceae   g__Fibrobacter
# bda6a27efc639988db56f2e6c4441722  f__[Paraprevotellaceae]         g__YRC22
# e0a14426bca3bb8af68555e7566cefd6  f__[Paraprevotellaceae]              g__
# 525eb1c79b690d4bc910e2675bedcd01       f__Ruminococcaceae              g__
# 42ff27a663b10f2bd6126d3f5a826540       f__Ruminococcaceae              g__
# 10044652a1468c62d16762aa7cd84484       f__Ruminococcaceae              g__
# b09c33bf7b6445099b6bbc92a4a24978                      f__              g__
# 495029d7f209a698f6d132e6597fb07a       f__Ruminococcaceae  g__Oscillospira
# e4480bd9faa243d4c52e81cb90afc115       f__Lachnospiraceae              g__
# 2b60fa16f41ecdb593e9f68fa98fb50f                  f__BS11              g__

###-- Genus level > nothing significant
fff = tax_glom(ff1,taxrank = 'Genus')
fff
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 83 taxa and 100 samples ]
# sample_data() Sample Data:       [ 100 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 83 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 83 tips and 82 internal nodes ]

fff_dds <- phyloseq_to_deseq2(fff, ~ Treatment + Time + Treatment:Time)
# calculate geometric means prior to estimate size factors
fffMeans = apply(DESeq2::counts(fff_dds), 1, gm_mean)
fffdds = DESeq2::estimateSizeFactors(fff_dds, geoMeans = fffMeans)

ddsf =  DESeq2::DESeq(fffdds, fitType="local")

### Results
res.f.primo = DESeq2::results(ddsf,name=c('Treatment_S_vs_R'))
res.f.primo = res.f.primo[order(res.f.primo$padj, na.last=NA), ]
sigtab.f.primo = res.f.primo[(res.f.primo$pvalue < alpha), ]
sigtab.f.primo

sigtab.f.primo = cbind(as(sigtab.f.primo, "data.frame"), as(tax_table(fff)[rownames(sigtab.f.primo), ], "matrix"))
posigtab.f.primo = sigtab.f.primo[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family","Genus")]
posigtab.f.primo

dim(posigtab.f.primo)
#[1] 5  9

posigtab.f.primo
#                                  baseMean log2FoldChange lfcSE padj             Phylum                   Class                Order
# 0e3ea307a2c3fa597d5abf697bf06195    37.82         -0.998 0.482    1      p__Firmicutes           c__Clostridia     o__Clostridiales
# 6829b119c3949ae0d6b9daba3215e1d8     9.46          2.766 1.211    1  p__Proteobacteria  c__Alphaproteobacteria              o__RF32
# cb022615eea1347c64eb13600b5a904c    20.64         -0.927 0.442    1  p__Proteobacteria  c__Alphaproteobacteria  o__Sphingomonadales
# 64c2b97f95351285b96d2123fb6a8696   326.97         -1.020 0.443    1   p__Bacteroidetes          c__Bacteroidia     o__Bacteroidales
# 346f0327c4443b74859284a317cf5e64   605.88         -0.515 0.226    1   p__Bacteroidetes          c__Bacteroidia     o__Bacteroidales
# Family             Genus
# 0e3ea307a2c3fa597d5abf697bf06195  f__[Mogibacteriaceae]  g__Mogibacterium
# 6829b119c3949ae0d6b9daba3215e1d8                    f__               g__
# cb022615eea1347c64eb13600b5a904c   f__Sphingomonadaceae   g__Sphingomonas
# 64c2b97f95351285b96d2123fb6a8696      f__Bacteroidaceae          g__BF311
# 346f0327c4443b74859284a317cf5e64                f__RF16               g__

## plot
keepPhyla.f = unique(posigtab.f.primo$Genus[!is.na(posigtab.f.primo)])
tabprimo.f = subset_taxa(fff, Genus %in% keepPhyla.f)

#pdf(file='./output/DifferentialFamilies_D6.pdf')
## plot

plot_abundance3(tabprimo.f,Facet='Genus',Color='Phylum',
                Varia='Treatment',val=c('R','S'))
#dev.off()

## Export
# write.csv(posigtab.o.primo,file='./DifferentialFamiliesBetwTrt_Pup.csv',quote=F,row.names=F)
# write.csv(posigtab.g.primo,file='./DifferentialGenusBetwTrt_Pup.csv',quote=F,row.names=F)


###-- Genus level > nothing significant
fff = tax_glom(ff1,taxrank = 'Family')
fff
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 83 taxa and 100 samples ]
# sample_data() Sample Data:       [ 100 samples by 10 sample variables ]
# tax_table()   Taxonomy Table:    [ 83 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 83 tips and 82 internal nodes ]

fff_dds <- phyloseq_to_deseq2(fff, ~ Treatment + Time + Treatment:Time)
# calculate geometric means prior to estimate size factors
fffMeans = apply(DESeq2::counts(fff_dds), 1, gm_mean)
fffdds = DESeq2::estimateSizeFactors(fff_dds, geoMeans = fffMeans)

ddsf =  DESeq2::DESeq(fffdds, fitType="local")

### Results
res.f.primo = DESeq2::results(ddsf,name=c('Treatment_S_vs_R'))
res.f.primo = res.f.primo[order(res.f.primo$padj, na.last=NA), ]
sigtab.f.primo = res.f.primo[(res.f.primo$padj < alpha), ]
sigtab.f.primo

sigtab.f.primo = cbind(as(sigtab.f.primo, "data.frame"), as(tax_table(fff)[rownames(sigtab.f.primo), ], "matrix"))
posigtab.f.primo = sigtab.f.primo[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family","Genus")]
posigtab.f.primo

dim(posigtab.f.primo)
#[1] 5  9

posigtab.f.primo
#                                  baseMean log2FoldChange lfcSE padj             Phylum                   Class                Order
# 0e3ea307a2c3fa597d5abf697bf06195    37.82         -0.998 0.482    1      p__Firmicutes           c__Clostridia     o__Clostridiales
# 6829b119c3949ae0d6b9daba3215e1d8     9.46          2.766 1.211    1  p__Proteobacteria  c__Alphaproteobacteria              o__RF32
# cb022615eea1347c64eb13600b5a904c    20.64         -0.927 0.442    1  p__Proteobacteria  c__Alphaproteobacteria  o__Sphingomonadales
# 64c2b97f95351285b96d2123fb6a8696   326.97         -1.020 0.443    1   p__Bacteroidetes          c__Bacteroidia     o__Bacteroidales
# 346f0327c4443b74859284a317cf5e64   605.88         -0.515 0.226    1   p__Bacteroidetes          c__Bacteroidia     o__Bacteroidales
#                                                  Family             Genus
# 0e3ea307a2c3fa597d5abf697bf06195  f__[Mogibacteriaceae]  g__Mogibacterium
# 6829b119c3949ae0d6b9daba3215e1d8                    f__               g__
# cb022615eea1347c64eb13600b5a904c   f__Sphingomonadaceae   g__Sphingomonas
# 64c2b97f95351285b96d2123fb6a8696      f__Bacteroidaceae          g__BF311
# 346f0327c4443b74859284a317cf5e64                f__RF16               g__

## plot
keepPhyla.f = unique(posigtab.f.primo$Genus[!is.na(posigtab.f.primo)])
tabprimo.f = subset_taxa(fff, Genus %in% keepPhyla.f)



######============================================================================
#######----------------------------- OTU analysis from Peachey et al. 2019
######============================================================================
rm(list =ls())

## Working on QIIME2 abundance tables
otu_table = read.table("./OTU_Peachey_minfreq5_minsamp2/dada2_output_mxe1_exported/feature-table.tsv",sep="\t",header=T)
colnames(otu_table)=gsub('X','',colnames(otu_table))
dim(otu_table)
#[1] 5233   95

otu_matrix = read.csv("./OTU_Peachey_minfreq5_minsamp2/dada2_output_mxe1_exported/biom-taxonomy.tsv",sep="\t",header=T)
otu_matrix = as.matrix(otu_matrix)

metadata = read.csv("./OTU_Peachey_minfreq5_minsamp2/sample-metadata_v6.csv", sep=",",header=T)
rownames(metadata) = metadata$SampleID_Sequencing
head(metadata)
#   Sample.id Label subject Include Infection Timepoint sex
# 1         1     1     ED1       1         0         a   F
# 2         3     3     BX1       1         0         a   F
# 3         4     4      F1       1         0         a   M
# 4         6     6     EA1       1         0         a   F
# 5         9     9     KD1       1         0         a   M
# 6        10    10     HL1       1         0         a   F

metadata$Sample.id = factor(metadata$Sample.id)
metadata$Infection = factor(metadata$Infection)
metadata$Timepoint = factor(metadata$Timepoint)
metadata$sex = factor(metadata$sex)
rownames(metadata)=metadata$Sample.id
## engineering to retrive true subject id
metadata$ID =  sapply(stringr::str_split(metadata$subject,' '), function(x) x[1])
metadata$ID = substr(metadata$ID, 1, nchar(metadata$ID)-1)

LOW = metadata$ID[metadata$Infection==0 & metadata$Timepoint=='a']
length(LOW)
#[1] 14
HIGH = metadata$ID[metadata$Infection==1 & metadata$Timepoint=='a']
length(HIGH)
#[1] 9

## Data struct
table(metadata$Infection,metadata$Timepoint)
#        a  b  c
# 0     14 15 14
# 1      9  9  9
# 2      7  6  6
# blank  0  0  0

merg = merge(otu_table,otu_matrix,by=c('OTUID'))
dim(merg)
#[1] 5233   97

## Split taxonomy info
merg2 = tidyr::separate(merg,taxonomy,c('Kingdom','Phylum','Class','Order','Family','Genus','Species'),';')
dim(merg2)
#[1] 5233  103

## New tables saved to csv // ALLLL names have to be similar in OTU, TAX ANNDDD phy tree !!!
otu = as.matrix(merg2[,2:(dim(metadata)[1]+1)])
row.names(otu) = merg2$OTUID 
#colnames(otu)=gsub('[.]','-',colnames(otu))

tax = as.matrix(merg2[,(dim(metadata)[1]+2):dim(merg2)[2]])
row.names(tax) = merg2$OTUID 

#####--------------------------------------------- Phyloseq data
OTU = otu_table(as.matrix(otu), taxa_are_rows = T)
TAX = tax_table(as.matrix(tax))
meta = sample_data(metadata)
phy_tree = read_tree("./OTU_Peachey_minfreq5_minsamp2/tree_out/tree_unrooted.1/tree.nwk")
peachey0 = phyloseq(OTU, TAX , meta, phy_tree)
peachey0
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5233 taxa and 94 samples ]
# sample_data() Sample Data:       [ 94 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 5233 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 5233 tips and 5215 internal nodes ]

### Subset to the 23 animals considered in study (infection 0 or 1)

peachey = subset_samples(peachey0, ID %in% c(HIGH,LOW) & Infection %in% c(0,1))
peachey
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5233 taxa and 67 samples ]
# sample_data() Sample Data:       [ 67 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 5233 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 5233 tips and 5215 internal nodes ]

#####===================================================================================================
#### PRE-PROCESSING & DATA EXPLORATION
#####===================================================================================================

#####--------------------------------------------- Seq Depth
sdt = data.table(as(sample_data(peachey), "data.frame"),
                 TotalReads = sample_sums(peachey), keep.rownames = TRUE)
#setnames(sdt, "rn", "sample.id")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth") + 
  annotation_logticks(scaled = T,sides="b") +
  scale_x_log10(breaks = c(30000,50000,100000,200000), 
                labels = c(30000,50000,100000,200000)) 

pSeqDepth + facet_grid( Infection ~ Timepoint)

aggregate(TotalReads ~ Infection + Timepoint, data= sdt, FUN=mean)
#   Infection Timepoint TotalReads
# 1         0         a      39412
# 2         1         a      34164
# 3         0         b      39202
# 4         1         b      34997
# 5         0         c      35592
# 6         1         c      37718

## no bias
summary.aov(lm(TotalReads ~ Infection + Timepoint,data=sdt))
#             Df   Sum Sq  Mean Sq F value Pr(>F)
# Infection    1 1.06e+08 1.06e+08    2.58   0.11
# Timepoint    2 1.11e+07 5.55e+06    0.13   0.87
# Residuals   63 2.59e+09 4.11e+07         

### Need to normalize sequencing depth across samples
summary(sdt$TotalReads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18561   33783   36942   37158   41358   51337 

#####--------------------------------------------- Taxa filtering based on prevalence and total counts
### Total counts
tdt = data.table(tax_table(peachey),
                 TotalCounts = taxa_sums(peachey),
                 OTU = taxa_names(peachey))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts")

# How many unknown Families ? 
length(which(is.na(tdt$Family)))
#[1] 163
# How many unknown Genus ? 
length(which(is.na(tdt$Genus)))
#[1] 430

# How many singletons?
tdt[(TotalCounts <= 0), .N]
#[1] 141
tdt[(TotalCounts <= 1), .N]
#[1] 141

# Min count = 5 
tdt[(TotalCounts < 5), .N]
# [1] 191

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
#- zoom-in
pCumSum + xlim(0, 20000)

###---------------------- Taxa prevalence

## Quicker (from evomics)
source('~/Documents/scripts/taxa_summary.R',local=TRUE)
mdt = fast_melt(peachey)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]

ggplot(prevdt, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Histogram of Taxa Prevalence")

## Singletons ?
prevdt[(Prevalence <= 1), .N]
#[1] 604
# How many doubletons?
prevdt[(Prevalence <= 2), .N]
# 1688

prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pPrevCumSum

### Prevalence vs. Total count
ggplot(prevdt, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

## Add Phylum information
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))

# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt2 <- addPhylum[prevdt]
showPhyla = prevdt2[, sum(TotalCounts), by = Phylum][order(-V1)][1:6]$Phylum ## only 6 phyla observed ??
setkey(prevdt2, Phylum)
ggplot(prevdt2[showPhyla], 
       mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_point(size = 4, alpha = 0.75) + 
  geom_vline(xintercept = 4) + geom_hline(yintercept = 5) +
  scale_y_log10()

### Define prevalence threshold as 5% of remaining samples
prevalenceThreshold = 0.05 * nsamples(peachey)
prevalenceThreshold
#[1] 3.4

### Filter entries with unidentified Phylum.
table(addPhylum$Phylum)
# # p__   p__Actinobacteria    p__Bacteroidetes    p__Cyanobacteria  p__Deferribacteres 
# # 18                  33                1857                   1                   1 
# # p__Elusimicrobia    p__Fibrobacteres       p__Firmicutes     p__Fusobacteria    p__Lentisphaerae 
# # 2                  41                2618                   1                   3 
# # p__Planctomycetes   p__Proteobacteria     p__Spirochaetes      p__Tenericutes              p__TM7 
# # 2                 119                 233                  64                   8 
# # p__Verrucomicrobia 
# # 218 
keepPhyla = table(addPhylum$Phylum)[-1]
keepPhyla
# p__Actinobacteria    p__Bacteroidetes    p__Cyanobacteria  p__Deferribacteres    p__Elusimicrobia 
#                33                1857                   1                   1                   2 
# p__Fibrobacteres       p__Firmicutes     p__Fusobacteria    p__Lentisphaerae   p__Planctomycetes 
#               41                2618                   1                   3                   2 
# p__Proteobacteria     p__Spirochaetes      p__Tenericutes              p__TM7  p__Verrucomicrobia 
#               119                 233                  64                   8                 218 

ff1 = subset_taxa(peachey, Phylum %in% names(keepPhyla))
ff1
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5201 taxa and 67 samples ]
# sample_data() Sample Data:       [ 67 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 5201 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 5201 tips and 5183 internal nodes ]

## Transform to relative abundance
ff1ra = transform_sample_counts(ff1, function(x){x / sum(x)})

sdt2 = data.table(as(sample_data(ff1), "data.frame"),
                  TotalReads = sample_sums(ff1), keep.rownames = TRUE)


#####--------------------------------------------- Abundances 
### plot rarefaction curve before Genus agglomeration
rc = vegan::rarecurve(t(otu_table(ff1)), step=50, cex=0.5,label=T)

rare <- lapply(rc, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

rare <- purrr::map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

rare$sample <- levels(sample_data(ff1)$Sample.id)[as.integer(rare$sample)]
rare$Infection = sample_data(ff1)$Infection[match(rare$sample,sample_data(ff1)$Sample.id)]

ggplot(data = rare,aes(x = raw.read, y = OTU, Day=sample,color=Infection))+
  geom_line() +
  scale_x_continuous(labels =  scales::scientific_format()) +
  ggtitle('Rarefaction curves')

#####===================================================================================================
#### Pre-processing & main observations
#####===================================================================================================

#####--------------------------------------------- Alpha diversity (not to be estimated on filtered data)
## Alpha div estimators account for variation in library sizes
## Shannon influenced by low abundance species
## InvSimpson influenced by major species

pAlpha2 = plot_richness(ff1,
                        shape = "Infection",
                        color = "Infection",
                        x = 'Timepoint',
                        measures = c("Shannon"),#c("Observed", "Chao1", "Shannon", "InvSimpson"),
                        title = "Alpha Diversity - Filtering minor OTUs")
pAlpha2 + geom_boxplot() #size = 3)

## By sample
pAlpha3 = plot_richness(ff1,
                        shape = "Infection",
                        color = "Infection",
                        x = 'Sample.id',
                        measures = c("Shannon"),
                        title = "Alpha Diversity - Filtering minor OTUs")
pAlpha3 + geom_point(size = 3)

###------------------------------------ Linear modeling of alpha diversity

####-------------- Unagglom. dataset

#### Shannon
ff1_alpha_div <- estimate_richness(ff1, split = TRUE, measure = "Shannon")
ff1_alpha_div$Sample.id <- rownames(ff1_alpha_div) %>%  as.factor()
ff1_samp <- sample_data(ff1) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ff1_alpha_div, by = "Sample.id") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")

# reorder's facet from lowest to highest diversity
diversity_means <- ff1_samp %>%
  group_by(Infection,Timepoint) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
# ff1_samp$Stage <- factor(ff1_samp$Stage,
#                          diversity_means$Stage)

alpha_div_model <- lm( alpha_diversity ~ Infection*Timepoint , data = ff1_samp)
summary.aov(alpha_div_model)
#                     Df Sum Sq Mean Sq F value  Pr(>F)    
# Infection            1   0.05   0.050    0.39 0.53425    
# Timepoint            2   2.20   1.099    8.52 0.00054 ***
# Infection:Timepoint  2   0.31   0.157    1.22 0.30270    
# Residuals           61   7.87   0.129      

summary(alpha_div_model)
# Call:
#   lm(formula = alpha_diversity ~ Infection * Timepoint, data = ff1_samp)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.1062 -0.1987  0.0376  0.2112  0.7006 
# 
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              5.431      0.096   56.59   <2e-16 ***
# Infection1              -0.250      0.153   -1.63    0.108    
# Timepointb               0.315      0.136    2.32    0.024 *  
# Timepointc               0.131      0.141    0.93    0.358    
# Infection1:Timepointb    0.309      0.217    1.43    0.159    
# Infection1:Timepointc    0.277      0.221    1.25    0.214    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.359 on 61 degrees of freedom
# Multiple R-squared:  0.246,	Adjusted R-squared:  0.184 
# F-statistic: 3.97 on 5 and 61 DF,  p-value: 0.00347 

# predictions and confidence intervals.
new_data <- expand.grid(Timepoint = levels(ff1_samp$Timepoint), Infection = levels(ff1_samp$Infection))
a = predict(alpha_div_model, newdata = new_data, interval="confidence")
a = data.frame(a)
new_data$pred = a$fit
new_data$lwr= a$lwr
new_data$upr= a$upr

# X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
#                   new_data[-ncol(new_data)])
# pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
# new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2

# fitted values, with error bars
shannon_plot=ggplot(ff1_samp %>% left_join(new_data)) +
  facet_wrap(~ Infection) +
  geom_errorbar(aes(x = Timepoint, ymin = lwr,
                    ymax = upr),
                col = "#858585", size = .5,width=.2) +
  geom_point(aes(x = Timepoint, y = alpha_diversity,
                 col = Infection), size = 4) +
  scale_y_continuous(limits = c(5, 7), breaks = seq(5, 7, .5)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Stage", y = "Shannon Diversity", color = "Stage") +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90),legend.position = 'bottom') 
shannon_plot

#####--------------------------------------------- Beta diversity

#####------------------ ALL OTUs: PCoA 
### UniFrac accounts for phylogenetic d
### Sequencing depth can greatly affects ordination
### Read counts should be normalized: either with DESeq2 or applying log transf.

## No transf.
DistBC = distance(ff1ra, method = "bray")
DistwUF = distance(ff1ra, method = "wUniFrac")
DistJ = distance(ff1ra, method = "jaccard")
DistuUF = distance(ff1ra, method = "UniFrac")

ordBCp = ordinate(ff1ra, method = "PCoA", distance = DistBC)
ordwUFp = ordinate(ff1ra, method = "PCoA", distance = DistwUF)
ordJp = ordinate(ff1ra, method = "PCoA", distance = DistJ)
orduUFp = ordinate(ff1ra, method = "PCoA", distance = DistuUF)
pBC = plot_ordination(ff1ra, ordBCp, type="samples", color="Timepoint", 
                      shape="Infection", label="Sample.id", title="PCoA - Bray Curtis - No norm.") +
  geom_point(size=2)
pUF = plot_ordination(ff1ra, orduUFp, type="samples", color="Timepoint", 
                      shape="Infection", label="Sample.id", title="PCoA - Unifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=CATEGORY),alpha=.2)
pJ = plot_ordination(ff1ra, ordJp, type="samples", color="Timepoint", 
                     shape="Infection", label="Sample.id", title="PCoA - Jaccard - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=CATEGORY),alpha=.2)
pwUF = plot_ordination(ff1ra, ordwUFp, type="samples", color="Timepoint", 
                       shape="Infection", label="Sample.id", title="PCoA - wUnifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=CATEGORY),alpha=.2)
gridExtra::grid.arrange(pBC,pUF,pJ,pwUF,ncol=2)

ordBCn = ordinate(ff1ra, method = "NMDS", distance = DistBC)
ordwUFn = ordinate(ff1ra, method = "NMDS", distance = DistwUF)
ordJn = ordinate(ff1ra, method = "NMDS", distance = DistJ)
orduUFn = ordinate(ff1ra, method = "NMDS", distance = DistuUF)
nBC = plot_ordination(ff1ra, ordBCn, type="samples", color="Timepoint", 
                      shape="Infection", label="Sample.id", title="NMDS - Bray Curtis - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
nUF = plot_ordination(ff1ra, orduUFn, type="samples", color="Timepoint", 
                      shape="Infection", label="Sample.id", title="NMDS - Unifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
nJ = plot_ordination(ff1ra, ordJn, type="samples", color="Timepoint", 
                     shape="Infection", label="Sample.id", title="NMDS - Jaccard - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
nwUF = plot_ordination(ff1ra, ordwUFn, type="samples", color="Timepoint", 
                       shape="Infection", label="Sample.id", title="NMDS - wUnifrac - No norm.") +
  geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
gridExtra::grid.arrange(nBC,nUF,nJ,nwUF,ncol=2)

rm(DistJ,DistwUF,DistuUF,ordBCp,ordwUFp,ordJp,orduUFp,ordBCn,ordwUFn,ordJn,orduUFn)

# ## 1st log transf.
# pslog <- transform_sample_counts(ff1ra, function(x) log(1 + x))
# DistBC = distance(pslog, method = "bray")
# DistwUF = distance(pslog, method = "wUniFrac")
# DistJ = distance(pslog, method = "jaccard")
# DistuUF = distance(pslog, method = "UniFrac")
# 
# ordBCp = ordinate(pslog, method = "PCoA", distance = DistBC)
# ordwUFp = ordinate(pslog, method = "PCoA", distance = DistwUF)
# ordJp = ordinate(pslog, method = "PCoA", distance = DistJ)
# orduUFp = ordinate(pslog, method = "PCoA", distance = DistuUF)
# 
# pBC = plot_ordination(pslog, ordBCp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Bray Curtis - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pUF = plot_ordination(pslog, orduUFp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Unifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pJ = plot_ordination(pslog, ordJp, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="PCoA - Jaccard - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pwUF = plot_ordination(pslog, ordwUFp, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="PCoA - wUnifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(pBC,pUF,pJ,pwUF,ncol=2)
# 
# ordBCn = ordinate(pslog, method = "NMDS", distance = DistBC)
# ordwUFn = ordinate(pslog, method = "NMDS", distance = DistwUF)
# ordJn = ordinate(pslog, method = "NMDS", distance = DistJ)
# orduUFn = ordinate(pslog, method = "NMDS", distance = DistuUF)
# nBC = plot_ordination(pslog, ordBCn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Bray Curtis - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nUF = plot_ordination(pslog, orduUFn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Unifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nJ = plot_ordination(pslog, ordJn, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="NMDS - Jaccard - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nwUF = plot_ordination(pslog, ordwUFn, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="NMDS - wUnifrac - Log norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(nBC,nUF,nJ,nwUF,ncol=2)
# 
# 
# ##----- 2nd rarefy
# rm(DistBC,DistJ,DistwUF,DistuUF,ordBCp,ordwUFp,ordJp,orduUFp,ordBCn,ordwUFn,ordJn,orduUFn)
# frr = rarefy_even_depth(ff1)
# frr
# DistBC = distance(frr, method = "bray")
# DistwUF = distance(frr, method = "wUniFrac")
# DistJ = distance(frr, method = "jaccard")
# DistuUF = distance(frr, method = "UniFrac")
# 
# ordBCp = ordinate(frr, method = "PCoA", distance = DistBC)
# ordwUFp = ordinate(frr, method = "PCoA", distance = DistwUF)
# ordJp = ordinate(frr, method = "PCoA", distance = DistJ)
# orduUFp = ordinate(frr, method = "PCoA", distance = DistuUF)
# 
# pBC = plot_ordination(frr, ordBCp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Bray Curtis - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pUF = plot_ordination(frr, orduUFp, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="PCoA - Unifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pJ = plot_ordination(frr, ordJp, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="PCoA - Jaccard - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# pwUF = plot_ordination(frr, ordwUFp, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="PCoA - wUnifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(pBC,pUF,pJ,pwUF,ncol=2)
# 
# ordBCn = ordinate(frr, method = "NMDS", distance = DistBC)
# ordwUFn = ordinate(frr, method = "NMDS", distance = DistwUF)
# ordJn = ordinate(frr, method = "NMDS", distance = DistJ)
# orduUFn = ordinate(frr, method = "NMDS", distance = DistuUF)
# 
# nBC = plot_ordination(frr, ordBCn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Bray Curtis - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nUF = plot_ordination(frr, orduUFn, type="samples", color="Day", 
#                       shape="Group", label="Indiv", title="NMDS - Unifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nJ = plot_ordination(frr, ordJn, type="samples", color="Day", 
#                      shape="Group", label="Indiv", title="NMDS - Jaccard - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# nwUF = plot_ordination(frr, ordwUFn, type="samples", color="Day", 
#                        shape="Group", label="Indiv", title="NMDS - wUnifrac - Rarefaction norm.") +
#   geom_point(size=2) #+ geom_polygon(aes(fill=Day),alpha=.2)
# gridExtra::grid.arrange(nBC,nUF,nJ,nwUF,ncol=2)

#####----------------- beta div test
#### On means
### PERMANOVA Bray-Curtis
pslog <- transform_sample_counts(ff1ra, function(x) log(1 + x))
vegan::adonis(DistBC  ~ sample_data(pslog)$Timepoint*sample_data(pslog)$Infection)

#                                                           Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)    
# sample_data(pslog)$Timepoint                               2      1.43   0.713   2.747 0.077  0.001 ***
# sample_data(pslog)$Infection                               1      0.71   0.712   2.743 0.039  0.001 ***
# sample_data(pslog)$Timepoint:sample_data(pslog)$Infection  2      0.51   0.255   0.981 0.028  0.533    
# Residuals                                                 61     15.83   0.260         0.857           
# Total                                                     66     18.48                 1.000            

#####------------- Differential families between control and infected
require(structSSI)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###------------------ Export Genus count table

## Aggregate counts by Genus
glom = tax_glom(ff1,taxrank = 'Genus')
## Export
gen_count = as(otu_table(glom), 'matrix')
OTUdf = as.data.frame(gen_count)
head(OTUdf)
OTUdf$otuid = rownames(OTUdf)

## Tax
tx = as(tax_table(ff1)[,6], 'matrix')
tx = as.data.frame(tx)
tx$otuid = rownames(tx)
tx$Genus = gsub('g__','',tx$Genus)

#Final
OTUdfe = merge(tx,OTUdf,by='otuid')
OTUdfe$otuid = NULL
head(OTUdfe)
##Rm unassigned Genus
OTUdfe = OTUdfe[OTUdfe$Genus!=' ',]
dim(OTUdfe)
#[1]  61 68
write.table(OTUdfe, file='./Genus_count_Peachey.tsv', sep='\t',quote=F, row.names = F)

###-------------------------- ASV level
ff_dds <- phyloseq_to_deseq2(ff1, ~ Infection + Timepoint + Infection:Timepoint)

ffMeans = apply(DESeq2::counts(ff_dds), 1, gm_mean)
ffdds = DESeq2::estimateSizeFactors(ff_dds, geoMeans = ffMeans)
#ffdds$Trt <- relevel(ffdds$Trt, ref = "TM")

dds.all =  DESeq2::DESeq(ffdds, fitType="local")
resultsNames(dds.all)
# [1] "Intercept"             "Infection_1_vs_0"      "Timepoint_b_vs_a"      "Timepoint_c_vs_a"      "Infection1.Timepointb"
# [6] "Infection1.Timepointc"

res.g.primo = DESeq2::results(dds.all,name=c('Infection_1_vs_0'))
alpha = 0.05
### Results
res.g.primo = res.g.primo[order(res.g.primo$padj, na.last=NA), ]
sigtab.g.primo = res.g.primo[(res.g.primo$padj < alpha), ]
sigtab.g.primo

sigtab.g.primo = cbind(as(sigtab.g.primo, "data.frame"), as(tax_table(ff1)[rownames(sigtab.g.primo), ], "matrix"))
posigtab.g.primo = sigtab.g.primo[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family","Genus",'Species')]

dim(posigtab.g.primo)
#[1] 63 10

###----------------------------------- Genus level
fff = tax_glom(ff1,taxrank = 'Genus')
fff
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 103 taxa and 67 samples ]
# sample_data() Sample Data:       [ 67 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 103 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 103 tips and 101 internal nodes ]

fff_dds <- phyloseq_to_deseq2(fff, ~ Infection + Timepoint + Infection:Timepoint)
# calculate geometric means prior to estimate size factors
fffMeans = apply(DESeq2::counts(fff_dds), 1, gm_mean)
fffdds = DESeq2::estimateSizeFactors(fff_dds, geoMeans = fffMeans)

ddsf =  DESeq2::DESeq(fffdds, fitType="local")

### Results
res.f.primo = DESeq2::results(ddsf,name=c('Infection_1_vs_0'))
res.f.primo = res.f.primo[order(res.f.primo$padj, na.last=NA), ]
sigtab.f.primo = res.f.primo[(res.f.primo$padj < alpha), ]
sigtab.f.primo

sigtab.f.primo = cbind(as(sigtab.f.primo, "data.frame"), as(tax_table(fff)[rownames(sigtab.f.primo), ], "matrix"))
posigtab.f.primo = sigtab.f.primo[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family","Genus")]

dim(posigtab.f.primo)
#[1] 2 9 

posigtab.f.primo
#                                  baseMean log2FoldChange lfcSE     padj            Phylum           Class             Order            Family
# 0a94fec23e738748a57d865e9a033993     11.9          -7.00  1.55 0.000631    p__Tenericutes          c__RF3      o__ML615J-28               f__
# c4db80a3b68faec7cc21146f1f416afd     24.6          -5.05  1.47 0.030331  p__Bacteroidetes  c__Bacteroidia  o__Bacteroidales  f__Rikenellaceae

## plot
keepPhyla.f = unique(posigtab.f.primo$Family[!is.na(posigtab.f.primo)])
tabprimo.f = subset_taxa(fff, Family %in% keepPhyla.f)

# #pdf(file='./output/DifferentialFamilies_D6.pdf')
# plot_abundance3 = function(physeq, ylabn = "",
#                            Facet = "Genus",
#                            Color = "Genus",Varia = 'Group',val = 'D6_CP',NCol=3){
#   mphyseq = psmelt(physeq)
#   mphyseq <- subset(mphyseq, Abundance > 0)
#   ggplot(data = mphyseq[match.fun("%in%")(mphyseq[[Varia]], val), ],
#          mapping = aes_string(x = Varia, y = "Abundance",
#                               color = Color, fill = Color)) +
#     geom_violin(fill = NA) +
#     geom_point(size = 2, alpha = 0.3,
#                position = position_jitter(width = 0.3)) +
#     facet_wrap(facets = Facet, ncol=NCol) + ylab(ylabn) +
#     scale_y_log10() +
#     theme(text = element_text(size=16),strip.text = element_text(size = 14),
#           axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position='bottom')
# }
# plot_abundance3(tabprimo.f,Facet='Family',Color='Phylum',Varia='Infection',val=c('LOW','HIGH'))
# #dev.off()

## Export
# write.csv(posigtab.o.primo,file='./DifferentialFamiliesBetwTrt_Pup.csv',quote=F,row.names=F)
# write.csv(posigtab.g.primo,file='./DifferentialGenusBetwTrt_Pup.csv',quote=F,row.names=F)

#   labs(col = "Group", x = "Axis1", y = "Axis2")
