#' There are two series in the GSE73518 dataset:
#' GPL13534	Illumina HumanMethylation450 BeadChip (HumanMethylation450_15017482)
#' GPL16876	Agilent-020382 Human Custom Microarray 44k (Feature Number version)
#' They comapred methylation info with Chip-seq hg19, so I guess they also had methylation in hg19.

#' Directory with data
data = './'

#' Download information about data as a data.frame and save it as a RDS object
library(GEOquery) #http://genomicsclass.github.io/book/pages/GEOquery.html
options('download.file.method'='curl')
gse <- getGEO("GSE73518", GSEMatrix = TRUE) # https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73518/matrix/

# 
# # GSE73518-GPL13534_series
GPL13534=pData(phenoData(gse[[1]]))
# # GSE73518-GPL16876_series
# GPL16876=pData(phenoData(gse[[2]]))
# 
GPL13534$geo_accession = as.character(GPL13534$geo_accession)
# Convert data.frame columns from factors to characters
GPL13534[] <- lapply(GPL13534, as.character)

#' 
# Take only some columns of the original data.frame
characteristics=c(
  'characteristics_ch1',
  'characteristics_ch1.1',
  'characteristics_ch1.2',
  'characteristics_ch1.3',
  'characteristics_ch1.4',
  'characteristics_ch1.5',
  'characteristics_ch1.6',
  'characteristics_ch1.7'
)
GPL13534.subcl = GPL13534[,c('geo_accession',
                             'source_name_ch1',
                             'supplementary_file',
                             'supplementary_file.1',
                             characteristics[1:6])]

GPL13534.subcl$stage = GPL13534.subcl[,5]
GPL13534.subcl$risk = GPL13534.subcl[,6]
GPL13534.subcl$age_diagnosis = GPL13534.subcl[,7]
GPL13534.subcl$mycn_status = GPL13534.subcl[,8]
GPL13534.subcl$status_1p = GPL13534.subcl[,9]
GPL13534.subcl$status_11q = GPL13534.subcl[,10]
GPL13534.subcl$file1 = basename(GPL13534.subcl$supplementary_file)
GPL13534.subcl$file2 =basename(GPL13534.subcl$supplementary_file.1)
GPL13534.subcl$name = sapply(GPL13534.subcl$file1, function(x) paste(strsplit( x, '_')[[1]][1:3], collapse="_") )
GPL13534.subcl$id = 1:nrow(GPL13534.subcl)

NB.vars = c('mycn_status',
            'stage',
            'risk',
            'age_diagnosis',
            'status_1p',
            'status_11q')

process_columns = function(y){
  sapply( y,
          function(x) substring(strsplit(x, ':')[[1]][2], 2) )
}


newvars_list = lapply(NB.vars, function(feat){
  
  process_columns( GPL13534.subcl[,feat] )
  
} )
newvars_cols = do.call("cbind", newvars_list)
colnames(newvars_cols) = NB.vars
rownames(newvars_cols) <- NULL

GPL13534.subcl[,NB.vars] <- NULL
GPL13534.subcl[,NB.vars] <- newvars_cols

#' Read files in the idat format
# https://github.com/crazyhottommy/DNA-methylation-analysis
files.idat = '/fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Data/AccessoryData/Henrichetal2016/GSE73518/idat/'
#'
library(minfi)
idatFiles <- list.files(files.idat, pattern = "idat.gz$", full = TRUE)
#sapply(idatFiles, gunzip, overwrite = TRUE)
rgSet <- read.metharray.exp(files.idat)

# Calculate beta values and do
# subset-Quantile Within Array Normalisation  (SWAN) For Illumina Infinium HumanMethylation450 BeadChips
library(missMethyl)
# # Subset-quantile within array normalization (SWAN)
mSet <- preprocessRaw(rgSet)
set.seed(12345) #ensure that the normalized intensities will be reproducible; for preprocessSWAN()
mSetSw <- SWAN(mSet,verbose=TRUE)
 
# Filter out poor quality probes
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]
 
# Extracting Beta and M-values
meth <- getMeth(mSetSw) 
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSw)
 
NB.Mval.beta = list(beta=beta,
                     Mval=Mval,
                     mset=mSetSw)
# beta vales
NB.beta = NB.Mval.beta$beta
# M values
NB.Mval = NB.Mval.beta$Mval
# a MethylSet object from info package
NB.mset = NB.Mval.beta$mset

GPL13534.subcl.vars = GPL13534.subcl[,NB.vars]

# Make sure that table with info about feature and a matrix with
# methylation values have the same names of rows
NB.beta.rownames = sapply(colnames(NB.beta), function(x) strsplit(x, "_")[[1]][1]  )




### Load gene expression matched to methylation data
expr.arrays.path="/data/akalin/Projects/BIH_Neuroblastoma/Project/Data/AccessoryData/Henrichetal2016/GSE73517/"
#idatFiles <- list.files(paste0("/idat/",expr.arrays.path), pattern = "txt.gz$", full = TRUE)
#library( GEOquery)
#sapply(idatFiles, gunzip, overwrite = TRUE)
#idatFiles <- list.files(expr.arrays.path, pattern = "txt", full = TRUE)
# Labeling and hybridization was performed following the manufacturer’s protocol and 
# data were background corrected using normexp method saddle-point approximation to maximum likelihood (3) 
# and normalized using the quantile algorithm from limma (4).
descr.expr.arrrays = read.table(paste0(expr.arrays.path, "descr.txt"), sep="\t")
descr.expr.arrrays = descr.expr.arrrays[[1]]
names(descr.expr.arrrays) = paste0("sample",1:length(descr.expr.arrrays))

GSE73517_quantile = read.table(paste0(expr.arrays.path,"GSE73517_quantile.csv"), 
                               sep=",", header=TRUE, row.names=1)
d=data.table::fread(paste0(expr.arrays.path,"A-GEOD-16876_comments.txt"), header=TRUE)
d2=mclapply(rownames(GSE73517_quantile), function(x) which(grepl(paste0(x, "$"), d[["Comment[ProbeName]"]])),
            mc.cores=10)
match.indx=unlist(d2)
agilent020382.matched=d[match.indx,]
agilent020382.matched[["Comment[ChrLoc]"]] = as.numeric(agilent020382.matched[["Comment[ChrLoc]"]])
agilent020382.matched$gene = agilent020382.matched[["Comment[GeneSymbol]"]]
gene.expr = GSE73517_quantile
gene.expr$gene=agilent020382.matched$gene
a=aggregate(gene.expr, list(gene.expr$gene), mean)
rownames(a) = a$Group.1; a$Group.1<-NULL; a$gene<-NULL
gene.expr = a
colnames(gene.expr) = rownames(GPL13534.subcl)




# Annotate probes with genomic locations
# get the 450k annotation data
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

matched.indx=match( rownames(NB.beta), rownames(ann450k))
ann450k.matched = ann450k[matched.indx,]
ann450k.nb = ann450k.matched[,c('chr'     ,  'pos'      ,'pos'    ,    'Name')]
colnames(ann450k.nb) = c("chr", "start", "end", "Name")


# Calculate differentially methylated regions
# 1. between high risk group with MYCN amplification vs intermediate and low risk group of patients
# 2. between high risk group without MYCN amplification vs intermediate and low risk group of patients

mna.indx=which(GPL13534.subcl$risk=="high-risk" & GPL13534.subcl$mycn_status=="MYCN-amplified") 
hr_nmna.indx=which(GPL13534.subcl$risk=="high-risk" & GPL13534.subcl$mycn_status=="MYCN-nonamplified")
control.indx=which(GPL13534.subcl$risk %in% c("intermediate-risk","low-risk"))

#rownames(GPL13534.subcl) == sapply(strsplit(colnames(NB.beta), "_"), function(x) x[[1]][1])
GPL13534.subcl$status=NA
GPL13534.subcl$status[mna.indx]="MNA"
GPL13534.subcl$status[hr_nmna.indx]="HR_nMNA"
GPL13534.subcl$status[control.indx]="control"

require(DMRcate)
dmrs.dmrcate.1 = mclapply(c("MNA", "HR_nMNA"), function(x){

  indx.1 = which(GPL13534.subcl$status %in% c(x,"control" ))
  GPL13534.subcl.1=GPL13534.subcl[indx.1,]
  NB.beta.1=NB.beta[,indx.1]

  group <- factor(GPL13534.subcl.1$status,levels=c(x,"control"))
  design <- model.matrix(~group)

  myannotation <- cpg.annotate("array", NB.beta.1, arraytype = "450K",
                             what="Beta",
                             analysis.type="differential",
                             design=design, coef=2)
  dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, min.cpgs=3)
  results.ranges <- DMRcate::extractRanges(dmrcoutput, genome = "hg19")
  results.ranges

}, mc.cores=2)
#dmrs.dmrcate.1=list(results.ranges.mna,results.ranges )
names(dmrs.dmrcate.1) = c("MNA", "HR_nMNA")
dmrs.dmrcate=lapply(dmrs.dmrcate.1, function(x) x[abs(x$meandiff)>.25])


### Load enhancers from NB cell lines
boevadir="/data/akalin/Projects/BIH_Neuroblastoma/Project/Data/AccessoryData/Boevaetal2017/"
superenhBoeva = readGeneric(paste0(boevadir, 'ng.3921-S3_LiftOverhg38.csv'),
                            sep = ",", 
                            keep.all.metadata = TRUE, 
                            skip=1)

h=c(
  #'Chromosome',
  #'Start',
  #'End',
  'Tmp1',
  'Tmp2',
  'Type (NB or hNCC)',
  'Gene(s)',
  'TF(s)',
  'Number of NB cell lines',
  '#cell lines with SE (group I)',
  '#cell lines with SE (group II)',
  'Median Score over 25 cell lines',
  'Median score group I',
  'Median score group II',
  'Median score hNCC',
  'FC Score Group I over Group II',
  'Wilcoxon p-value (two sided test)',
  'Loadings to PC1',
  'Loadings to PC2',
  'Has 15% overlap with hNCC SEs',
  'Pearson Corr. with gene expression')
colnames(mcols(superenhBoeva)) = h
superenhBoeva.nb = superenhBoeva[as.character(superenhBoeva$'Type (NB or hNCC)')== 'NB']
superenhBoeva.hncc = superenhBoeva[as.character(superenhBoeva$'Type (NB or hNCC)')== 'NCC']


# Rad acet. marks to call enhancers by myself.
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE90683&format=file")
library(utils)
untar("~/GSE90683_RAW.tar")
boevapeaks = list.files(path = "./GSE90683_RAW/", ".K27ac.rep3_regions.bed.gz", full.names = TRUE)

ncc = boevapeaks[grep("NCC", boevapeaks, fixed=T) ]
pdx = boevapeaks[grep('GSM2664372|GSM2664373|GSM2664374|GSM2664375|GSM2664376|GSM2664377', boevapeaks) ]
celllines = setdiff(boevapeaks, c(ncc, pdx))

boeva.grl = lapply(list(celllines=celllines, ncc=ncc, pdx=pdx), function(x){
  mlist = lapply(x, readGeneric)
  names(mlist) = basename(x)
  mlist
} )


MYCN.ampl.cl=c("IGR-NB8","IGR-NB835","NB011","MAP-GR-A99-NB-1","MAP-GR-B25-NB-1",
               "SNJB6","SNJB8", "CLB-CAR", "CLBPE", "IMR32", "LAN1", "N206",
               "SKNBE2C", "SKNDZ","TR14","CLBBER", "CLBMA", "NB69", "CHP212")
noMYCN.ampl.cl=c("MAP-IC-A23-NB-1","CLBGA","SJNB1","SKNFI","SHSY5Y","NBEB",
                 "SKNAS","SKNSH", "SJNB12","GICAN", "SHEP")

boeva.grl$pdx.mycnamp = boeva.grl$pdx[grep("MAP-IC-A23-NB-1", names(boeva.grl$pdx))]
boeva.grl$pdx.nomycnamp = boeva.grl$pdx[-grep("MAP-IC-A23-NB-1", names(boeva.grl$pdx))]
boeva.grl$celllines.mycnamp = boeva.grl$celllines[unlist(sapply(MYCN.ampl.cl, 
                                                                function(x) grep(x, names(boeva.grl$celllines))))]
boeva.grl$celllines.nomycnamp =boeva.grl$celllines[unlist(sapply(noMYCN.ampl.cl, 
                                                                 function(x) grep(x, names(boeva.grl$celllines))))]



# take union of peaks and then only those that are in at least 3 celllines
union.boeva = function(mygrl, n.cl=3, ignore.strand=TRUE){
  
  #https://stackoverflow.com/questions/5516070/multiple-unions
  red.u.grl=Reduce(union, mygrl)
  #red.u.grl$score = countOverlaps(red.u.grl, ignore.strand=ignore.strand)
  
  tmp = lapply( 1:length(mygrl), function(i) countOverlaps(red.u.grl, mygrl[[i]]) )
  ma = do.call(cbind,tmp)
  ma[ma>1]=1
  red.u.grl$score = rowSums(ma)
  
  red.u.grl[ red.u.grl$score >= n.cl ]
  
}


peaks.cl.mycnamp = union.boeva(boeva.grl$celllines.mycnamp)
peaks.cl.nomycnamp = union.boeva(boeva.grl$celllines.nomycnamp)

peaks.ncc = union.boeva(boeva.grl$ncc, n.cl=2) # there are only 2 datasets
peaks.pdx.nomycnamp = union.boeva(boeva.grl$pdx.nomycnamp, n.cl=2)
peaks.pdx.mycnamp = union.boeva(boeva.grl$pdx.mycnamp, n.cl=2)



enhancers = list(
  enhancer.nb.cl.mycnamp=peaks.cl.mycnamp,
  enhancer.nb.cl.nomycnamp=peaks.cl.nomycnamp,
  enhancer.hncc.cl=peaks.ncc,
  enhancer.pdx.mycnamp=peaks.pdx.mycnamp,
  enhancer.pdx.nomycnamp=peaks.pdx.nomycnamp
)



GPL13534pheno=GPL13534.subcl
save(dmrs.dmrcate,enhancers,GPL13534pheno,NB.beta,gene.expr,
     file="/home/kwreczy/projects/motifActivity/data/motifActivity_DMRs_enhancers_GPL13534pheno_NBbeta_gene.expr.RData")


