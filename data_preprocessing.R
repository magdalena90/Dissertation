
setwd('~/MSc/Dissertation')
source('utils.R')

library(lumi)
library(limma)
library(QCGWAS)
library(sva)

######################################################################################################
# LOAD DATA

load('Data/LumiBatch_Colantuoni.RData')
exprs_c = exprs(LumiBatch_c)
pData_c = pData(LumiBatch_c)

load('Data/LumiBatch_Voineagu.RData')
exprs_v = exprs(LumiBatch_v)
pData_v = pData(LumiBatch_v)

# keep only 32 frontal cortex samples
LumiBatch_v = LumiBatch_v[,colnames(LumiBatch_v) %in% pData_v$GSM[pData_v$Cortex.area!='cerebellum']]
pData(LumiBatch_v) = pData_v[pData_v$Cortex.area != 'cerebellum',]
exprs_v = exprs(LumiBatch_v)

######################################################################################################
# PCA Check:
exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL
labels = rbind(data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni'),
               data.frame('ID' = rownames(pData_v), 'Source' = 'Voineagu'))
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu original data')

remove(exprs_c_v, labels, pca)
######################################################################################################
######################################################################################################
# FILTER PROBES

remove_probes = c()
######################################################################################################
# 1. VOINEAGU FILTER: remove probes that are not robustly expressed
#    (detection p value >0.05 for more than half of the samples)
not_robustly_expressed_probes = function(det, threshold=0.05){
  det_corrected = 1-det
  perc_above_threshold_by_row = apply(det_corrected, 1, function(x) sum(x>threshold)/ncol(det))
  print(hist(perc_above_threshold_by_row, breaks=100))
  filter_rows = perc_above_threshold_by_row > 0.9
  filter_rows = rownames(det[filter_rows,])
  return(filter_rows)
}

detection_v = detection(LumiBatch_v)
filter_v = not_robustly_expressed_probes(detection_v)

remove_probes = append(remove_probes, filter_v)

remove(detection_v, filter_v, not_robustly_expressed_probes)
######################################################################################################
# 2. COLANTUONI FILTERS:
# 2.1 Containing polymorphisms with minor allele frequency >0.01 according to HapMap in either 
#     YRI or CEU populations


# 2.2 Outliers (>6 mean average deviations from the age-appropriate/disease status linear fit)
#     Probes (rows) are the observations and subjects (columns) the variables
detect_outliers = function(exprs, pData){
  class_name = colnames(pData)[2]
  classes = unique(pData[[class_name]])
  outliers_tot = c()
  for(class in classes){
    class_exprs = exprs[,colnames(exprs) %in% pData$ID[pData[[class_name]] == class]]
    l_fit = lm(data.frame(class_exprs))
    outliers = names(rstandard(l_fit))[abs(rstandard(l_fit))>6]
    outliers_tot = append(outliers_tot, outliers)
    print(paste0('Class: ',class, ', outliers: ',length(outliers)))
  }
  outliers_tot = unique(outliers_tot)
  print(paste0('Total outliers: ', length(outliers_tot)))
  return(outliers_tot)
}

pData = data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group)
outliers_c = detect_outliers(exprs_c, pData)

pData = data.frame('ID' = as.character(rownames(pData_v)), 'disease_status' = pData_v$Disease.status)
outliers_v = detect_outliers(exprs_v, pData)

outliers = unique(append(outliers_c, outliers_v))

remove_probes = unique(append(remove_probes, outliers))

LumiBatch_c = LumiBatch_c[!rownames(LumiBatch_c) %in% remove_probes,]
exprs_c = exprs(LumiBatch_c)
LumiBatch_v = LumiBatch_v[!rownames(LumiBatch_v) %in% remove_probes,]
exprs_v = exprs(LumiBatch_v)

remove(pData, outliers_c, outliers_v, outliers, detect_outliers)
######################################################################################################
# PCA Check:
exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL
labels = rbind(data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni'),
               data.frame('ID' = rownames(pData_v), 'Source' = 'Voineagu'))
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu after probe filtering')

remove(exprs_c_v, labels, pca, filter_probes)
######################################################################################################
######################################################################################################
# FILTER SAMPLES

filter_samples_c = c()
filter_samples_v = c()
######################################################################################################
# 1. With low inter-array correlation (Pearson correlation coeffs < 0.85)
IAC_c = cor(exprs_c, use = 'p')
meanIAC_c = apply(IAC_c, 2, mean)

filter_samples_c = append(filter_samples_c, colnames(exprs_c)[meanIAC_c<0.81])

IAC_v = cor(exprs_v, use = 'p')
meanIAC_v = apply(IAC_v, 2, mean)

filter_samples_v = append(filter_samples_v, colnames(exprs_v)[meanIAC_v<0.81])

# 2. Outliers (based on mean inter-array correlation and hierarchical clustering)
numbersd_c = (meanIAC_c-mean(meanIAC_c))/sd(meanIAC_c)
plot(numbersd_c)
abline(h=-2)

filter_samples_c = unique(append(filter_samples_c, rownames(pData_c)[numbersd_c < -2]))

numbersd_v = (meanIAC_v-mean(meanIAC_v))/sd(meanIAC_v)
plot(numbersd_v)
abline(h=-2)

filter_samples_v = unique(append(filter_samples_v, rownames(pData_v)[numbersd_v < -2]))


LumiBatch_c = LumiBatch_c[,!colnames(LumiBatch_c) %in% filter_samples_c]
exprs_c = exprs(LumiBatch_c)
pData_c = pData_c[rownames(pData_c) %in% colnames(exprs_c),]
pData(LumiBatch_c) = pData_c

LumiBatch_v = LumiBatch_v[,!colnames(LumiBatch_v) %in% filter_samples_v]
exprs_v = exprs(LumiBatch_v)
pData_v = pData_v[rownames(pData_v) %in% colnames(exprs_v),]
pData(LumiBatch_v) = pData_v

remove(IAC_c, meanIAC_c, numbersd_c, filter_c, IAC_v, meanIAC_v, numbersd_v, filter_v)
######################################################################################################
# PCA Check:
exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL
labels = rbind(data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni'),
               data.frame('ID' = rownames(pData_v), 'Source' = 'Voineagu'))
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu after sample filtering')

remove(exprs_c_v, labels, pca)
######################################################################################################
######################################################################################################
# Background correction

lumiB(LumiBatch_c, method='bgAdjust') # The data has already been background adjusted!
lumiB(LumiBatch_v, method='bgAdjust') # There is no control probe information in the LumiBatch object!

######################################################################################################
######################################################################################################
# TRANSFORM SIGNAL TO RATIOS

# Colantuoni
reference = create_Colantuonis_raw_data_df('SG_Mean', 'SG_Mean', '.')
reference = reference[!duplicated(reference$HEEBO_ID),]
rownames(reference) = reference$HEEBO_ID
reference$HEEBO_ID = NULL
colnames(reference) = sapply(colnames(reference), function(x) unlist(strsplit(x, '\\.'))[1])

heebo_ilmn_map = data.frame(map_HEEBO_ILMN_ids()[2])
remaining_genes = heebo_ilmn_map[heebo_ilmn_map$ILMN_ID %in% rownames(exprs_c),]
heebo_to_ilmn = remaining_genes$ILMN_ID
names(heebo_to_ilmn) = remaining_genes$HEEBO_ID

reference = reference[rownames(reference) %in% names(heebo_to_ilmn),]
rownames(reference) = heebo_to_ilmn[rownames(reference)]

ratio_c = (exprs_c+1)/(reference+1)


exprs(LumiBatch_c) = as.matrix(ratio_c[,colnames(ratio_c) %in% colnames(exprs_c)])
exprs_c = exprs(LumiBatch_c)

# Voineagu
row_mean = apply(exprs_v, 1, mean)
ratio_v = (exprs_v+1)/(row_mean+1)

exprs(LumiBatch_v) = as.matrix(ratio_v[,colnames(ratio_v) %in% colnames(exprs_v)])
exprs_v = exprs(LumiBatch_v)

remove(reference, heebo_ilmn_map, ratio_c, row_mean, ratio_v)
######################################################################################################
# PCA Check:
exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL
labels = rbind(data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni'),
               data.frame('ID' = rownames(pData_v), 'Source' = 'Voineagu'))
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu ratios')

labels = rbind(data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group),
               data.frame('ID' = rownames(pData_v), 'age_group' = pData_v$age_group))
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu ratios by age')

labels = data.frame('ID' = rownames(pData_v), 'disease_status' = pData_v$Disease.status)
pca = perform_pca(exprs_v, labels, 'Voineagu ratio data by disease status')

remove(exprs_c_v, labels, pca)
######################################################################################################
######################################################################################################
# Variance stabilisation
LumiBatch_c = lumiT(LumiBatch_c, method = 'log2', ifPlot = TRUE)
LumiBatch_v = lumiT(LumiBatch_v, method = 'log2', ifPlot = TRUE)

exprs_c = exprs(LumiBatch_c)
exprs_v = exprs(LumiBatch_v)

######################################################################################################
# PCA Check:
exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL
labels = rbind(data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni'),
               data.frame('ID' = rownames(pData_v), 'Source' = 'Voineagu'))
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu log2 data')

remove(exprs_c_v, labels, pca)
######################################################################################################
######################################################################################################
# NORMALISATION

LumiBatch_c = lumiN(LumiBatch_c, method='rsn')
exprs_c = exprs(LumiBatch_c)

LumiBatch_v = lumiN(LumiBatch_v, method='rsn')
exprs_v = exprs(LumiBatch_v)

######################################################################################################
# PCA Check:
exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL
labels = rbind(data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni'),
               data.frame('ID' = rownames(pData_v), 'Source' = 'Voineagu'))
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu normalised data')

remove(exprs_c_v, labels, pca)
######################################################################################################
######################################################################################################
# COMBAT BATCH CORRECTION
# http://genomicsclass.github.io/book/pages/svacombat.html
# http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html

exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL

pData_c = pData_c[,c('geo_accession','age:ch1','Sex:ch1','postmortem interval (pmi):ch1','age_group')]
pData_c$Cortex.area = 'frontal'
pData_c$Disease.status = 'control'
pData_c$Source = 'Colantuoni'
colnames(pData_c) = c('geo_accession', 'age', 'sex', 'pmi', 'age_group', 'Cortex.area',
                      'Disease.status', 'Source')
pData_v = pData_v[,c('GSM', 'AGE', 'SEX', 'PMI', 'age_group', 'Cortex.area', 'Disease.status')]
pData_v$Source = 'Voineagu'
colnames(pData_v) = c('geo_accession', 'age', 'sex', 'pmi', 'age_group', 'Cortex.area',
                      'Disease.status', 'Source')
pData_c_v = rbind(pData_c, pData_v)
pData_c_v = pData_c_v[rownames(pData_c_v) %in% colnames(exprs_c_v),]
pData_c_v = pData_c_v[match(colnames(exprs_c_v),rownames(pData_c_v)),]

modcombat = model.matrix(~1, data = pData_c_v)
batch = pData_c_v$Source

exprs_c_v = ComBat(exprs_c_v, batch, modcombat)

######################################################################################################
# PCA Check:
labels = data.frame('ID' = rownames(pData_c_v), 'disease_status' = pData_c_v$Source)
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu ComBat data')

labels = data.frame('ID' = rownames(pData_c_v), 'age_group' = pData_c_v$age_group)
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu ComBat data by age')

remove(exprs_c_v, labels, pca)



######################################################################################################
######################################################################################################
######################################################################################################
# SVA: https://www.bioconductor.org/packages/devel/bioc/vignettes/sva/inst/doc/sva.pdf

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3307112/

# Use combat for each of the confounding variables one at a time?

# SVA should be used with caution for creating networks or to illustrate patterns (SVA):
# https://support.bioconductor.org/p/42615/

library(sva)
library(limma)

pheno = pData(gse)
edata = exprs(gse)

# Change continuous age to age groups
age_groups = cut(as.numeric(pheno$`age:ch1`), c(-0.5, 0, 0.6, 10, 20, 30, 40, 50, 60, 70, 80),
                 labels=c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s'))
pheno$`age:ch1` = age_groups

# SVA
full_model = model.matrix(~as.factor(`age:ch1`), data=pheno)
null_model = model.matrix(~1, data=pheno)
n_sv = num.sv(edata, full_model)   # 32
svobj = sva(edata, full_model, null_model, n.sv=n_sv)

# Adjusting for surrogate vars using f.pvalue
pValues = f.pvalue(edata, full_model, null_model)
qValues = p.adjust(pValues, method = 'BH')

full_model_Sv = cbind(full_model, svobj$sv)
null_model_Sv = cbind(null_model, svobj$sv)

design(full_model_Sv) = ~ V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24+V25+V26+V27+V28+
  V29+V30+V31+V32+V33+V34+V35+V36+V37+V38+V39+V40+V41+V42+
  "(Intercept)"+"as.factor(`age:ch1`)Infant"+"as.factor(`age:ch1`)Child"+ 
  "as.factor(`age:ch1`)10-20"+"as.factor(`age:ch1`)20s"+"as.factor(`age:ch1`)30s"+  
  "as.factor(`age:ch1`)40s"+"as.factor(`age:ch1`)50s"+"as.factor(`age:ch1`)60s"+   
  "as.factor(`age:ch1`)70s" 


pValuesSv = f.pvalue(edata, full_model_Sv, null_model_Sv)
qValuesSv = p.adjust(pValuesSv, method ='BH')

# Adjusting for surrogate vars using limma package
fit = lmFit(edata, full_model_Sv)

contrast.matrix <- cbind('C1' =c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, rep(0, svobj$n.sv)), 
                         'C2' =c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, rep(0, svobj$n.sv)),
                         'C3' =c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, rep(0, svobj$n.sv)),
                         'C4' =c(0, 0, 0, -1, 1, 0, 0, 0, 0, 0, rep(0, svobj$n.sv)),
                         'C5' =c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, rep(0, svobj$n.sv)),
                         'C6' =c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, rep(0, svobj$n.sv)),
                         'C7' =c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0, rep(0, svobj$n.sv)),
                         'C8' =c(0, 0, 0, 0, 0, 0, 0, -1, 1, 0, rep(0, svobj$n.sv)),
                         'C9' =c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1, rep(0, svobj$n.sv)),
                         'C10'=c(1, 0, 0, 0, 0, 0, 0, 0, 0, -1, rep(0, svobj$n.sv)))

fitContrasts = contrasts.fit(fit, contrast.matrix)
eb = eBayes(fitContrasts)
topTableF(eb, adjust='BH')

#assay_data2 = read.delim('Data/Colantuoni/GSE30272_ExprsMtxCleanedN269_31SVN.txt')





