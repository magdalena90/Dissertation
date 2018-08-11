
setwd('~/MSc/Dissertation')
source('utils.R')

library(sva)
library(WGCNA)

plot_pca = TRUE
######################################################################################################
# LOAD DATA

folder = 'max_var_between' # avg, max_var, max_var_between

load(paste('Data/RDatas', folder, 'LumiBatch_Colantuoni.RData', sep='/'))
exprs_c = exprs(LumiBatch_c)
pData_c = pData(LumiBatch_c)

load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu.RData', sep='/'))
exprs_v = exprs(LumiBatch_v)
pData_v = pData(LumiBatch_v)

# keep only 32 cortex samples
# LumiBatch_v = LumiBatch_v[,colnames(LumiBatch_v) %in% pData_v$GSM[pData_v$Cortex.area!='cerebellum']]
# pData(LumiBatch_v) = pData_v[pData_v$Cortex.area != 'cerebellum',]
# exprs_v = exprs(LumiBatch_v)

######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Original data', folder, '1_by_src')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Original data', folder, '1_by_age',
                           by = 'age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# FILTER PROBES (rows)

# Outliers (>6 mean average deviations from the age-appropriate/disease status linear fit)
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

# Update LumiBatch objects
LumiBatch_c = LumiBatch_c[!rownames(LumiBatch_c) %in% outliers,]
exprs_c = exprs(LumiBatch_c)
LumiBatch_v = LumiBatch_v[!rownames(LumiBatch_v) %in% outliers,]
exprs_v = exprs(LumiBatch_v)

remove(pData, outliers_c, outliers_v, outliers, detect_outliers)
######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering probes', folder, '2_by_src')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering probes', folder, '2_by_age',
                          by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# FILTER SAMPLES

# 1. With low inter-array correlation (Pearson correlation coeffs < 0.85)
detect_low_IAC = function(exprs, pData, threshold=0.85){
  class_name = colnames(pData)[2]
  classes = unique(pData[[class_name]])
  outliers_tot = c()
  for(class in classes){
    class_exprs = exprs[,colnames(exprs) %in% pData$ID[pData[[class_name]] == class]]
    IAC = cor(class_exprs, use='p')
    mean_IAC = apply(IAC, 2, mean)
    outliers = colnames(class_exprs)[mean_IAC<threshold]
    outliers_tot = append(outliers_tot, outliers)
    print(paste0('Class: ', class, ', outliers: ', length(outliers), '/', ncol(class_exprs)))
  }
  outliers_tot = unique(outliers_tot)
  print(paste0('Total outliers: ', length(outliers_tot)))
  return(outliers_tot)
}

pData = data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group)
filter_samples_c = detect_low_IAC(exprs_c, pData, 0.8)

pData = data.frame('ID' = rownames(pData_v), 'disease_status' = pData_v$Disease.status)
filter_samples_v = detect_low_IAC(exprs_v, pData, 0.8)

# 2. Outliers (based on mean inter-array correlation and hierarchical clustering)
detect_sample_outliers = function(exprs, pData, n_sd=2){
  class_name = colnames(pData)[2]
  classes = unique(pData[[class_name]])
  outliers_tot = c()
  for(class in classes){
    class_exprs = exprs[,colnames(exprs) %in% pData$ID[pData[[class_name]] == class]]
    IAC = cor(class_exprs, use='p')
    mean_IAC = apply(IAC, 2, mean)
    number_sd = (mean_IAC - mean(mean_IAC))/sd(mean_IAC)
    #print(plot(number_sd))
    #abline(h=-2)
    outliers = colnames(class_exprs)[number_sd < -n_sd]
    outliers_tot = append(outliers_tot, outliers)
    print(paste0('Class: ', class, ', outliers: ', length(outliers), '/', ncol(class_exprs)))
  }
  outliers_tot = unique(outliers_tot)
  print(paste0('Total outliers: ', length(outliers_tot)))
  return(outliers_tot)
}

pData = data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group)
filter_samples_c = unique(append(filter_samples_c, detect_sample_outliers(exprs_c, pData)))

pData = data.frame('ID' = rownames(pData_v), 'disease_status' = pData_v$Disease.status)
filter_samples_v = unique(append(filter_samples_v, detect_low_IAC(exprs_v, pData)))

# Update LumiBatch objects
LumiBatch_c = LumiBatch_c[,!colnames(LumiBatch_c) %in% filter_samples_c]
exprs_c = exprs(LumiBatch_c)
pData_c = pData_c[rownames(pData_c) %in% colnames(exprs_c),]
pData(LumiBatch_c) = pData_c

LumiBatch_v = LumiBatch_v[,!colnames(LumiBatch_v) %in% filter_samples_v]
exprs_v = exprs(LumiBatch_v)
pData_v = pData_v[rownames(pData_v) %in% colnames(exprs_v),]
pData(LumiBatch_v) = pData_v

remove(detect_low_IAC, pData, filter_samples_c, filter_samples_v, detect_sample_outliers)
######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering samples', folder, '3_by_src')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering samples', folder, '3_by_age',
                           by='age_group')
  remove(pca)
}

######################################################################################################
######################################################################################################
# TRANSFORM SIGNAL TO RATIOS

# Colantuoni
# reference = create_Colantuonis_raw_data_df('SG_Mean', 'SG_Mean', '.')
reference = read.csv('Data/Colantuoni/reference_signal.csv') # Saved df for faster load
reference = reference[!duplicated(reference$HEEBO_ID),]
rownames(reference) = reference$HEEBO_ID
reference$HEEBO_ID = NULL
colnames(reference) = sapply(colnames(reference), function(x) unlist(strsplit(x, '\\.'))[1])

heebo_ilmn_map = data.frame(map_HEEBO_ILMN_ids()[1])
remaining_genes = heebo_ilmn_map[heebo_ilmn_map$ILMN_ID %in% rownames(exprs_c),]

# Collapse probes
reference_genes = remaining_genes$Gene_Symbol[match(rownames(reference), remaining_genes$HEEBO_ID)]
collapse_rows = collapseRows(reference, rowGroup=reference_genes, rowID=rownames(reference), 
                             method = 'Average')
reference = as.data.frame(collapse_rows$datETcollapsed)
reference_ILMN_ids = remaining_genes$ILMN_ID[match(rownames(reference), remaining_genes$Gene_Symbol)]
rownames(reference) = reference_ILMN_ids

ratio_c = (exprs_c+1)/(reference+1)

# Voineagu
row_mean = apply(exprs_v, 1, mean)
ratio_v = (exprs_v+1)/(row_mean+1)

# Update LumiBatch objects
exprs(LumiBatch_c) = as.matrix(ratio_c[,colnames(ratio_c) %in% colnames(exprs_c)])
exprs_c = exprs(LumiBatch_c)

exprs(LumiBatch_v) = as.matrix(ratio_v[,colnames(ratio_v) %in% colnames(exprs_v)])
exprs_v = exprs(LumiBatch_v)

remove(heebo_ilmn_map, ratio_c, row_mean, ratio_v, remaining_genes, reference, collapse_rows)
######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Transforming signal to ratios', 
                           folder, '4_by_src')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Transforming signal to ratios',
                           folder, '4_by_age', by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# Variance stabilisation
LumiBatch_c = lumiT(LumiBatch_c, method = 'log2', ifPlot = TRUE)
LumiBatch_v = lumiT(LumiBatch_v, method = 'log2', ifPlot = TRUE)

# Update expression matrices
exprs_c = exprs(LumiBatch_c)
exprs_v = exprs(LumiBatch_v)

######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Variance stabilisation', folder, 
                           '5_by_src')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Variance stabilisation', folder, 
                           '5_by_age', by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# NORMALISATION

LumiBatch_c = lumiN(LumiBatch_c, method='rsn')
LumiBatch_v = lumiN(LumiBatch_v, method='rsn')

# Update expression matrices
exprs_c = exprs(LumiBatch_c)
exprs_v = exprs(LumiBatch_v)

######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Normalised data', folder, '6_by_src')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Normalised data', folder, '6_by_age',
                           by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# COMBAT BATCH CORRECTION
# http://genomicsclass.github.io/book/pages/svacombat.html
# http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html

exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL

# Transform pDatas to contain the same columns
cols = c('geo_accession', 'age', 'sex', 'pmi', 'age_group', 'Cortex.area', 'Disease.status', 'Source')
pData_c = pData_c[,c('geo_accession','age:ch1','Sex:ch1','postmortem interval (pmi):ch1','age_group')]
pData_c$Cortex.area = 'frontal'
pData_c$Disease.status = 'control'
pData_c$Source = 'Colantuoni'
colnames(pData_c) = cols

pData_v = pData_v[,c('GSM', 'AGE', 'SEX', 'PMI', 'age_group', 'Cortex.area', 'Disease.status')]
pData_v$Source = 'Voineagu'
colnames(pData_v) = cols

pData_c_v = rbind(pData_c, pData_v)
pData_c_v = pData_c_v[rownames(pData_c_v) %in% colnames(exprs_c_v),]
pData_c_v = pData_c_v[match(colnames(exprs_c_v),rownames(pData_c_v)),]

# ComBat
modcombat = model.matrix(~1, data = pData_c_v)
batch = pData_c_v$Source
exprs_c_v = ComBat(exprs_c_v, batch, modcombat)

# Update expression matrices
exprs_c = exprs_c_v[,colnames(exprs_c_v) %in% rownames(pData_c)]
exprs_v = exprs_c_v[,colnames(exprs_c_v) %in% rownames(pData_v)]

######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Batch Correction', folder, '7_by_src')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Batch Correction', folder, '7_by_age', 
                           by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# SAVE PREPROCESSED DATA

exprs(LumiBatch_c) = as.matrix(exprs_c)
exprs(LumiBatch_v) = as.matrix(exprs_v)

save(LumiBatch_c, file=paste0('Data/RDatas/',folder,'/LumiBatch_Colantuoni_preprocessed.RData'))
save(LumiBatch_v, file=paste0('Data/RDatas/',folder,'/LumiBatch_Voineagu_preprocessed.RData'))


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# PROCESS TEMPORAL AND FRONTAL SAMPLES SEPARATELY FROM CEREBELLUM

# Filter probes and samples as in preprocessing of full data:
folder = 'max_var_between'
load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed.RData', sep='/'))
exprs_v = exprs(LumiBatch_v)
pData_v = pData(LumiBatch_v)
load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu.RData', sep='/'))
LumiBatch_v = LumiBatch_v[rownames(exprs(LumiBatch_v)) %in% rownames(exprs_v),
                          colnames(exprs(LumiBatch_v)) %in% rownames(pData_v)]
pData(LumiBatch_v) = pData(LumiBatch_v)[rownames(pData(LumiBatch_v)) %in% rownames(pData_v),]
exprs_v = exprs(LumiBatch_v)
pData_v = pData(LumiBatch_v)


LumiBatch = LumiBatch_v[,colnames(exprs_v) %in% rownames(pData_v)[pData_v$Cortex.area!='cerebellum']]
exprs = exprs(LumiBatch)
pData = pData(LumiBatch)

# Transform Signal to Ratio
row_mean = apply(exprs, 1, mean)
ratio = (exprs+1)/(row_mean+1)
exprs(LumiBatch) = as.matrix(ratio[,colnames(ratio) %in% colnames(exprs)])

# Variance Stabilisation
LumiBatch = lumiT(LumiBatch, method = 'log2', ifPlot = TRUE)

# Normalisation
LumiBatch = lumiN(LumiBatch, method='rsn')

LumiBatch_ft = LumiBatch
save(LumiBatch_ft, file=paste0('Data/RDatas/',folder,'/LumiBatch_Voineagu_preprocessed_ft.RData'))


remove(exprs_exprs_c, exprs_v, LumiBatch, LumiBatch_c, LumiBatch_v, LumiBatch_ft, pData, pData_c, 
       pData_v, phenoData, row_mean)
