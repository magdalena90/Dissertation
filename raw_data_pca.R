
setwd('~/MSc/Dissertation')
source('utils.R')

######################################################################################################
# Load data

load('Data/LumiBatch_Colantuoni.RData')
exprs_c = exprs(LumiBatch_c)
pData_c = pData(LumiBatch_c)

load('Data/LumiBatch_Voineagu.RData')
exprs_v = exprs(LumiBatch_v)
pData_v = pData(LumiBatch_v)


remove(LumiBatch_c, LumiBatch_v)

######################################################################################################
# PCA

# Colantuoni labelled by age
labels = data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group)
pca = perform_pca(exprs_c, labels, 'Colantuoni raw data by age')

# Voineagu labelled by brain region
labels = data.frame('ID' = rownames(pData_v), 'cortex_area' = pData_v$Cortex.area)
pca = perform_pca(exprs_v, labels, 'Voineagu raw data by cortex area')

# Voineagu filtered by frontal cortex and labelled by disease status
exprs_fc_v = exprs_v[,colnames(exprs_v) %in% pData_v$GSM[pData_v$Cortex.area=='frontal']]
labels = data.frame('ID' = rownames(pData_v), 'disease_status' = pData_v$Disease.status)
pca = perform_pca(exprs_fc_v, labels, 'Voineagu frontal cortex data')

# Colantuoni vs Voineagu frontal cortex
exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
rownames(exprs_c_v) = exprs_c_v$Row.names
exprs_c_v$Row.names = NULL
labels_c = data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni')
labels_v = data.frame('ID' = rownames(pData_v), 'Source' = pData_v$Disease.status)
labels = rbind(labels_c, labels_v)
pca = perform_pca(exprs_c_v, labels, 'Colantuoni vs Voineagu data')

remove(labels, labels_c, labels_v, pca)
