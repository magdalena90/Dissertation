
setwd('~/MSc/Dissertation')
source('utils.R')

library(lumi)
library(sva)

plot_pca = TRUE
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
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Original data')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Original data', by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# FILTER PROBES (rows)

# 1 Drop low intensity values and filter probes with <0.5 fetal or postnatal data
exprs_c[exprs_c<39.4] = NA
detect_high_NA_perc = function(exprs, pData, perc=0.5){
  class_name = colnames(pData)[2]
  classes = unique(pData[[class_name]])
  drop_tot = c()
  for(class in classes){
    class_exprs = exprs[,colnames(exprs) %in% pData$ID[pData[[class_name]] == class]]
    perc_NA = apply(class_exprs, 1, function(x) sum(is.na(x))/length(x))
    drop = names(perc_NA)[perc_NA > perc]
    drop_tot = append(drop_tot, drop)
    print(paste0('Class: ', class, ', drop: ', length(drop)))
  }
  drop_tot = unique(drop_tot)
  print(paste0('Total drops: ', length(drop_tot)))
  return(drop_tot)
}

pData = data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group)
pData$age_group = as.character(pData$age_group)
pData$age_group[pData$age_group != 'Fetal'] = 'Postnatal'
high_NA_perc = detect_high_NA_perc(exprs_c, pData)

exprs(LumiBatch_c) = exprs_c
LumiBatch_c = LumiBatch_c[!rownames(LumiBatch_c) %in% high_NA_perc,]
exprs_c = exprs(LumiBatch_c)
LumiBatch_v = LumiBatch_v[!rownames(LumiBatch_v) %in% high_NA_perc,]
exprs_v = exprs(LumiBatch_v)

# 2 Outliers (>6 mean average deviations from the age-appropriate/disease status linear fit)
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

remove(pData, outliers_c, outliers_v, high_NA_perc, outliers, detect_high_NA_perc, detect_outliers)
######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering probes')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering probes', by='age_group')
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
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering samples')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Filtering samples', by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# # Background correction
# 
# lumiB(LumiBatch_c, method='bgAdjust') # The data has already been background adjusted!
# lumiB(LumiBatch_v, method='bgAdjust') # There is no control probe information in the LumiBatch object!
# 
# # Dataset comparison
# exprs_v_control=exprs_c[,colnames(exprs_v) %in% rownames(pData_v[pData_v$Disease.status=='control',])]
# mean_exprs_v = data.frame('ID'=rownames(exprs_v), 'mean_v'=apply(exprs_v, 1, mean))
# 
# exprs_c_non_fetal = exprs_c[,colnames(exprs_c) %in% rownames(pData_c[pData_c$age_group != 'Fetal',])]
# mean_exprs_c = data.frame('ID'=rownames(exprs_c_non_fetal), 'mean_c'=apply(exprs_c_non_fetal, 1, mean))
# mean_exprs_c = mean_exprs_c[mean_exprs_c$mean_c > min(mean_exprs_v$mean_v)*2,]
# 
# mean_exprs = merge(mean_exprs_c, mean_exprs_v, by='ID')
# rownames(mean_exprs) = mean_exprs$ID
# mean_exprs$ID = NULL
# 
# library(ggplot2)
# ggplot(log(mean_exprs), aes(x=mean_c, y=mean_v)) + geom_point(alpha=0.5) +
#   geom_smooth(method='lm', formula=y~x) + theme_minimal()

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

heebo_ilmn_map = data.frame(map_HEEBO_ILMN_ids()[2])
remaining_genes = heebo_ilmn_map[heebo_ilmn_map$ILMN_ID %in% rownames(exprs_c),]
heebo_to_ilmn = remaining_genes$ILMN_ID
names(heebo_to_ilmn) = remaining_genes$HEEBO_ID

reference = reference[rownames(reference) %in% names(heebo_to_ilmn),]
rownames(reference) = heebo_to_ilmn[rownames(reference)]

ratio_c = (exprs_c+1)/(reference+1)

# Voineagu
row_mean = apply(exprs_v, 1, mean)
ratio_v = (exprs_v+1)/(row_mean+1)

# Update LumiBatch objects
exprs(LumiBatch_c) = as.matrix(ratio_c[,colnames(ratio_c) %in% colnames(exprs_c)])
exprs_c = exprs(LumiBatch_c)

exprs(LumiBatch_v) = as.matrix(ratio_v[,colnames(ratio_v) %in% colnames(exprs_v)])
exprs_v = exprs(LumiBatch_v)

remove(heebo_ilmn_map, ratio_c, row_mean, ratio_v, heebo_to_ilmn, remaining_genes, reference)
######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Transforming signal to ratios')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Transforming signal to ratios',
                           by='age_group')
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
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Variance stabilisation')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Variance stabilisation',
                           by='age_group')
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
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Normalised data')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Normalised data',
                           by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# # Impute missing data
# 
# # Data seems to distribute ~Normal, so impute data with rnorm(row mean, row sd)
# impute_data = function(exprs, pData){
#   class_name = colnames(pData)[2]
#   classes = unique(pData[[class_name]])
#   is_na = is.na(exprs)
#   for(class in classes){
#     class_exprs = exprs[, colnames(exprs) %in% pData$ID[pData[[class_name]] == class]]
#     class_is_na = is_na[, colnames(is_na) %in% pData$ID[pData[[class_name]] == class]]
#     na_counts = apply(class_is_na, 1, sum)
#     names(na_counts) = rownames(class_is_na)
#     rows_with_nas = names(na_counts)[na_counts>0]
#     for(row in rows_with_nas){
#       row_data = class_exprs[row,]
#       if(sum(!is.na(row_data))==0) {
#         if(class == 'Fetal'){
#           print('No fetal samples')
#         } else {
#           print('No samples')
#           print(class)
#           row_data = exprs[row, colnames(exprs) %in% pData$ID[pData$age_group!='Fetal']]
#           rand_samps = rnorm(sum(is.na(row_data)), mean(row_data, na.rm=T), sd(row_data, na.rm=T)/2)
#           rand_samps = rep(min(row_data, na.rm=T), sum(is.na(row_data)))
#           #rand_samps = rep(NA, sum(is.na(row_data)))
#           cols = names(row_data)[is.na(row_data)]
#           exprs[row, cols] = rand_samps
#         }
#       } else {
#         if(sum(!is.na(row_data))==1){
#           if(class == 'Fetal'){
#             print('1 fetal sample')
#           } else {
#             print('1 sample')
#             print(class)
#             rand_samps = rep(min(row_data, na.rm=T), sum(is.na(row_data)))
#             cols = names(row_data)[is.na(row_data)]
#             exprs[row, cols] = rand_samps
#           }
#         } else {
#           rand_samps = rnorm(sum(is.na(row_data)), mean(row_data, na.rm=T), sd(row_data, na.rm=T))
#           rand_samps = rep(min(row_data, na.rm=T), sum(is.na(row_data)))
#           #rand_samps = rep(NA, sum(is.na(row_data)))
#           cols = names(row_data)[is.na(row_data)]
#           exprs[row, cols] = rand_samps
#         }
#       }
#     }
#   }
#   return(exprs)
# }
# 
# na_counts = apply(is.na(exprs_c), 1, sum)
# names(na_counts) = rownames(exprs_c)
# with_nas = names(na_counts)[na_counts>10]
# prueba = exprs_c[with_nas[1:10],]
# 
# pData = data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group)
# prueba2 = impute_data(exprs_c, pData)
# 
# # Prueba con matriz chica
# a = data.frame('a'=c(1,2,3,4,5),'b'=c(5,6,7,8,9))
# a$a1 = a$a; a$a2 = a$a; a$a3 = a$a; a$a4 = a$a
# a$b1 = a$b; a$b2 = a$b; a$b3 = a$b; a$b4 = a$b
# a[3, sample(ncol(a),5)] = NA
# pData_a = data.frame('ID'=seq(1,5), 'class'=c('a','b','a','a','a','a','b','b','b','b'))
# a_impute = impute_data(a, pData_a)
# all.equal(a_impute, a)
# 
# exprs_c2 = exprs_c
# k = which(is.na(exprs_c2), arr.ind=TRUE)
# exprs_c2[k] <- rowMeans(exprs_c2, na.rm=TRUE)[k[,1]]
# 
# exprs(LumiBatch_c) = exprs_c

######################################################################################################
# PCA Check:
if(plot_pca){
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Imputing missing data')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Imputing missing data',
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
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Batch Correction')
  pca = prepare_pca_c_vs_v(exprs_c, exprs_v, pData_c, pData_v, 'Batch Correction', by='age_group')
  remove(pca)
}
######################################################################################################
######################################################################################################
# SAVE PREPROCESSED DATA

exprs(LumiBatch_c) = as.matrix(exprs_c)
exprs(LumiBatch_v) = as.matrix(exprs_v)

save(LumiBatch_c, file='Data/LumiBatch_Colantuoni_preprocessed.RData')
save(LumiBatch_v, file='Data/LumiBatch_Voineagu_preprocessed.RData')
