
setwd('~/MSc/Dissertation')
source('utils.R')

library(lumi)
library(GEOquery)

######################################################################################################
# COLANTUONI

# Load data for exprs and se.exprs dataframes
assayData = create_Colantuonis_raw_data_df(cols = c('SR_Mean','SR_S.Dev'), 
                                           colnames = c('exprs','se'), sep = '.')
assayData = assayData[!duplicated(assayData$HEEBO_ID),]
rownames(assayData) = assayData$HEEBO_ID
assayData$HEEBO_ID = NULL

# Filter probes to one-to-one relation to ilmn ids
heebo_ilmn_map = map_HEEBO_ILMN_ids()
heebo_ilmn_map = data.frame(heebo_ilmn_map[2])    # without control samples
unique_gene_mapping = heebo_ilmn_map[!duplicated(heebo_ilmn_map$Gene_Symbol) & 
                                     !duplicated(heebo_ilmn_map$Gene_Symbol, fromLast = TRUE),]
assayData = assayData[rownames(assayData) %in% unique_gene_mapping$HEEBO_ID,]

# Change IDs to ILMN notation
heebo_to_ilmn = unique_gene_mapping$ILMN_ID
names(heebo_to_ilmn) = unique_gene_mapping$HEEBO_ID
rownames(assayData) = heebo_to_ilmn[rownames(assayData)]

# Create exprs and se.exprs dataframes
exprs = assayData[,grepl('exprs', colnames(assayData))]
se.exprs = assayData[,grepl('se', colnames(assayData))]
colnames(exprs) = sapply(colnames(exprs), function(x) unlist(strsplit(x, '\\.'))[1])
colnames(se.exprs) = sapply(colnames(se.exprs), function(x) unlist(strsplit(x, '\\.'))[1])

# Phenotype data
gse = getGEO(filename = 'Data/Colantuoni/GSE30272_series_matrix.txt.gz', getGPL = FALSE)
phenoData = pData(gse)
age_group = cut(as.numeric(phenoData$`age:ch1`), c(-0.5, 0, 0.6, 10, 20, 30, 40, 50, 60, 70, 80),
                labels=c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s'))
phenoData$age_group = as.character(age_group)

# LumiBatch object
LumiBatch_c = new('LumiBatch', exprs = as.matrix(exprs), se.exprs = as.matrix(se.exprs))
pData(LumiBatch_c) = phenoData


save(LumiBatch_c, file='Data/LumiBatch_Colantuoni.RData')


remove(assayData, exprs, se.exprs, gse, phenoData, heebo_ilmn_map, heebo_to_ilmn, age_group)
######################################################################################################
# VOINEAGU

# AssayData data
LumiBatch_v = lumiR('Data/Voineagu/GSE28521_non-normalized_data.txt.gz')

# Phenotype data
phenoData = create_Voineagus_pData_df()
phenoData = phenoData[match(colnames(exprs(LumiBatch_v)), phenoData$chip_array),]
rownames(phenoData) = phenoData$GSM

# Change subject IDs from chip_array to GSM format
chipArray_to_gsm = phenoData$GSM
names(chipArray_to_gsm) = phenoData$chip_array
colnames(LumiBatch_v) = chipArray_to_gsm[colnames(LumiBatch_v)]

# Filter probes to one-to-one relation to ilmn ids
LumiBatch_v = LumiBatch_v[rownames(LumiBatch_v) %in% unique_gene_mapping$ILMN_ID,]

# LumiBatch object
pData(LumiBatch_v) = phenoData


save(LumiBatch_v, file='Data/LumiBatch_Voineagu.RData')


remove(phenoData, chipArray_to_gsm, unique_gene_mapping)
