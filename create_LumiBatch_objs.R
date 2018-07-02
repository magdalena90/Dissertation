
setwd('~/MSc/Dissertation')
source('utils.R')

library(lumi)
library(GEOquery)

######################################################################################################
# COLANTUONI

# Load data for exprs and se.exprs dataframes
exprs_c = create_Colantuonis_raw_data_df(cols = c('SR_Mean','SR_S.Dev'), 
                                         colnames = c('exprs','se'), sep = '.')
exprs_c = exprs_c[!duplicated(exprs_c$HEEBO_ID),]
rownames(exprs_c) = exprs_c$HEEBO_ID
exprs_c$HEEBO_ID = NULL

# Filter probes to one-to-one relation to ilmn ids
heebo_ilmn_map = map_HEEBO_ILMN_ids()
heebo_ilmn_map = data.frame(heebo_ilmn_map[2])    # without control samples
unique_gene_mapping = heebo_ilmn_map[!duplicated(heebo_ilmn_map$Gene_Symbol) & 
                                     !duplicated(heebo_ilmn_map$Gene_Symbol, fromLast=T),]
exprs_c = exprs_c[rownames(exprs_c) %in% unique_gene_mapping$HEEBO_ID,]

# Change IDs to ILMN notation
heebo_to_ilmn = unique_gene_mapping$ILMN_ID
names(heebo_to_ilmn) = unique_gene_mapping$HEEBO_ID
rownames(exprs_c) = heebo_to_ilmn[rownames(exprs_c)]

# Create exprs and se.exprs dataframes
exprs = exprs_c[,grepl('exprs', colnames(exprs_c))]
se.exprs = exprs_c[,grepl('se', colnames(exprs_c))]
colnames(exprs) = sapply(colnames(exprs), function(x) unlist(strsplit(x, '\\.'))[1])
colnames(se.exprs) = sapply(colnames(se.exprs), function(x) unlist(strsplit(x, '\\.'))[1])

# Phenotype data
gse_c = getGEO(filename='Data/Colantuoni/GSE30272_series_matrix.txt.gz', getGPL = FALSE)
samples_c = pData(gse_c)

# LumiBatch object
LumiBatch_c = new('LumiBatch', exprs=as.matrix(exprs), se.exprs=as.matrix(se.exprs), pData=samples_c)


save(LumiBatch_c, file='Data/LumiBatch_Colantuoni.RData')


remove(exprs_c, exprs, se.exprs, gse_c, samples_c, heebo_ilmn_map, heebo_to_ilmn)
######################################################################################################
# VOINEAGU

# Expression data
lumiBatch_v = lumiR('Data/Voineagu/GSE28521_non-normalized_data.txt.gz')

# Phenotype data
samples_v = create_Voineagus_pData_df()
samples_v = samples_v[match(colnames(exprs(lumiBatch_v)), samples_v$chip_array),]
rownames(samples_v) = samples_v$chip_array

# LumiBatch object
pData(lumiBatch_v) = samples_v


save(LumiBatch_v, file='Data/LumiBatch_Voineagu.RData')


remove(samples_v, unique_gene_mapping)
