
setwd('~/MSc/Dissertation')
source('utils.R')

library(lumi)
library(GEOquery)
library(WGCNA)

######################################################################################################
######################################################################################################
# COLANTUONI

# Load data for exprs and se.exprs dataframes
assayData = create_Colantuonis_raw_data_df(cols = c('SR_Mean','SR_S.Dev'), 
                                           colnames = c('exprs','se'), sep = '.')
assayData = assayData[!duplicated(assayData$HEEBO_ID),]
rownames(assayData) = assayData$HEEBO_ID
assayData$HEEBO_ID = NULL

# Get ID-gene relation
gpl = getGEO(filename='Data/Colantuoni/GPL4611_family_short.soft', getGPL = FALSE)
heebo_gene_map = data.frame('ID' = Table(gpl)[,'ID'], 'Gene_Symbol' = Table(gpl)[,'Gene_Symbol'])
heebo_gene_map = heebo_gene_map[!heebo_gene_map$Gene_Symbol %in% c('','EMPTY','##noname##'),]
heebo_gene_map$ID = as.character(heebo_gene_map$ID)
assayData = assayData[rownames(assayData) %in% heebo_gene_map$ID,]

# Average probe values corresponding to a same gene
assayData_genes = heebo_gene_map$Gene_Symbol[match(rownames(assayData), heebo_gene_map$ID)]
collapse_rows = collapseRows(assayData, rowGroup=assayData_genes, rowID=rownames(assayData), 
                             method = 'Average')
assayData = as.data.frame(collapse_rows$datETcollapsed)

# Filter probes corresponding to genes not present in Voineagu's experiment
heebo_ilmn_map = data.frame(map_HEEBO_ILMN_ids()[1])
assayData = assayData[rownames(assayData) %in% heebo_ilmn_map$Gene_Symbol,]

# Change IDs to ILMN notation (keep only 1st ILMN ID for each gene)
gene_ilmn_map = unique(heebo_ilmn_map[,c('Gene_Symbol','ILMN_ID')])
gene_ilmn_map = gene_ilmn_map[!duplicated(gene_ilmn_map$Gene_Symbol),]
assayData_ILMN_ids = gene_ilmn_map$ILMN_ID[match(rownames(assayData), gene_ilmn_map$Gene_Symbol)]
rownames(assayData) = assayData_ILMN_ids

# # Filter probes to one-to-one relation to ilmn ids
# heebo_ilmn_map = map_HEEBO_ILMN_ids()
# heebo_ilmn_map = data.frame(heebo_ilmn_map[2])    # without control samples
# unique_gene_mapping = heebo_ilmn_map[!duplicated(heebo_ilmn_map$Gene_Symbol) & 
#                                      !duplicated(heebo_ilmn_map$Gene_Symbol, fromLast = TRUE),]
# assayData = assayData[rownames(assayData) %in% unique_gene_mapping$HEEBO_ID,]
# 
# # Change IDs to ILMN notation
# heebo_to_ilmn = unique_gene_mapping$ILMN_ID
# names(heebo_to_ilmn) = unique_gene_mapping$HEEBO_ID
# rownames(assayData) = heebo_to_ilmn[rownames(assayData)]

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

save(LumiBatch_c, file='Data/LumiBatch_Colantuoni_big.RData')


remove(collapse_rows, gpl, heebo_gene_map, assayData, exprs, se.exprs, gse, assayData_genes, 
       assayData_ILMN_ids, phenoData, age_group)
######################################################################################################
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

# Filter probes to match to Colantuoni's genes and average probes corresponding to the same gene
assayData = data.frame(exprs(LumiBatch_v))
assayData = assayData[rownames(assayData) %in% heebo_ilmn_map$ILMN_ID,]
assayData_genes = heebo_ilmn_map$Gene_Symbol[match(rownames(assayData), heebo_ilmn_map$ILMN_ID)]
collapse_rows = collapseRows(assayData, rowGroup=assayData_genes, rowID=rownames(assayData), 
                             method = 'Average')
assayData = as.data.frame(collapse_rows$datETcollapsed)
assayData_ILMN_ids = gene_ilmn_map$ILMN_ID[match(rownames(assayData), gene_ilmn_map$Gene_Symbol)]
rownames(assayData) = assayData_ILMN_ids

# Update LumiBatch object
LumiBatch_v = LumiBatch_v[rownames(LumiBatch_v) %in% assayData_ILMN_ids,]
exprs(LumiBatch_v) = as.matrix(assayData)

# # Filter probes to one-to-one relation to ilmn ids
# LumiBatch_v = LumiBatch_v[rownames(LumiBatch_v) %in% unique_gene_mapping$ILMN_ID,]

# LumiBatch object
pData(LumiBatch_v) = phenoData


save(LumiBatch_v, file='Data/LumiBatch_Voineagu_big.RData')


remove(assayData, collapse_rows, gene_ilmn_map, heebo_ilmn_map, assayData_genes, 
       assayData_ILMN_ids, phenoData, chipArray_to_gsm)
