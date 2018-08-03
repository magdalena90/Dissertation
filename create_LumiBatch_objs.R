
setwd('~/MSc/Dissertation')
source('utils.R')

library(lumi)
library(GEOquery)
library(WGCNA)

collapse_method = 'max_var_between' # 'avg', 'max_var', 'max_var_between'
######################################################################################################
######################################################################################################
# COLANTUONI

# Phenotype data
gse = getGEO(filename = 'Data/Colantuoni/GSE30272_series_matrix.txt.gz', getGPL = FALSE)
phenoData = pData(gse)
age_group = cut(as.numeric(phenoData$`age:ch1`), c(-0.5, 0, 0.6, 10, 20, 30, 40, 50, 60, 70, 80),
                labels=c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s'))
phenoData$age_group = as.character(age_group)

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

# Drop probes with entries <=1 that other repeats of the same gene and age group don't have
assayData = remove_faulty_probes(assayData, phenoData, heebo_gene_map)

# Collapse probe values corresponding to a same gene
assayData_genes = heebo_gene_map$Gene_Symbol[match(rownames(assayData), heebo_gene_map$ID)]
if(collapse_method == 'avg'){
  collapse_rows = collapseRows(assayData, rowGroup=assayData_genes, rowID=rownames(assayData),
                               method = 'Average')
  assayData = as.data.frame(collapse_rows$datETcollapsed)
}else if(collapse_method == 'max_var'){
  collapse_rows = collapseRows(assayData, rowGroup=assayData_genes, rowID=rownames(assayData),
                               method = 'maxRowVariance')
  assayData = as.data.frame(collapse_rows$datETcollapsed)
} else {
  exprs_all = assayData[,grepl('exprs', colnames(assayData))]
  exprs_aggr = t(aggregate(t(exprs_all), list(phenoData$age_group), mean))
  exprs_aggr = exprs_aggr[-1,]
  collapse_rows = collapseRows(exprs_aggr, rowGroup=assayData_genes, rowID=rownames(exprs_aggr), 
                               method = 'maxRowVariance')
  assayData = assayData[collapse_rows$selectedRow,]
  rownames(assayData) = unique(assayData_genes)
  
}

# Filter probes corresponding to genes not present in Voineagu's experiment
heebo_ilmn_map = data.frame(map_HEEBO_ILMN_ids()[1])
assayData = assayData[rownames(assayData) %in% heebo_ilmn_map$Gene_Symbol,]

# Change IDs to ILMN notation (keep only 1st ILMN ID for each gene)
gene_ilmn_map = unique(heebo_ilmn_map[,c('Gene_Symbol','ILMN_ID')])
gene_ilmn_map = gene_ilmn_map[!duplicated(gene_ilmn_map$Gene_Symbol),]
assayData_ILMN_ids = gene_ilmn_map$ILMN_ID[match(rownames(assayData), gene_ilmn_map$Gene_Symbol)]
rownames(assayData) = assayData_ILMN_ids

# Create exprs and se.exprs dataframes
exprs = assayData[,grepl('exprs', colnames(assayData))]
se.exprs = assayData[,grepl('se', colnames(assayData))]
colnames(exprs) = sapply(colnames(exprs), function(x) unlist(strsplit(x, '\\.'))[1])
colnames(se.exprs) = sapply(colnames(se.exprs), function(x) unlist(strsplit(x, '\\.'))[1])

# LumiBatch object
LumiBatch_c = new('LumiBatch', exprs = as.matrix(exprs), se.exprs = as.matrix(se.exprs))
pData(LumiBatch_c) = phenoData

save(LumiBatch_c, file=paste0('Data/RDatas/',collapse_method,'/LumiBatch_Colantuoni.RData'))


remove(collapse_rows, gpl, heebo_gene_map, assayData, exprs, se.exprs, gse, assayData_genes, 
       assayData_ILMN_ids, phenoData, age_group, exprs_aggr, exprs_all)
######################################################################################################
######################################################################################################
# VOINEAGU

# AssayData data
LumiBatch_v = lumiR('Data/Voineagu/GSE28521_non-normalized_data.txt.gz')

# Phenotype data
phenoData = create_Voineagus_pData_df()
phenoData = phenoData[match(colnames(exprs(LumiBatch_v)), phenoData$chip_array),]
phenoData$SEX[phenoData$SEX=='5'] = 'F'
rownames(phenoData) = phenoData$GSM

# Change subject IDs from chip_array to GSM format
chipArray_to_gsm = phenoData$GSM
names(chipArray_to_gsm) = phenoData$chip_array
colnames(LumiBatch_v) = chipArray_to_gsm[colnames(LumiBatch_v)]

# Filter probes to match to Colantuoni's genes and average probes corresponding to the same gene
assayData = data.frame(exprs(LumiBatch_v))
assayData = assayData[rownames(assayData) %in% heebo_ilmn_map$ILMN_ID,]
assayData_genes = heebo_ilmn_map$Gene_Symbol[match(rownames(assayData), heebo_ilmn_map$ILMN_ID)]
collapse_rows = collapseRows(assayData, rowGroup=assayData_genes, rowID=rownames(assayData))
assayData = as.data.frame(collapse_rows$datETcollapsed)
assayData_ILMN_ids = gene_ilmn_map$ILMN_ID[match(rownames(assayData), gene_ilmn_map$Gene_Symbol)]
rownames(assayData) = assayData_ILMN_ids

# Update LumiBatch object
LumiBatch_v = LumiBatch_v[rownames(LumiBatch_v) %in% assayData_ILMN_ids,]
exprs(LumiBatch_v) = as.matrix(assayData)

# LumiBatch object
pData(LumiBatch_v) = phenoData

save(LumiBatch_v, file=paste0('Data/',collapse_method,'/LumiBatch_Voineagu.RData'))


remove(assayData, collapse_rows, gene_ilmn_map, heebo_ilmn_map, assayData_genes, 
       assayData_ILMN_ids, phenoData, chipArray_to_gsm)
