
setwd('~/MSc/Dissertation')
source('utils.R')

library(GEOquery)
library(lumi)

######################################################################################################
# Load data

# Colantuoni
raw_df_c = create_Colantuonis_raw_data_df('SR_Mean')
gse_c = getGEO(filename='Data/Colantuoni/GSE30272_series_matrix.txt.gz', getGPL = FALSE)

# Voineagu
gse_v = lumiR('Data/Voineagu/GSE28521_non-normalized_data.txt.gz')
raw_df_v = data.frame(exprs(gse_v))
sample_info_v = data.frame(read.csv('~/MSc/Dissertation/Data/Voineagu/samples_full_metadata.csv'))

# Colantuoni-Voineagu mapping
heebo_ilmn_maps = map_HEEBO_ILMN_ids()
heebo_ilmn_map = data.frame(heebo_ilmn_maps[1])

remove(gse_v, heebo_ilmn_maps)

######################################################################################################
# PCA COLANTUONI LABELLED BY AGE

# Filter dataframe
# 1. IDs corresponding to genes present in ILMN dataset:
# filtered_raw_df_c = raw_df_c[raw_df_c$HEEBO_ID %in% heebo_ilmn_map$HEEBO_ID,]
# 2. One-to-one relation on both datasets:
heebo_ids_unique_genes = heebo_ilmn_map$HEEBO_ID[!duplicated(heebo_ilmn_map$Gene_Symbol) & 
                                                 !duplicated(heebo_ilmn_map$Gene_Symbol, fromLast=T)]
filtered_raw_df_c = raw_df_c[raw_df_c$HEEBO_ID %in% heebo_ids_unique_genes,]

# Rename columns with probe IDs
new_colnames = sapply(colnames(filtered_raw_df_c), function(x) unlist(strsplit(x, '_'))[1])
colnames(filtered_raw_df_c) = as.character(new_colnames)
colnames(filtered_raw_df_c)[1] = 'HEEBO_ID'

# Create labels dataframe
sample_info = data.frame('ID'=gse_c$geo_accession, 'age'=as.numeric(as.character(gse_c$`age:ch1`)))
sample_info$age_group = as.factor(cut(sample_info$age, c(-0.5, 0, 0.6, 10, 20, 30, 40, 50, 60, 70, 80),
                        labels=c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')))
sample_info = sample_info[,c('ID','age_group')]

pca = perform_pca(filtered_raw_df_c, sample_info, 'Colantuoni raw data by age')

remove(heebo_ids_unique_genes, new_colnames, sample_info, pca)
######################################################################################################
# PCA VOINEAGU LABELLED BY BRAIN REGION

# Filter dataframe
# 1. IDs corresponding to genes present in ILMN dataset:
# filtered_raw_df_v = raw_df_v[rownames(raw_df_v) %in% heebo_ilmn_map$ILMN_ID,]
# 2. One-to-one relation on both datasets:
ilmn_ids_unique_genes = heebo_ilmn_map$ILMN_ID[!duplicated(heebo_ilmn_map$Gene_Symbol) & 
                                                 !duplicated(heebo_ilmn_map$Gene_Symbol, fromLast=T)]
filtered_raw_df_v = raw_df_v[rownames(raw_df_v) %in% ilmn_ids_unique_genes,]

# Rename columns with probe IDs
colnames_dict = as.character(sample_info_v$GSM)
names(colnames_dict) = sample_info_v$chip_array
new_colnames = sapply(colnames(filtered_raw_df_v), function(x) colnames_dict[substr(x,2,nchar(x))])
colnames(filtered_raw_df_v) = as.character(new_colnames)
filtered_raw_df_v = cbind(rownames(filtered_raw_df_v), filtered_raw_df_v)
colnames(filtered_raw_df_v)[1] = 'ILMN_ID'

# Create labels dataframe
sample_info = sample_info_v[,c('GSM','Cortex.area')]
colnames(sample_info)[1] = 'ID'

pca = perform_pca(filtered_raw_df_v, sample_info, 'Voineagu raw data')

remove(ilmn_ids_unique_genes, colnames_dict, new_colnames, sample_info)
######################################################################################################
# PCA VOINEAGU FILTERED BY BRAIN REGION AND LABELLED BY DISEASE STATUS

# Filter by brain region
brain_region = 'frontal'
filtered_raw_df_v_br = cbind(filtered_raw_df_v[,1],filtered_raw_df_v[,colnames(filtered_raw_df_v) %in% 
                                  sample_info_v$GSM[sample_info_v$Cortex.area==brain_region]])
colnames(filtered_raw_df_v_br)[1] = 'ILMN_ID'

sample_info = sample_info_v[,c('GSM','Disease.status')]
colnames(sample_info)[1] = 'ID'

pca = perform_pca(filtered_raw_df_v_br, sample_info, paste0('Voineagu raw data from ',brain_region,
                                                         ' region'))

######################################################################################################
# PCA COLANTUONI VS VOINEAGU

# Filter autism samples from Voineagu
filtered_raw_df_v = cbind(filtered_raw_df_v_br[,1],filtered_raw_df_v_br[,colnames(filtered_raw_df_v_br)
                          %in% sample_info_v$GSM[sample_info_v$Disease.status=='control']])
colnames(filtered_raw_df_v)[1] = 'ILMN_ID'

# Combine dataframes
heebo_gene = merge(filtered_raw_df_c, unique(heebo_ilmn_map[,c('Gene_Symbol','HEEBO_ID')]), 
                   by = 'HEEBO_ID')
heebo_gene$HEEBO_ID = NULL

heebo_gene_ilmn = merge(heebo_gene, unique(heebo_ilmn_map[,c('Gene_Symbol','ILMN_ID')]),
                        by = 'Gene_Symbol')
heebo_gene_ilmn = merge(heebo_gene_ilmn, filtered_raw_df_v, by='ILMN_ID')
heebo_gene_ilmn$ILMN_ID = NULL

#  Create labels dataframe
sample_info = data.frame('ID'=colnames(heebo_gene_ilmn), 'source'='Colantuoni')
sample_info$ID = as.character(sample_info$ID)
sample_info$source = as.character(sample_info$source)
sample_info = sample_info[-1,]
sample_info$source[sample_info$ID %in% colnames(filtered_raw_df_v)] = 'Voineagu'

pca = perform_pca(heebo_gene_ilmn, sample_info, 'Voineagu vs. Colantuoni raw data')
