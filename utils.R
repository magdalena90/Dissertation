
library(GEOquery)
library(AnnotationDbi)
library(illuminaHumanv4.db)
library(tibble)
library(data.table)
library(ggfortify)

map_HEEBO_ILMN_ids = function(){
  # Output: Dataframe with all the genes present in both datasets with their corresponding ids
  #         Dataframe columns: Gene_Symbol, HEEBO_ID and ILMN_ID
  
  # COLANTUONI (HEEBO):
  gpl = getGEO(filename='Data/Colantuoni/GPL4611_family_short.soft', getGPL = FALSE)
  heebo_gene_map = data.frame(Table(gpl)[,c('ID','Gene_Symbol')])
  
  # Remove NAs and numeric values in genes:
  heebo_gene_map = heebo_gene_map[complete.cases(heebo_gene_map),]
  heebo_gene_map = heebo_gene_map[is.na(as.numeric(heebo_gene_map$Gene_Symbol)),]
  colnames(heebo_gene_map)[1] = 'HEEBO_ID'
  
  # VOINEAGU (ILMN):
  gse = read.delim('~/MSc/Dissertation/Data/Voineagu/GSE28521_non-normalized_data.txt.gz')
  ilmn_gene_map = data.frame(mapIds(illuminaHumanv4.db, key=as.character(gse$PROBE_ID), 
                                    column=c('SYMBOL'), keytype='PROBEID'))
  ilmn_gene_map = rownames_to_column(ilmn_gene_map, var='ILMN_ID')
  colnames(ilmn_gene_map)[2] = 'Gene_Symbol'
  
  # INNER MERGE:
  heebo_ilmn_map = merge(heebo_gene_map, ilmn_gene_map, by='Gene_Symbol')
  heebo_ilmn_map = heebo_ilmn_map[heebo_ilmn_map$Gene_Symbol!='',]
  heebo_ilmn_map_wo_ctrls = heebo_ilmn_map[!grepl('CTRL', heebo_ilmn_map$HEEBO_ID),]
  
  return(list(heebo_ilmn_map,heebo_ilmn_map_wo_ctrls))
}

create_Colantuonis_raw_data_df = function(cols){
  # Input: list of column names to extract
  # Output: Dataframe with {sample}_{col} as columns and probes as rows
  
  folder_dir = 'Data/Colantuoni/GSE30272_RAW/'
  files = list.files(folder_dir)
  
  # Initialise dataframe
  sample_data = fread(paste0('zcat ', folder_dir, files[1]), sep='\t')
  df = data.frame(HEEBO_ID=sample_data$PlatePos)
  
  # Fill input columns for each sample
  for(i in 1:length(files)){
    sample_data = fread(paste0('zcat ', folder_dir, files[i]), sep='\t')
    col_name_prefix = substr(files[i], 1, nchar(files[i])-7)
    for(j in 1:length(cols)){
      col_name = paste(col_name_prefix, cols[j], sep='_')
      df[[col_name]] = sample_data[[cols[j]]]
    }
  }
  return(df)
}

perform_pca = function(df, labels, title){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs an ID column
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: pca object and plot of 2 principal components coloured by label
  
  # Transpose and set colnames of dataframe
  pca_data = data.frame(t(df))
  colnames(pca_data) = as.character(unlist(pca_data[1,]))
  pca_data = pca_data[-1,]
  pca_data[] = lapply(lapply(pca_data, as.character), as.numeric)
  pca_data = log(pca_data+1)
  pca_data = pca_data[!duplicated(pca_data),]
  
  # Merge with labels df
  pca_data$ID = rownames(pca_data)
  pca_data = merge(pca_data, labels, by='ID')
  
  # Perform pca and plot
  pca_data_vals = pca_data[,!names(pca_data) %in% colnames(labels)]
  pca = prcomp(pca_data_vals, center = TRUE, scale. = TRUE)
  print(autoplot(pca, data=pca_data, colour=colnames(labels)[2]) + theme_minimal() + ggtitle(title))
  
  return(pca)
}
