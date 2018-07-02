
library(GEOquery)
library(AnnotationDbi)
library(illuminaHumanv4.db)
library(tibble)
library(data.table)
library(ggfortify)

map_HEEBO_ILMN_ids = function(){
  # Output: Dataframe with all the genes present in both datasets with their corresponding ids
  #         Dataframe columns: Gene_Symbol, HEEBO_ID and ILMN_ID
  # Note: IDs that map to more than one gene are filtered
  
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
                                    column=c('SYMBOL'), keytype='PROBEID', multiVals='filter'))
  ilmn_gene_map = rownames_to_column(ilmn_gene_map, var='ILMN_ID')
  colnames(ilmn_gene_map)[2] = 'Gene_Symbol'
  
  # INNER MERGE:
  heebo_ilmn_map = merge(heebo_gene_map, ilmn_gene_map, by='Gene_Symbol')
  heebo_ilmn_map = heebo_ilmn_map[heebo_ilmn_map$Gene_Symbol!='',]
  heebo_ilmn_map_wo_ctrls = heebo_ilmn_map[!grepl('CTRL', heebo_ilmn_map$HEEBO_ID),]
  
  return(list(heebo_ilmn_map,heebo_ilmn_map_wo_ctrls))
}

create_Colantuonis_raw_data_df = function(cols, colnames, sep='_'){
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
      col_name = paste(col_name_prefix, colnames[j], sep=sep)
      df[[col_name]] = sample_data[[cols[j]]]
    }
  }
  return(df)
}

create_Voineagus_pData_df = function(){
  # Output: Dataframe with pData of Voineagu's raw data
 
  # BrainBank, chip_array
  Voineagu_suppMat = read.csv('Data/Voineagu/aux/Voineagu_suppMat.csv')
  Voineagu_suppMat$chip_array = paste(Voineagu_suppMat$Chip, Voineagu_suppMat$Array, sep='_')
  Voineagu_suppMat = Voineagu_suppMat[,c('Brain.Bank.Case.Id','Cortex.area','Disease.status',
                                         'chip_array')]
  Voineagu_suppMat$Cortex.area = as.character(Voineagu_suppMat$Cortex.area)
  
  # MATCH ALL FRONTAL + TEMPORAL IDS:
  # GSM, BrainBank
  ids_mapping = read.delim('Data/Voineagu/aux/ids_mapping.txt', header=FALSE)
  colnames(ids_mapping)[1:2] = c('GSM','Brain.Bank.Plus')
  id_split = strsplit(as.character(ids_mapping$Brain.Bank.Plus),'_')
  ids_mapping$Disease.status = sapply(id_split, `[`, 1)
  ids_mapping$Brain.Bank.Case.Id = sapply(id_split, `[`, 2)
  ids_mapping$Cortex.area = sapply(id_split, `[`, 3)
  
  cortex_area = list('C'='cerebellum','T'='temporal','F'='frontal')
  ids_mapping$Cortex.area[ids_mapping$Cortex.area %in% names(cortex_area)] = 
    cortex_area[ids_mapping$Cortex.area[ids_mapping$Cortex.area %in% names(cortex_area)]]
  
  disease_status = list('A'='autism', 'C'='control')
  ids_mapping$Disease.status[ids_mapping$Disease.status %in% names(disease_status)] = 
    disease_status[ids_mapping$Disease.status[ids_mapping$Disease.status %in% names(disease_status)]]
  
  # BrainBank, metadata
  samples_metadata = read.delim2('Data/Voineagu/aux/samples_metadata.tsv')
  
  # Frontal and temporal info:
  ids_mapping_f_t = merge(Voineagu_suppMat, ids_mapping, 
                          by=c('Brain.Bank.Case.Id','Disease.status','Cortex.area'))
  
  # MATCH ALL CEREBELLUM IDS:
  # GSM, BrainBank, chip_array (cerebellum)
  ids_mapping_c = read.delim('Data/Voineagu/aux/cerebellum_ids_mapping.txt', header=FALSE)
  colnames(ids_mapping_c) = c('GSM','Brain.Bank.Plus','chip_array')
  
  id_split = strsplit(as.character(ids_mapping_c$Brain.Bank.Plus),'_')
  ids_mapping_c$Disease.status = sapply(id_split, `[`, 1)
  ids_mapping_c$Brain.Bank.Case.Id = sapply(id_split, `[`, 2)
  ids_mapping_c$Cortex.area = sapply(id_split, `[`, 3)
  
  ids_mapping_c$Cortex.area[ids_mapping_c$Cortex.area %in% names(cortex_area)] = 
    cortex_area[ids_mapping_c$Cortex.area[ids_mapping_c$Cortex.area %in% names(cortex_area)]]
  
  ids_mapping_c$Disease.status[ids_mapping_c$Disease.status %in% names(disease_status)] = 
  disease_status[ids_mapping_c$Disease.status[ids_mapping_c$Disease.status %in% names(disease_status)]]
  
  # JOIN FRONTAL+TEMPORAL INFO WITH CEREBELLUM AND ADD METADATA:
  ids_mapping_all = rbind(ids_mapping_f_t, ids_mapping_c)
  ids_mapping_all = ids_mapping_all[,c('GSM','chip_array','Brain.Bank.Case.Id','Cortex.area',
                                       'Disease.status')]
  
  samples_full_data = merge(ids_mapping_all, samples_metadata, by='Brain.Bank.Case.Id')
  samples_full_data$Cortex.area = as.character(samples_full_data$Cortex.area)
  
  # Add age group data
  age_group = cut(as.numeric(samples_full_data$AGE), c(-0.5, 0, 0.6, 10, 20, 30, 40, 50, 60, 70, 80),
                  labels=c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s'))
  samples_full_data$age_group = as.character(age_group)
  
  return(samples_full_data)
}

perform_pca = function(df, labels, title){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs an ID column
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: pca object and plot of 2 principal components coloured by label
  
  # Transpose and set colnames of dataframe
  pca_data = data.frame(t(df))
  pca_data[] = lapply(lapply(pca_data, as.character), as.numeric)
  pca_data = log(pca_data+1)
  pca_data = pca_data[!duplicated(pca_data),]   # Check
  
  # Merge with labels df
  pca_data$ID = rownames(pca_data)
  pca_data = merge(pca_data, labels, by='ID')
  if(colnames(labels[2])=='age_group'){
    ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
    pca_data = transform(pca_data, age_group=factor(age_group, levels = ordered_levels))
  }
  
  # Perform pca and plot
  pca_data_vals = pca_data[,!names(pca_data) %in% colnames(labels)]
  pca = prcomp(pca_data_vals, center = TRUE, scale. = TRUE)
  print(autoplot(pca, data = pca_data, colour = colnames(labels)[2], alpha = 0.7) + 
        theme_minimal() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)))
  
  return(pca)
}
