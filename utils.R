
library(GEOquery)
library(AnnotationDbi)
library(illuminaHumanv4.db)
library(tibble)
library(data.table)
library(ggfortify)
library(MASS)        # LDA
library(plsgenomics) # PLS.LDA
library(caret)
library(samr)

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
  Voineagu_suppMat = read.csv('Data/Voineagu/aux/supplementary_table_2.csv')
  Voineagu_suppMat$chip_array = paste(Voineagu_suppMat$Chip, Voineagu_suppMat$Array, sep='_')
  Voineagu_suppMat = Voineagu_suppMat[,c('Brain.Bank.Case.Id','Cortex.area','Disease.status',
                                         'chip_array')]
  Voineagu_suppMat$Cortex.area = as.character(Voineagu_suppMat$Cortex.area)
  
  # MATCH ALL FRONTAL + TEMPORAL IDS:
  # GSM, BrainBank (relation obtained from the series matrix phenotype information)
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
  samples_metadata = read.delim2('Data/Voineagu/aux/supplementary_table_1.tsv')
  
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

prepare_visualisation_data = function(df, labels, vis='PCA'){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  #        - vis:    Visualisation (PCA, MDS)
  # Output: data ready for visualisation analysis
  
  # Remove rows with NAs
  df = na.omit(df)
  labels = labels[labels$ID %in% colnames(df),]
  
  # Visualisation specific transformations
  if(vis=='MDS'){
    pearson_cor = cor(exprs_c, method = 'pearson')
    distance_mat = 1 - pearson_cor
    rownames(distance_mat) = colnames(distance_mat)
    mds = cmdscale(distance_mat, k = 2, eig = TRUE)
    
    vis_data = data.frame(mds$points)
    rownames(vis_data) = rownames(distance_mat)
  } else {
    vis_data = data.frame(t(df))
    vis_data[] = lapply(lapply(vis_data, as.character), as.numeric)
  }
  
  # BoxCox, center and rescale data
  trans = preProcess(vis_data, c('BoxCox', 'center', 'scale'))
  vis_data = data.frame(predict(trans, vis_data))
  
  # Merge vis_data with labels
  vis_data$ID = rownames(vis_data)
  vis_data = merge(vis_data, labels, by='ID')
  vis_data$ID = NULL
  if(colnames(labels[2])=='age_group'){
    ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
    vis_data = transform(vis_data, age_group=factor(age_group, levels = ordered_levels))
  }
  
  return(vis_data)
}

prepare_pca_c_vs_v = function(exprs_c, exprs_v, pData_c, pData_v, title, by='Source'){
  exprs_c_v = merge(exprs_c, exprs_v, by='row.names')
  rownames(exprs_c_v) = exprs_c_v$Row.names
  exprs_c_v$Row.names = NULL
  if(by == 'Source'){
    labels = rbind(data.frame('ID' = rownames(pData_c), 'Source' = 'Colantuoni'),
                   data.frame('ID' = rownames(pData_v), 'Source' = 'Voineagu'))
  } else {
    labels = rbind(data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group),
                   data.frame('ID' = rownames(pData_v), 'age_group' = pData_v$age_group))
  }
  pca = perform_pca(exprs_c_v, labels, title)
  
  return(pca)
}

perform_pca = function(df, labels, title){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: pca object and plot of 2 principal components coloured by label
  
  # Prepare data
  pca_data = prepare_visualisation_data(df, labels, vis='PCA')
  
  # Perform pca and plot
  pca_data_vals = pca_data[,!names(pca_data) %in% colnames(labels)]
  pca = prcomp(pca_data_vals, center = TRUE, scale. = TRUE)
  print(autoplot(pca, data = pca_data, colour = colnames(labels)[2], alpha = 0.7) + 
        theme_minimal() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)))
  
  return(pca)
}

perform_mds = function(df, labels, title){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: mds object and 2D plot coloured by label
  
  # Prepare data
  mds_data = prepare_visualisation_data(df, labels, vis = 'MDS')

  # Plot MDS
  print(ggplot(mds_data, aes(x=X1, y=X2, colour=get(colnames(labels)[2]))) + geom_point() + 
        theme_minimal() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
        labs(colour = colnames(labels)[2]))
  
  return(mds_data)
}

perform_lda = function(df, labels, title, pls=FALSE){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: lda object and 2D plot coloured by label
  
  # Prepare data
  lda_data = prepare_visualisation_data(df, labels, vis = 'LDA')
  
  # Perform LDA
  if(pls){
    pls_fit = pls.lda(lda_data[,-ncol(lda_data)], lda_data[,ncol(lda_data)], ncom=1:4, nruncv=10)
    lda_fit = pls_fit$lda.out
    lda_points = as.matrix(lda_data[,-ncol(lda_data)]) %*% pls_fit$pls.out$R %*% lda_fit$scaling
  } else {
    lda_fit = lda(age_group ~ ., data = lda_data)
    lda_points = as.matrix(lda_data[,-ncol(lda_data)]) %*% lda_fit$scaling 
  }
  lda_points = data.frame(cbind(lda_points[,1:2], as.character(pData$age_group)))
  colnames(lda_points) = c('LD1','LD2',colnames(labels)[2])
  
  if(colnames(labels[2])=='age_group'){
    ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
    lda_points = transform(lda_points, age_group=factor(age_group, levels = ordered_levels))
  }
  
  # Plot LDA
  xlab = paste0('LD1 (', round(100*lda_fit$svd[1]^2/sum(lda_fit$svd^2),1),'%)') # between group var
  ylab = paste0('LD2 (', round(100*lda_fit$svd[2]^2/sum(lda_fit$svd^2),1),'%)')
  
  print(ggplot(lda_points, aes(x=LD1, y=LD2, colour=get(colnames(labels)[2]))) + geom_point() + 
          ggtitle(title) + labs(colour = colnames(labels)[2]) + xlab(xlab) + ylab(ylab) +
          theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank(), 
                axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
                axis.ticks.y=element_blank()))
  
  return(lda_fit)
}

perform_tsne = function(df, labels, perplexity, title){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: data frame with x,y coordinates and 2D plot
  
  # Prepare data
  tsne_data = prepare_visualisation_data(df, labels, vis='TSNE')
  
  # Perform tsne and plot
  tsne_data_vals = tsne_data[,!names(tsne_data) %in% colnames(labels)]
  ptm <- proc.time()
  tsne_points = tsne(tsne_data_vals, perplexity = perplexity)
  (proc.time() - ptm)/60
  
  tsne_points = data.frame(cbind(tsne_points, as.character(tsne_data[,ncol(tsne_data)])))
  
  colnames(tsne_points) = c('x','y',colnames(labels)[2])
  
  if(colnames(labels[2])=='age_group'){
    ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
    tsne_points = transform(tsne_points, age_group=factor(age_group, levels = ordered_levels))
  }
  
  print(ggplot(tsne_points, aes(x=x, y=y, colour=get(colnames(labels)[2]))) + geom_point() + 
          ggtitle(title) +  labs(colour = colnames(labels)[2]) +
          theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
                axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
                axis.text.y=element_blank(), axis.ticks.y=element_blank()) + theme_minimal())
  
  return(tsne_points)  
}

identify_significant_genes = function(exprs, labels, source, fdr = 0.05){
  if(source == 'Colantuoni'){
    obj.var = 'age_group'
    resp.type = 'Multiclass'
    sex_col = '`Sex:ch1`'
  } else {
    obj.var = 'Disease.status'
    resp.type = 'Two class unpaired'
    sex_col = 'SEX'
  }
  
  # Perform SAM and extract significant genes
  SAM_fit = SAM(exprs, as.factor(labels[[obj.var]]), resp.type=resp.type, 
                geneid=rownames(exprs), random.seed=1234, logged2=TRUE, fdr.output=fdr)
  signif_genes = data.frame(SAM_fit$siggenes.table$genes.up)
  signif_genes = signif_genes[as.numeric(as.character(signif_genes$Fold.Change))>1.3,]
  
  # Detect significant genes for gender identification and remove them from list
  SAM_fit_sex = SAM(exprs, as.factor(pData_v[[sex_col]]), resp.type=resp.type, 
                    geneid=rownames(exprs), random.seed=1234, logged2=TRUE, fdr.output=fdr)
  signif_genes_sex = data.frame(SAM_fit_sex$siggenes.table$genes.up)
  signif_genes_sex = signif_genes_sex[as.numeric(as.character(signif_genes_sex$Fold.Change))>1.3,]
  
  signif_genes = signif_genes[!signif_genes$Gene.Name %in% signif_genes_sex$Gene.Name,]
  
  return(signif_genes)
}


