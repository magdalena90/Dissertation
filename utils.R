
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
library(topGO)

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

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

remove_faulty_probes = function(assayData, phenoData, heebo_gene_map){
  # Input:  assayData with its phenotype data and the mapping between heebo IDs and genes
  # Output: assayData with faulty probes removed
  
  # Get expression data and filter probes with repeats for their gene
  exprs_all = assayData[,grepl('exprs', colnames(assayData))]
  rptd_genes = names(table(heebo_gene_map$Gene_Symbol)[table(heebo_gene_map$Gene_Symbol)>1])
  exprs = exprs_all[rownames(exprs_all) %in% 
                      heebo_gene_map$ID[heebo_gene_map$Gene_Symbol %in% rptd_genes],]
  
  # Indicate if value is <1, group by age (horizontal sum)
  exprs_ind = data.frame(matrix(0, ncol = ncol(exprs), nrow = nrow(exprs)))
  exprs_ind[exprs <= 1] = 1
  exprs_ind_grp = data.frame(t(aggregate(t(exprs_ind), list(phenoData$age_group), sum)))
  colnames(exprs_ind_grp) = exprs_ind_grp[1,]
  exprs_ind_grp = exprs_ind_grp[-1,]
  exprs_ind_grp = data.frame(apply(exprs_ind_grp, 2, function(x) as.numeric(as.character(x))))
  rownames(exprs_ind_grp) = rownames(exprs)
  
  # Group by gene (vertical sum)
  sorted_heebo_gene_map = heebo_gene_map[match(rownames(exprs_ind_grp), heebo_gene_map$ID),]
  exprs_ind_grp_genes = aggregate(exprs_ind_grp, list(sorted_heebo_gene_map$Gene_Symbol), sum)
  rownames_ = exprs_ind_grp_genes[,1]
  exprs_ind_grp_genes = exprs_ind_grp_genes[,-1]
  exprs_ind_grp_genes = data.frame(apply(exprs_ind_grp_genes, 2, 
                                         function(x) as.numeric(as.character(x))))
  rownames(exprs_ind_grp_genes) = rownames_
  max_rpt_by_gene = apply(exprs_ind_grp_genes, 1, max)
  error_genes = rownames(exprs_ind_grp_genes)[max_rpt_by_gene==1]
  
  # Keep exprs of error genes, find error probes
  candidate_error_probes = heebo_gene_map$ID[heebo_gene_map$Gene_Symbol %in% error_genes]
  exprs_ind_grp_cep = exprs_ind_grp[rownames(exprs_ind_grp) %in% candidate_error_probes,]
  max_rpt_by_probe = apply(exprs_ind_grp_cep, 1, max)
  error_probes = rownames(exprs_ind_grp_cep)[max_rpt_by_probe==1]
  
  # Filter assayData
  assayData = assayData[!rownames(assayData) %in% error_probes,]
  
  return(assayData)
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
    pearson_cor = cor(df, method = 'pearson')
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

prepare_pca_c_vs_v=function(exprs_c, exprs_v, pData_c, pData_v, title, folder, file_name, by='Source'){
  
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
  pca = perform_pca(exprs_c_v, labels, title, folder, file_name)
  
  return(pca)
}

perform_pca = function(df, labels, title, folder, file_name){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: pca object and plot of 2 principal components coloured by label
  
  # Prepare data
  pca_data = prepare_visualisation_data(df, labels, vis='PCA')
  
  # Perform pca and plot
  pca_data_vals = pca_data[,!names(pca_data) %in% colnames(labels)]
  pca = prcomp(pca_data_vals, center = TRUE, scale. = TRUE)
  pca_plot = autoplot(pca, data = pca_data, colour = colnames(labels)[2], alpha = 0.7) + 
    theme_minimal() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5), 
                                             legend.position='bottom') #+
    #scale_color_manual(values=c('#7324A6','#f43e6c','#FFDB00')) # Brain region
    #scale_color_manual(values=c('#009999','#006666','#99cc00')) # Autism
  print(pca_plot)
  ggsave(paste0('Plots/PCA/',folder,'/',file_name,'.png'), pca_plot, width=5.54, height=3.66,
        units = 'in')
  
  return(pca)
}

perform_mds = function(df, labels, title){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - labels: Dataframe with columns ID and label with each row corresponding to a sample
  # Output: mds object and 2D plot coloured by label
  
  # Prepare data
  mds_data = prepare_visualisation_data(df, labels, vis = 'MDS')

  # Plot MDS
  mds_plot = ggplot(mds_data, aes(x=X1, y=X2, colour=get(colnames(labels)[2]))) + theme_minimal() + 
    geom_point(alpha = 0.7) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5), 
    legend.position='bottom') + labs(colour = colnames(labels)[2])
  print(mds_plot) 
  
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
    lda_fit = lda(x=lda_data[,-ncol(lda_data)], grouping = lda_data[,ncol(lda_data)])
    lda_points = as.matrix(lda_data[,-ncol(lda_data)]) %*% lda_fit$scaling
  }
  lda_points = data.frame(cbind(lda_points[,1:2], as.character(lda_data[,ncol(lda_data)])))
  colnames(lda_points) = c('LD1','LD2',colnames(labels)[2])
  
  if(colnames(labels[2])=='age_group'){
    ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
    lda_points = transform(lda_points, age_group=factor(age_group, levels = ordered_levels))
  }
  
  # Plot LDA
  xlab = paste0('LD1 (', round(100*lda_fit$svd[1]^2/sum(lda_fit$svd^2),1),'%)') # between group var
  ylab = paste0('LD2 (', round(100*lda_fit$svd[2]^2/sum(lda_fit$svd^2),1),'%)')
  
  lda_plot = ggplot(lda_points, aes(x=LD1, y=LD2, colour=get(colnames(labels)[2]))) + geom_point() + 
    ggtitle(title) + labs(colour = colnames(labels)[2]) + xlab(xlab) + ylab(ylab) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank(), 
          legend.position='bottom',
          axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  print(lda_plot)
  
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
  
  tsne_plot = ggplot(tsne_points, aes(x=x, y=y, colour=get(colnames(labels)[2]))) + geom_point() + 
    ggtitle(title) + labs(colour = colnames(labels)[2]) + theme(plot.title = element_text(hjust=0.5),
      axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
      axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    theme_minimal()
  print(tsne_plot)

  return(tsne_points)  
}

identify_significant_genes = function(exprs, pData, source, obj.var, fdr = 0.05){
  # Input: - df:     Dataframe with probes as rows and samples as columns. Needs rownames
  #        - pData:  pData dataframe
  # Output: data frame with statistically significant genes
  
  if(source == 'Colantuoni'){
    if(obj.var=='age_group'){
      resp.type = 'Multiclass'
      obj.var_vals = as.factor(pData[[obj.var]])
      if(length(unique(pData$age_group)) == 10){
        ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
      } else {
        #ordered_levels = c('Fetal','Infant_Child','10-20','20s_50s','60s_70s')
        ordered_levels = c('Fetal','Infant_Child','10-20','20s_30s','40s_50s','60s_70s') 
      }
      obj.var_vals = factor(obj.var_vals, ordered = TRUE, levels = ordered_levels)
      } else {
        resp.type = 'Quantitative'
        obj.var_vals = pData[[obj.var]]
      }
  } else {
    resp.type = 'Two class unpaired'
    obj.var_vals = as.factor(pData[[obj.var]])
  }
  
  # Perform SAM and extract significant genes
  SAM_fit = SAM(exprs, obj.var_vals, resp.type = resp.type, geneid = rownames(exprs),
                        random.seed = 1234, logged2 = TRUE, fdr.output = fdr)
  signif_genes = data.frame(SAM_fit$siggenes.table$genes.up)
  print(paste0(SAM_fit$siggenes.table$ngenes.up, ' upregulated genes'))
  print(paste0(SAM_fit$siggenes.table$ngenes.do, ' downregulated genes'))
  if('Fold.Change' %in% colnames(signif_genes)){
    signif_genes = signif_genes[as.numeric(as.character(signif_genes$Fold.Change))>1.3,]  
  }
  
  signif_genes$Score.d. = as.numeric(as.character(signif_genes$Score.d.))
  print(paste0(nrow(signif_genes), ' statistically significant genes found'))
  
  return(signif_genes)
}

GO_enrichment_analysis = function(genes, modules, module, ontology='BP'){
  # Input: - genes:     vector with the ILMN IDs of all the genes
  #        - modules:   vector of modules each gene belongs to
  #        - module:    module on which perform the enrichment analysis
  #        - ontology:  BP: Biological process / CC: Cellular component / MF: Molecular function
  # Output:  data frame with top GO terms for the selected ontology and some details and statistics
  
  # Create topGO object
  mod_n_genes = function (allModules) { return(allModules == 1) }
  geneList = as.numeric(modules==module)
  names(geneList) = genes
  GOdata = new('topGOdata', ontology = ontology, allGenes = geneList, geneSel = mod_n_genes,
                nodeSize = 10, annot = annFUN.db, affyLib = 'illuminaHumanv4.db')
  
  # Perform statistical tests
  weight01_fisher = runTest(GOdata, statistic = 'fisher')
  classic_fisher = runTest(GOdata, algorithm = 'classic', statistic = 'fisher')
  elim_fisher = runTest(GOdata, algorithm = 'elim', statistic = 'fisher')
  
  all_res = GenTable(GOdata, classicFisher = classic_fisher, weight01Fisher = weight01_fisher,
    elimFisher = elim_fisher, orderBy = 'elimFisher', ranksOf = 'classicFisher', topNodes = 10)
  
  # Plot significant genes in the GO tree
  print(showSigOfNodes(GOdata, score(elim_fisher), firstSigNodes = 5, useInfo = 'all'))
  
  return(list(all_res, GOdata))
  
}
