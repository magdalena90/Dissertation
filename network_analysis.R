
  setwd('~/MSc/Dissertation')
  
  library(WGCNA)
  library(lumi)
  library(xtable)
  
  folder = 'max_var_between'
  load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed_ft.RData', sep='/'))
  exprs = exprs(LumiBatch_v)
  pData = pData(LumiBatch_v)
  
  module_signif_test = function(MEs, obj_var, method='Welch'){
    if(method == 'Welch'){
      obj_var_vals = names(table(obj_var))
      p_vals = as.numeric(apply(MEs, 2, function(x){
        t.test(x[obj_var==obj_var_vals[1]], x[obj_var==obj_var_vals[2]], var.equal=FALSE)$p.value
      }))
    } else {
      p_vals = apply(MEs, 2, function(x) cor.test(x, obj_var)$p.value)
    }
    return(p_vals)
  }

######################################################################################################
######################################################################################################
# VOINEAGU'S DATASET 

######################################################################################################
# SEPARATE NETWORKS FOR DISEASE STATUS, FIND MODULES RELATED TO CORTEX AREA

atsm_results = data.frame(matrix(ncol=4, nrow=0))
ctrl_results = atsm_results
# colnames(atsm_results) = c('power','n_modules','signif_modules')

disease_status = c('autism','control')
disease_status_networks_module_significance = list()
for(i in 4:15){
  for(ds in disease_status){
    exprs_ds = exprs[,colnames(exprs) %in% rownames(pData[pData$Disease.status == ds,])]
    pData_ds = pData[pData$Disease.status == ds,]
    
    # !!!!! This function does not work if AnnotationDBi is loaded .......
    #best_power = pickSoftThreshold(t(exprs)) # autism: 4, control: 4 / a: 4 c: 4 / a: 4 c: 4
    best_power = i #best_power$powerEstimate
    print(best_power)
    network_ds = blockwiseModules(t(exprs_ds), power=best_power, numericLabels=TRUE,
                                  minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
    
    signif_cortex_area = module_signif_test(network_ds$MEs, pData_ds$Cortex.area)
    signif_sex = module_signif_test(network_ds$MEs, pData_ds$SEX)
    signif_age = module_signif_test(network_ds$MEs, pData_ds$AGE, 'Pearson')
    
    signif_ds = data.frame('Cortex_area'=signif_cortex_area, 'Gender'=signif_sex, 'Age'=signif_age)
    rownames(signif_ds) = colnames(network_ds$MEs)
    signif_ds[signif_ds>0.05/nrow(signif_ds)] = 'ns'
    
    n_tot = nrow(signif_ds)
    n_sig = nrow(signif_ds[signif_ds$Cortex_area!='ns' & signif_ds$Gender=='ns' & signif_ds$Age=='ns',])
    
    if(ds=='autism'){
      atsm_results = rbind(atsm_results, c(i, n_tot, n_sig))
    } else {
      ctrl_results = rbind(ctrl_results, c(i, n_tot, n_sig))
    }

    disease_status_networks_module_significance[[ds]] = signif_ds
  }
}

lapply(disease_status_networks_module_significance, function(x) xtable(x))

######################################################################################################
# TOGETHER

best_power = pickSoftThreshold(t(exprs)) # 4
network_5 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 1, numericLabels=TRUE, 
                           minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_6 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 2, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_7 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 3, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_8 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 4, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_9 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 5, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_10 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 6, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_11 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 7, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_12 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 8, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_13 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 9, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_14 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 10, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_15 = blockwiseModules(t(exprs), power=best_power$powerEstimate + 11, numericLabels=TRUE, 
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)

library(AnnotationDbi)
library(illuminaHumanv4.db)

networks = list(network_4, network_5, network_6, network_7, network_8, network_9, network_10, 
network_11, network_12, network_13, network_14, network_15)

together_results = data.frame(matrix(ncol=4, nrow=0))
for(i in 1:12){
  network = networks[[i]]
  gene_module_relation = data.frame('probe_ID'=rownames(exprs), 'Module'=network$colors)
  gene_module_relation$Gene = mapIds(illuminaHumanv4.db, key=as.character(rownames(exprs)),
                                     column=c('SYMBOL'), keytype='PROBEID', multiVals='filter')
  colnames(gene_module_relation) = c('probe_ID', 'Module', 'Gene')
  
  write.csv(gene_module_relation, file = 'Data/Voineagu/genes_by_module_p4.csv', row.names = FALSE)
  
  signif_disease_status = module_signif_test(network$MEs, pData$Disease.status)
  signif_age = module_signif_test(network$MEs, pData$AGE, 'Pearson')
  signif_sex = module_signif_test(network$MEs, pData$SEX)
  signif_pmi = module_signif_test(network$MEs, as.numeric(pData$PMI), 'Pearson')
  
  signif = data.frame('Disease status'=signif_disease_status, 'Age'=signif_age, 
                      'Gender'=signif_sex, 'PMI'=signif_pmi)
  rownames(signif) = colnames(network$MEs)
  signif[signif>0.05/nrow(signif)] = 'ns'
  
  n_tot = nrow(signif)
  n_sig = nrow(signif[signif$Disease.status!='ns' & signif$Gender=='ns' & signif$Age=='ns',])
  together_results = rbind(together_results, c(i+3, n_tot, n_sig))
}

xtable(signif)

genes_p4 = rownames(exprs)[network_4$colors==14]
genes_p5 = rownames(exprs)[network_5$colors==14]
genes_p6 = rownames(exprs)[network_6$colors %in% c(11,4)]
genes_p7 = rownames(exprs)[network_7$colors %in% c(11,14)]
genes_p8 = rownames(exprs)[network_8$colors %in% c(15,7)]
genes_p9 = rownames(exprs)[network_9$colors==16]
genes_p10 = rownames(exprs)[network_10$colors==10]
genes_p11 = rownames(exprs)[network_11$colors==12]
genes_p12 = rownames(exprs)[network_12$colors==14]
genes_p13 = rownames(exprs)[network_13$colors==16]
genes_p14 = rownames(exprs)[network_14$colors==16]
genes_p15 = rownames(exprs)[network_15$colors==14]
genes_per_p = list(genes_p4, genes_p5, genes_p6, genes_p7, genes_p8, genes_p9, genes_p10,
                   genes_p11, genes_p12, genes_p13, genes_p14, genes_p15)

signif_genes_mat = data.frame(matrix(ncol=12, nrow=nrow(exprs)))
for(i in 1:12){
  signif_genes_mat[,i] = as.numeric(rownames(exprs) %in% genes_per_p[[i]])
}
rownames(signif_genes_mat) = rownames(exprs)
colnames(signif_genes_mat) = c('P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14','P15')

dist_mat = as.matrix(dist(t(signif_genes_mat)))

dist_mat_sq = dist_mat^2

library(ComplexHeatmap)
rg_palette = colorRampPalette(c('green', 'green', 'black', 'red', 'red'))(n = 100)
Heatmap(dist_mat, column_title='Similarity of Significant Modules', col = rg_palette, 
        show_heatmap_legend = FALSE, show_column_dend = FALSE)

signif_genes = unique(c(genes_p11, genes_p12, genes_p13, genes_p14, genes_p15))

map_id_gene = data.frame(map_HEEBO_ILMN_ids()[1])
map_id_gene = map_id_gene[,-2]
map_id_gene = unique(map_id_gene)

signif_genes = map_id_gene$Gene_Symbol[map_id_gene$ILMN_ID %in% signif_genes]
signif_genes = c(signif_genes, c('','','',''))

a = data.frame('a' = signif_genes[1:31],'b'=signif_genes[32:62],'c'=signif_genes[63:93],
               'd'=signif_genes[94:124],'e'=signif_genes[125:155])

######################################################################################################
# Enrichment Analysis

# power = 4: M14
source('utils.R')
write.csv(signif_genes, file='Data/WCGNA_signif_genes_DS.csv', row.names = FALSE)
mod_14_ea = GO_enrichment_analysis(rownames(exprs), network$colors, 14)
mod_14_GenTable = data.frame(mod_20_ea[1])


remove(LumiBatch_v, exprs, exprs_ds, network, network_ds, pData, pData_ds, signif, ds, signif_age, 
       signif_cortex_area, signif_sex, signif_pmi, signif_disease_status)
######################################################################################################
######################################################################################################
# COLANTUONI'S DATA

load(paste('Data/RDatas', folder, 'LumiBatch_Colantuoni_preprocessed.RData', sep='/'))
exprs = exprs(LumiBatch_c)
pData = pData(LumiBatch_c)
pData$`Sex:ch1`[pData$`Sex:ch1`=='5'] = 'F'

best_power = pickSoftThreshold(t(exprs)) # 8
network_8 = blockwiseModules(t(exprs), power=best_power$powerEstimate, numericLabels=TRUE,
                           minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_9 = blockwiseModules(t(exprs), power=best_power$powerEstimate+1, numericLabels=TRUE,
                           minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_10 = blockwiseModules(t(exprs), power=best_power$powerEstimate+2, numericLabels=TRUE,
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_11 = blockwiseModules(t(exprs), power=best_power$powerEstimate+3, numericLabels=TRUE,
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_12 = blockwiseModules(t(exprs), power=best_power$powerEstimate+4, numericLabels=TRUE,
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_13 = blockwiseModules(t(exprs), power=best_power$powerEstimate+5, numericLabels=TRUE,
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_14 = blockwiseModules(t(exprs), power=best_power$powerEstimate+6, numericLabels=TRUE,
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
network_15 = blockwiseModules(t(exprs), power=best_power$powerEstimate+7, numericLabels=TRUE,
                             minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)

networks = list(network_8, network_9, network_10, network_11, network_12, network_13, network_14,
                network_15)

library(AnnotationDbi)
library(illuminaHumanv4.db)

age_results = data.frame(matrix(ncol=4, nrow=0))
for(i in 1:8){
  network = networks[[i]]
  
  network = network_8
  
  gene_module_relation = data.frame('probe_ID'=rownames(exprs), 'Module'=network$colors)
  gene_module_relation$Gene = mapIds(illuminaHumanv4.db, key=as.character(rownames(exprs)),
                                     column=c('SYMBOL'), keytype='PROBEID', multiVals='filter')
  colnames(gene_module_relation) = c('probe_ID', 'Module', 'Gene')
  
  # write.csv(gene_module_relation, file = 'Data/Colantuoni/genes_by_module.csv', row.names = FALSE)
  
  signif_age = module_signif_test(network$MEs, as.numeric(pData$`age:ch1`), 'Pearson')
  signif_sex = module_signif_test(network$MEs, pData$`Sex:ch1`)
  signif_pmi = module_signif_test(network$MEs, as.numeric(pData$`postmortem interval (pmi):ch1`), 'Prs')
  
  signif = data.frame('Age'=signif_age, 'Gender'=signif_sex, 'PMI'=signif_pmi)
  rownames(signif) = colnames(network$MEs)
  signif[signif>0.05/nrow(signif)] = 'ns'
  
  n_tot = nrow(signif)
  n_sig = nrow(signif[signif$Age!='ns' & signif$Gender=='ns',])
  age_results = rbind(age_results, c(i+3, n_tot, n_sig))
}

genes_p8_1 = rownames(exprs)[network_8$colors==7] # 7, 26, 22, 12
genes_p8_2 = rownames(exprs)[network_8$colors==26]
genes_p8_3 = rownames(exprs)[network_8$colors==22]
genes_p8_4 = rownames(exprs)[network_8$colors==12]
genes_p9_1 = rownames(exprs)[network_9$colors==22] # 22 13 25 8
genes_p9_2 = rownames(exprs)[network_9$colors==13]
genes_p9_3 = rownames(exprs)[network_9$colors==25]
genes_p9_4 = rownames(exprs)[network_9$colors==8]
genes_p10_1 = rownames(exprs)[network_10$colors==14] # 14 20 22 10
genes_p10_2 = rownames(exprs)[network_10$colors==20]
genes_p10_3 = rownames(exprs)[network_10$colors==22]
genes_p10_4 = rownames(exprs)[network_10$colors==10]
genes_p11_1 = rownames(exprs)[network_11$colors==15] # 15 21 22 26 10
genes_p11_2 = rownames(exprs)[network_11$colors==21]
genes_p11_3 = rownames(exprs)[network_11$colors==22]
genes_p11_4 = rownames(exprs)[network_11$colors==26]
genes_p11_5 = rownames(exprs)[network_11$colors==10]
genes_p12_1 = rownames(exprs)[network_12$colors==10] # 10 20 12 19
genes_p12_2 = rownames(exprs)[network_12$colors==20]
genes_p12_3 = rownames(exprs)[network_12$colors==12]
genes_p12_4 = rownames(exprs)[network_12$colors==19]
genes_p13_1 = rownames(exprs)[network_13$colors==11] # 11 16 19
genes_p13_2 = rownames(exprs)[network_13$colors==16]
genes_p13_3 = rownames(exprs)[network_13$colors==19]
genes_p14_1 = rownames(exprs)[network_14$colors==10] # 10 16 18
genes_p14_2 = rownames(exprs)[network_14$colors==16]
genes_p14_3 = rownames(exprs)[network_14$colors==18]
genes_p15_1 = rownames(exprs)[network_15$colors==10] # 10 14 16
genes_p15_2 = rownames(exprs)[network_15$colors==14]
genes_p15_3 = rownames(exprs)[network_15$colors==16]

genes_per_p = list(genes_p8_1,genes_p8_2,genes_p8_3,genes_p8_4,genes_p9_1,genes_p9_2,genes_p9_3,
                   genes_p9_4,genes_p10_1,genes_p10_2,genes_p10_3,genes_p10_4,genes_p11_1,
                   genes_p11_2,genes_p11_3,genes_p11_4,genes_p11_5,genes_p12_1,genes_p12_2,
                   genes_p12_3,genes_p12_4,genes_p13_1,genes_p13_2,genes_p13_3,genes_p14_1,
                   genes_p14_2,genes_p14_3,genes_p15_1,genes_p15_2,genes_p15_3)

signif_genes_mat = data.frame(matrix(ncol=30, nrow=nrow(exprs)))
for(i in 1:30){
  signif_genes_mat[,i] = as.numeric(rownames(exprs) %in% genes_per_p[[i]])
}
rownames(signif_genes_mat) = rownames(exprs)
colnames(signif_genes_mat) = c('P8_1','P8_2','P8_3','P8_4','P9_1','P9_2','P9_3','P9_4','P10_1',
  'P10_2','P10_3','P10_4','P11_1','P11_1','P11_2','P11_3','P11_4','P12_1','P12_2','P12_3','P12_4',
  'P13_1','P13_2','P13_3','P14_1','P14_2','P14_3','P15_1','P15_2','P15_3')

dist_mat = as.matrix(dist(t(signif_genes_mat)))
dist_mat_sq = dist_mat^2

library(ComplexHeatmap)
rg_palette = colorRampPalette(c('green', 'green', 'black', 'red', 'red'))(n = 100)
Heatmap(dist_mat, column_title='Similarity of Significant Modules', col = rg_palette, 
        show_heatmap_legend = FALSE, show_column_dend = FALSE)

signif_genes_1=unique(c(genes_p12_4, genes_p13_2, genes_p14_2, genes_p15_2, genes_p11_1, genes_p10_2))
signif_genes_2=unique(c(genes_p8_2, genes_p9_3, genes_p11_3))
signif_genes_3=unique(c(genes_p11_2, genes_p12_2, genes_p10_3, genes_p14_3, genes_p15_3, genes_p13_3, 
                        genes_p8_3, genes_p9_1))
signif_genes_4=unique(c(genes_p15_1, genes_p14_1, genes_p13_1, genes_p12_2))

map_id_gene = data.frame(map_HEEBO_ILMN_ids()[1])
map_id_gene = map_id_gene[,-2]
map_id_gene = unique(map_id_gene)

signif_genes_1 = map_id_gene$Gene_Symbol[map_id_gene$ILMN_ID %in% signif_genes_1]
signif_genes_2 = map_id_gene$Gene_Symbol[map_id_gene$ILMN_ID %in% signif_genes_2]
signif_genes_3 = map_id_gene$Gene_Symbol[map_id_gene$ILMN_ID %in% signif_genes_3]
signif_genes_4 = map_id_gene$Gene_Symbol[map_id_gene$ILMN_ID %in% signif_genes_4]

save(signif_genes_1,signif_genes_2,signif_genes_3,signif_genes_4, file='Data/signif_genes_age.RData')
