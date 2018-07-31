
  setwd('~/MSc/Dissertation')
  
  library(WGCNA)
  library(lumi)
  library(xtable)
  
  folder = 'collapsed_probes_max'
  load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed.RData', sep='/'))
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

disease_status = c('autism','control')
disease_status_networks_module_significance = list()
for(ds in disease_status){
  exprs_ds = exprs[,colnames(exprs) %in% rownames(pData[pData$Disease.status == ds,])]
  pData_ds = pData[pData$Disease.status == ds,]
  
  # !!!!! This function does not work if AnnotationDBi is loaded .......
  best_power = pickSoftThreshold(t(exprs))
  network_ds = blockwiseModules(t(exprs_ds), power=best_power$powerEstimate+10, numericLabels=TRUE, 
                                minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)
  
  signif_cortex_area = module_signif_test(network_ds$MEs, pData_ds$Cortex.area)
  signif_sex = module_signif_test(network_ds$MEs, pData_ds$SEX)
  signif_age = module_signif_test(network_ds$MEs, pData_ds$AGE, 'Pearson')
  
  signif_ds = data.frame('Cortex_area'=signif_cortex_area, 'Gender'=signif_sex, 'Age'=signif_age)
  rownames(signif_ds) = colnames(network_ds$MEs)
  signif_ds[signif_ds>0.05/nrow(signif_ds)] = 'ns'
  
  disease_status_networks_module_significance[[ds]] = signif_ds
  # 3 no, 4 no, 5 no, 6 si, 7 no, 8 no, 9 no, 10 si, 11 no, 12 no, 13 no
}

lapply(disease_status_networks_module_significance, function(x) xtable(x))

######################################################################################################
# TOGETHER

best_power = pickSoftThreshold(t(exprs)) # 3
network = blockwiseModules(t(exprs), power=best_power$powerEstimate, numericLabels=TRUE, 
                           minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)

library(AnnotationDbi)
library(illuminaHumanv4.db)
gene_module_relation = data.frame('probe_ID'=rownames(exprs), 'Module'=network$colors)
gene_module_relation$Gene = mapIds(illuminaHumanv4.db, key=as.character(rownames(exprs)),
                                   column=c('SYMBOL'), keytype='PROBEID', multiVals='filter')
colnames(gene_module_relation) = c('probe_ID', 'Module', 'Gene')

write.csv(gene_module_relation, file = 'Data/Voineagu/genes_by_module.csv', row.names = FALSE)

signif_disease_status = module_signif_test(network$MEs, pData$Disease.status)
signif_age = module_signif_test(network$MEs, pData$AGE, 'Pearson')
signif_sex = module_signif_test(network$MEs, pData$SEX)
signif_pmi = module_signif_test(network$MEs, as.numeric(pData$PMI), 'Pearson')

signif = data.frame('Disease status'=signif_disease_status, 'Age'=signif_age, 'Gender'=signif_sex,
                    'PMI'=signif_pmi)
rownames(signif) = colnames(network$MEs)
signif[signif>0.05/nrow(signif)] = 'ns'

xtable(signif)

# ME20, ME18, ME12

# ME20:
source('utils.R')
mod_20_ea = GO_enrichment_analysis(rownames(exprs), network$colors, 20)
mod_20_GenTable = data.frame(mod_20_ea[1])

mod_18_ea = GO_enrichment_analysis(rownames(exprs), network$colors, 18)
mod_18_GenTable = data.frame(mod_18_ea[1])

mod_12_ea = GO_enrichment_analysis(rownames(exprs), network$colors, 12)
mod_12_GenTable = data.frame(mod_12_ea[1])

remove(LumiBatch_v, exprs, exprs_ds, network, network_ds, pData, pData_ds, signif, ds, signif_age, 
       signif_cortex_area, signif_sex, signif_pmi, signif_disease_status)
######################################################################################################
######################################################################################################
# COLANTUONI'S DATA

load(paste('Data/RDatas', folder, 'LumiBatch_Colantuoni_preprocessed.RData', sep='/'))
exprs = exprs(LumiBatch_c)
pData = pData(LumiBatch_c)
pData$`Sex:ch1`[pData$`Sex:ch1`=='5'] = 'F'

best_power = pickSoftThreshold(t(exprs))
network = blockwiseModules(t(exprs), power=best_power$powerEstimate, numericLabels=TRUE, 
                           minModuleSize=40, mergeCutHeight=0.1, minCoreKME=0.7)

library(AnnotationDbi)
library(illuminaHumanv4.db)
gene_module_relation = data.frame('probe_ID'=rownames(exprs), 'Module'=network$colors)
gene_module_relation$Gene = mapIds(illuminaHumanv4.db, key=as.character(rownames(exprs)),
                                   column=c('SYMBOL'), keytype='PROBEID', multiVals='filter')
colnames(gene_module_relation) = c('probe_ID', 'Module', 'Gene')

write.csv(gene_module_relation, file = 'Data/Colantuoni/genes_by_module.csv', row.names = FALSE)

signif_age = module_signif_test(network$MEs, as.numeric(pData$`age:ch1`), 'Pearson')
signif_sex = module_signif_test(network$MEs, pData$`Sex:ch1`)
signif_pmi = module_signif_test(network$MEs, as.numeric(pData$`postmortem interval (pmi):ch1`), 'Prs')

signif = data.frame('Age'=signif_age, 'Gender'=signif_sex, 'PMI'=signif_pmi)
rownames(signif) = colnames(network$MEs)
signif[signif>0.05/nrow(signif)] = 'ns'

# Signif modules for age: ME2, ME5, ME5, ME14, ME10

xtable(signif)
