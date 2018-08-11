
setwd('~/MSc/Dissertation')
source('utils.R')

library(ComplexHeatmap)
library(circlize)
library(xtable)

folder = 'max_var_between' # avg, max_var, max_var_between
load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed.RData', sep='/'))
exprs = exprs(LumiBatch_v)
pData = pData(LumiBatch_v)

######################################################################################################
# STATISTICALLY SIGNIFICANT GENES BY REGION

improve_df = function(df, map_id_gene){
  genes = map_id_gene$Gene_Symbol[match(df$Gene.Name, map_id_gene$ILMN_ID)]
  df$Gene.Name = genes
  df = df[,-1]
  if(ncol(df) == 6){
    colnames(df) = c('Gene','Score','Numerator','Denomin','Fold.Change','q.val')
  } else {
    colnames(df) = c('Gene','Score','Numerator','Denomin','q.val')
  }
  
  return(df)
}

map_id_gene = data.frame(map_HEEBO_ILMN_ids()[1])
map_id_gene = map_id_gene[,c('ILMN_ID','Gene_Symbol')]
map_id_gene = unique(map_id_gene)

filter_cereb = pData$Cortex.area=='cerebellum'
# 983 significant genes / 971 / 718
cortex_signif_genes = identify_significant_genes(exprs[,!filter_cereb], pData[!filter_cereb,], 
                                                 'Voineagu', 'Disease.status')
cortex_signif_genes_clean = improve_df(cortex_signif_genes, map_id_gene)
xtable(cortex_signif_genes_clean[1:20,])

# 29 significant genes / 67 / 4
cerebellum_signif_genes = identify_significant_genes(exprs[,filter_cereb], pData[filter_cereb,], 
                                                     'Voineagu', 'Disease.status')
cerebellum_signif_genes_clean = improve_df(cerebellum_signif_genes, map_id_gene)
xtable(cerebellum_signif_genes_clean[1:20,])

# 181 significant genes / 125 / 61
cerebellum_signif_genes_rel = identify_significant_genes(exprs[,filter_cereb],pData[filter_cereb,],
                                                     'Voineagu', 'Disease.status', 0.25)
cerebellum_signif_genes_rel_clean = improve_df(cerebellum_signif_genes_rel, map_id_gene)
xtable(cerebellum_signif_genes_clean[1:20,])

remove(filter_cereb, cerebellum_signif_genes)
######################################################################################################
# HIERARCHICAL CLUSTERING AND HEATMAP FOR 200 DE GENES

# 200 DE genes
top_200_DE_genes = cortex_signif_genes[with(cortex_signif_genes,order(-Score.d.)),]
top_200_DE_genes = top_200_DE_genes[1:200,]
cortex_samples = rownames(pData)[pData$Cortex.area!='cerebellum']
exprs_top_200 = data.frame(exprs[rownames(exprs) %in% top_200_DE_genes$Gene.Name, 
                                 colnames(exprs) %in% cortex_samples])
pData = pData[pData$Cortex.area!='cerebellum',]
pData$PMI = as.numeric(as.character(pData$PMI))

# Color palettes
rg_palette = colorRampPalette(c('green', 'black', 'red'))(n = 100)

colormap_disease = c('#009999', '#99cc00')
names(colormap_disease) = c('autism','control')

colormap_cortex_area = gg_color_hue(3)
names(colormap_cortex_area) = c('frontal','temporal','cerebellum')

colormap_age = gg_color_hue(10)
names(colormap_age) = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')

colormap_sex = gg_color_hue(2)
names(colormap_sex)= c('F', 'M')

# Heatmap
annotation_df = data.frame('Disease'=pData$Disease.status)
cols_annotation_df = list('Disease'=colormap_disease)
hm_top_annotation = HeatmapAnnotation(annotation_df, col=cols_annotation_df)

annotation_df = data.frame('Cortex_area'=pData$Cortex.area, 'Age_group'=pData$age_group,
                           'Sex'=pData$SEX, 'PMI'=pData$PMI)
cols_annotation_df = list('Cortex_area'=colormap_cortex_area, 'Age_group'=colormap_age,
                          'Sex'=colormap_sex, 'PMI'=circlize::colorRamp2(c(200, 400),
                           c('#FF0080', '#00FF80')))
hm_bottom_annotation = HeatmapAnnotation(annotation_df, col = cols_annotation_df)

Heatmap(exprs_top_200, name='Expression levels', column_title='Samples', row_title='Genes',
        show_row_names=F, show_column_names=F, cluster_rows=T, show_row_dend=F, col=rg_palette, 
        show_heatmap_legend=F, top_annotation=hm_top_annotation, bottom_annotation=hm_bottom_annotation,
        top_annotation_height=unit(2, 'mm'), bottom_annotation_height=unit(8, 'mm') )

remove(annotation_df, cols_annotation_df, hm_bottom_annotation, hm_top_annotation)
######################################################################################################
# LINEAR REGRESSION VS AGE, PMI AN SEX TO ASSESS SIGNIFICANCE IN RESIDUALS

exprs_cortex = exprs[,!filter_cereb]
pData_cortex = pData[!filter_cereb,]

residual_exprs = exprs_cortex
residual_exprs[] = 0

for(i in 1:nrow(residual_exprs)){
  df = data.frame('age' = pData_cortex$AGE, 'sex' = as.factor(pData_cortex$SEX), 
                  'PMI' = as.numeric(as.character(pData_cortex$PMI)), 'exprs' = exprs_cortex[i,])
  df$PMI[is.na(df$PMI)] = mean(df$PMI, na.rm=T)
  lm_fit = lm(exprs ~ age+sex+PMI+1, df)
  residual_exprs[i,] = lm_fit$residuals
}

# 133 / 125 / 322
res_sign_genes = identify_significant_genes(residual_exprs, pData_cortex, 'Voineagu', 'Disease.status')
res_sign_genes_clean = improve_df(res_sign_genes, map_id_gene)
xtable(res_sign_genes_clean[1:20,])

sum(res_sign_genes$Gene.Name %in% cortex_signif_genes$Gene.Name)/nrow(res_sign_genes)

# scatter plot of scores
a = data.frame('NewScore' = res_sign_genes$Score.d., 'Gene' = res_sign_genes$Gene.Name)
b = data.frame('Pos_before' = seq(1,nrow(cortex_signif_genes)), 'Gene' = cortex_signif_genes$Gene.Name)
c = merge(a,b,by='Gene')

ggplot(c, aes(x=Pos_before, y=NewScore, colour=NewScore)) + geom_bar(stat='identity') + theme_minimal()+ 
  theme(plot.title = element_text(hjust = 0.5)) + ggtitle('Scores of remaining DE genes') + 
  xlab('Gene ranking in original scoring') +  ylab('Score after correction')


ggplot(c, aes(x = Pos_before, y = Pos_after)) + theme_minimal() + geom_point(alpha = 0.7) + 
  geom_smooth(method = 'lm', formula = y ~ x) + theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Scores before and after Correcting by Confounding Vars') + 
  xlab('Score before correction') +  ylab('Score after correction')

######################################################################################################
######################################################################################################
# DE BY BRAIN REGION FOR AUTISM AND CONTROL SAMPLES SEPARATELY

folder = 'avg'
load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed_ft.RData', sep='/'))
exprs = exprs(LumiBatch_ft)
pData = pData(LumiBatch_ft)

exprs_atsm = exprs[,colnames(exprs) %in% rownames(pData)[pData$Disease.status == 'autism']]
pData_atsm = pData[pData$Disease.status == 'autism',]
exprs_ctrl = exprs[,colnames(exprs) %in% rownames(pData)[pData$Disease.status == 'control']]
pData_ctrl = pData[pData$Disease.status == 'control',]

# 0 DE genes / 2 / 0
atsm_signif_genes = identify_significant_genes(exprs_atsm, pData_atsm, 'Voineagu', 'Cortex.area')
atsm_signif_genes_clean = improve_df(atsm_signif_genes, map_id_gene)
xtable(atsm_signif_genes_clean)

# 29 DE genes / 33 / 44
ctrl_signif_genes = identify_significant_genes(exprs_ctrl, pData_ctrl, 'Voineagu', 'Cortex.area')
ctrl_signif_genes_clean = improve_df(ctrl_signif_genes, map_id_gene)
xtable(ctrl_signif_genes_clean)

# Bartlett's test for homocedasticity of top 20 DE ctrl genes
top_20_DE_ctrl = ctrl_signif_genes[with(ctrl_signif_genes,order(-Score.d.)),]
top_20_DE_ctrl = top_20_DE_ctrl[1:20,]

homocedasticity = c()
for(gene in top_20_DE_ctrl$Gene.Name){
  x = list(exprs_ctrl[rownames(exprs_ctrl) == gene,], exprs_atsm[rownames(exprs_atsm) == gene,])
  bt = bartlett.test(x)
  homocedasticity = c(homocedasticity, bt$p.value)
}
names(homocedasticity) = top_20_DE_ctrl$Gene.Name

library(ggplot2)
plot_df = data.frame('ILMN_ID'=names(homocedasticity), 'Homocedasticity'=homocedasticity)
plot_df = merge(plot_df, map_id_gene, by='ILMN_ID', all.x = TRUE)
plot_df$significance = plot_df$Homocedasticity>(0.05/nrow(plot_df)) # Bonferroni correction
plot_df$Gene_Symbol = factor(plot_df$Gene_Symbol, levels=rev(plot_df$Gene_Symbol))

ggplot(plot_df, aes(x=Gene_Symbol, y=Homocedasticity, fill=significance)) + geom_bar(stat='identity') + 
  geom_hline(yintercept=0.0025, color='gray') + theme_minimal() + coord_flip() + xlab('Genes') + 
  ylab('Bartlett\'s test p-value') + theme(legend.position="none") + scale_y_sqrt() + 
  ggtitle('Homocedasticity test between autism and control')

xtable(plot_df[,c('Gene_Symbol','Homocedasticity')])

######################################################################################################
######################################################################################################
# DE BY AGE USING COLANTUONI DATA

folder = 'avg'
load(paste('Data/RDatas', folder, 'LumiBatch_Colantuoni_preprocessed.RData', sep='/'))
exprs = exprs(LumiBatch_c)
pData = pData(LumiBatch_c)

age_group_signif_genes = identify_significant_genes(exprs, pData, 'Colantuoni', 'age_group')
gene_idx = substring(age_group_signif_genes$Gene.ID, 2)
age_group_signif_genes$Gene.Name = rownames(exprs[as.numeric(gene_idx),])

pData_new_age_grps = pData
pData_new_age_grps$age_group[pData$age_group %in% c('Infant','Child')] = 'Infant_Child'
pData_new_age_grps$age_group[pData$age_group %in% c('20s','30s')] = '20s_30s'
pData_new_age_grps$age_group[pData$age_group %in% c('40s','50s')] = '40s_50s'
pData_new_age_grps$age_group[pData$age_group %in% c('60s','70s')] = '60s_70s'
age_group_signif_genes=identify_significant_genes(exprs, pData_new_age_grps, 'Colantuoni', 'age_group')
age_group_signif_genes=identify_significant_genes(exprs, pData_new_age_grps, 'Colantuoni',
                                                  'age_group', 0.005)
age_group_signif_genes=identify_significant_genes(exprs, pData_new_age_grps, 'Colantuoni',
                                                  'age_group', 0.0005)

pData_new_age_grps$age_group[pData_new_age_grps$age_group %in% c('20s_30s','40s_50s')] = '20s_50s'
age_group_signif_genes=identify_significant_genes(exprs, pData_new_age_grps, 'Colantuoni',
                                                  'age_group', 0.005)
gene_idx = substring(age_group_signif_genes$Gene.ID, 2)
age_group_signif_genes$Gene.Name = rownames(exprs[as.numeric(gene_idx),])

# 1837 / 1150 / 639
pData$age = pData$`age:ch1`
age_signif_genes = identify_significant_genes(exprs, pData, 'Colantuoni', 'age')
age_signif_genes_clean = improve_df(age_signif_genes, map_id_gene)

sum(age_signif_genes$Gene.Name %in% age_group_signif_genes$Gene.Name)/nrow(age_signif_genes)
nrow(age_group_signif_genes)/ nrow(exprs)

# score_relations = merge(data.frame('Gene' = age_group_signif_genes$Gene.Name, 'age_group_rank' = 
#   rownames(age_group_signif_genes)), data.frame('Gene' = age_signif_genes$Gene.Name, 
#   'age_rank' = rownames(age_signif_genes)), by='Gene')
# score_relations$age_group_rank = as.numeric(as.character(score_relations$age_group_rank))
# score_relations$age_rank = as.numeric(as.character(score_relations$age_rank))
# 
# plot(score_relations$age_group_rank, score_relations$age_rank)
# 
# cor(score_relations$age_group_rank, score_relations$age_rank)
# 
# ggplot((score_relations), aes(x=age_rank, y=age_group_rank)) + geom_point(alpha=0.25) +
#   geom_smooth(method='lm', formula=y~x) + theme_minimal() +
#   ggtitle('Gene ranking relation using age and age group') +
#   theme(plot.title = element_text(hjust=0.5)) + xlab('Age Rank') +
#   ylab('Age Group Rank')

######################################################################################################
# Sex and PMI correction

# Sex and PMI
folder = 'avg'
load(paste('Data/RDatas', folder, 'LumiBatch_Colantuoni_preprocessed.RData', sep='/'))
exprs = exprs(LumiBatch_c)
pData = pData(LumiBatch_c)

residual_exprs = exprs
residual_exprs[] = 0

for(i in 1:nrow(residual_exprs)){
  df = data.frame('sex' = as.factor(pData$`Sex:ch1`), 'PMI' = pData$`postmortem interval (pmi):ch1`,
                  'exprs' = exprs[i,])
  lm_fit = lm(exprs ~ sex+PMI+1, df)
  residual_exprs[i,] = lm_fit$residuals
}

# 77 / 0 / 0
pData$age = pData$`age:ch1`
res_sign_genes = identify_significant_genes(residual_exprs, pData, 'Colantuoni', 'age')
res_sign_genes_clean = improve_df(res_sign_genes, map_id_gene)

xtable(res_sign_genes_clean[1:20,])

sum(res_sign_genes$Gene.Name %in% age_signif_genes$Gene.Name)/nrow(res_sign_genes)
sum(age_signif_genes$Gene.Name %in% res_sign_genes$Gene.Name)/nrow(age_signif_genes)

# Sex
residual_exprs = exprs
residual_exprs[] = 0

for(i in 1:nrow(residual_exprs)){
  df = data.frame('sex' = as.factor(pData$`Sex:ch1`), 'exprs' = exprs[i,])
  lm_fit = lm(exprs ~ sex+1, df)
  residual_exprs[i,] = lm_fit$residuals
}

# 1753 / 1161 / 675
pData$age = pData$`age:ch1`
res_sign_genes = identify_significant_genes(residual_exprs, pData, 'Colantuoni', 'age')
res_sign_genes_clean = improve_df(res_sign_genes, map_id_gene)

sum(res_sign_genes$Gene.Name %in% age_signif_genes$Gene.Name)/nrow(res_sign_genes)
sum(age_signif_genes$Gene.Name %in% res_sign_genes$Gene.Name)/nrow(age_signif_genes)

# scatter plot of scores
a = data.frame('NewScore' = res_sign_genes$Score.d., 'Gene' = res_sign_genes$Gene.Name)
b = data.frame('Pos_before' = age_signif_genes$Score.d., 'Gene' = age_signif_genes$Gene.Name)
c = merge(a,b,by='Gene')

ggplot(c, aes(x=Pos_before, y=NewScore, colour=NewScore)) + geom_bar(stat='identity') + theme_minimal()+ 
  theme(plot.title = element_text(hjust = 0.5)) + ggtitle('Scores of remaining DE genes') + 
  xlab('Gene ranking in original scoring') +  ylab('Score after correction')

######################################################################################################
# HIERARCHICAL CLUSTERING AND HEATMAP FOR 200 DE GENES

# 200 DE genes
top_200_DE_genes = age_signif_genes[with(age_signif_genes,order(-Score.d.)),]
top_200_DE_genes = top_200_DE_genes[1:200,]
exprs_top_200 = data.frame(exprs[rownames(exprs) %in% top_200_DE_genes$Gene.Name,])
pData$PMI = as.numeric(pData$`postmortem interval (pmi):ch1`)

# Color palettes
rg_palette = colorRampPalette(c('green', 'green', 'black', 'red', 'red'))(n = 100)

colormap_age = gg_color_hue(10)
names(colormap_age) = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')

colormap_sex = c(gg_color_hue(2),'#ffffff')
names(colormap_sex)= c('F', 'M', '5')

# Heatmap
annotation_df = data.frame('Disease'=pData$age_group)
cols_annotation_df = list('Disease'=colormap_age)
hm_top_annotation = HeatmapAnnotation(annotation_df, col=cols_annotation_df)

annotation_df = data.frame('Sex'=pData$`Sex:ch1`, 'PMI'=pData$PMI)
cols_annotation_df = list('Sex'=colormap_sex, 'PMI'=circlize::colorRamp2(c(1, 10),
                          c('#FF0080', '#00FF80')))
hm_bottom_annotation = HeatmapAnnotation(annotation_df, col = cols_annotation_df)

Heatmap(exprs_top_200, name='Expression levels', column_title='Samples', row_title='Genes',
        show_row_names=F, show_column_names=F, cluster_rows=T, show_row_dend=F,col=rg_palette,
        show_heatmap_legend=F,
        top_annotation=hm_top_annotation, bottom_annotation=hm_bottom_annotation,
        top_annotation_height=unit(3, 'mm'), bottom_annotation_height=unit(6, 'mm') )

remove(annotation_df, cols_annotation_df, hm_bottom_annotation, hm_top_annotation)

######################################################################################################
######################################################################################################
# COMPARE RESULTS BETWEEN PROBE COLLAPSING METHODS

load('Data/DE_by_DS.RData')

join_methods_DE_genes = function(method_mvb, method_mv, method_avg){
  method_all=merge(merge(data.frame('Gene'=method_mvb$Gene,'Score_mvb'=method_mvb$Score),
                         data.frame('Gene'=method_mv$Gene,'Score_mv'=method_mv$Score),by='Gene',all=T),
                   data.frame('Gene'=method_avg$Gene,'Score_avg'=method_avg$Score), by='Gene', all=T)
  method_all$sum_nas = apply(method_all, 1, function(x) sum(is.na(x)))
  
  method_info_df = method_all
  method_info_df = method_info_df[order(-method_info_df$Score_mvb),]
  method_info_df$Rank = seq(1:nrow(method_info_df))
  method_info_df$Avg_new_score = rowMeans(method_info_df[,c('Score_mv','Score_avg')], na.rm = TRUE)
  method_info_df$Shared_with = 'One'
  method_info_df$Shared_with[method_info_df$sum_nas==0] = 'Both'
  
  return(method_info_df)
}

cortex_info_df = join_methods_DE_genes(cortex_mvb, cortex_mv, cortex_avg)

cortex_plot_df = cortex_info_df[!is.na(cortex_info_df$Score_mvb) & cortex_info_df$sum_nas!=2,]

colors = c('#ff3399','#00cc99')
ggplot(cortex_plot_df, aes(x=Rank, y=Avg_new_score, colour=Shared_with, fill=Shared_with)) + 
  geom_bar(stat='identity') + theme_minimal() + ggtitle('Average new score of original DE genes') + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = NULL) + 
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
  xlab('Rank by score in Max Var Between method') + ylab('Avg new scores')

cor(cortex_plot_df$Score_mvb[!is.na(cortex_plot_df$Score_mv)], 
    cortex_plot_df$Score_mv[!is.na(cortex_plot_df$Score_mv)])

cor(cortex_plot_df$Score_mvb[!is.na(cortex_plot_df$Score_avg)], 
    cortex_plot_df$Score_avg[!is.na(cortex_plot_df$Score_avg)])

