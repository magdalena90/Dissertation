
setwd('~/MSc/Dissertation')
source('utils.R')

library(ComplexHeatmap)
library(circlize)

folder = 'collapsed_probes_max'
load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed_w_cerebellum.RData', sep='/'))
exprs = exprs(LumiBatch_v)
pData = pData(LumiBatch_v)

######################################################################################################
# STATISTICALLY SIGNIFICANT GENES BY REGION

filter_cereb = pData$Cortex.area=='cerebellum'
# 1073 significant genes
cortex_signif_genes = identify_significant_genes(exprs[,!filter_cereb], pData[!filter_cereb,], 
                                                 'Voineagu', 'Disease.status')

# 48 significant genes
cerebellum_signif_genes = identify_significant_genes(exprs[,filter_cereb], pData[filter_cereb,], 
                                                     'Voineagu', 'Disease.status')

# 178 significant genes
cerebellum_signif_genes = identify_significant_genes(exprs[,filter_cereb],pData[filter_cereb,],
                                                     'Voineagu', 'Disease.status', 0.25)

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
rg_palette = colorRampPalette(c('green', 'green', 'black', 'red', 'red'))(n = 100)

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
                          'Sex'=colormap_sex, 'PMI'=circlize::colorRamp2(c(200, 600),
                           c('#FF0080', '#00FF80')))
hm_bottom_annotation = HeatmapAnnotation(annotation_df, col = cols_annotation_df)

Heatmap(exprs_top_200, name='Expression levels', column_title='Samples', row_title='Genes',
        show_row_names=F, show_column_names=F, cluster_rows=F, col=rg_palette, show_heatmap_legend=F,
        top_annotation=hm_top_annotation, bottom_annotation=hm_bottom_annotation,
        top_annotation_height=unit(3, 'mm'), bottom_annotation_height=unit(12, 'mm') )

remove(annotation_df, cols_annotation_df, hm_bottom_annotation, hm_top_annotation)
######################################################################################################
# LINEAR REGRESSION VS AGE AN SEX TO ASSESS SIGNIFICANCE IN RESIDUALS

exprs_cortex = exprs[,!filter_cereb]
pData_cortex = pData[!filter_cereb,]

residual_exprs = exprs_cortex
residual_exprs[] = 0

for(i in 1:nrow(residual_exprs)){
  df = data.frame('age' = pData_cortex$AGE, 'sex' = as.factor(pData_cortex$SEX), 
                  'exprs' = exprs_cortex[i,])
  lm_fit = lm(exprs ~ age+sex+1, df)
  residual_exprs[i,] = lm_fit$residuals
}

res_sign_genes = identify_significant_genes(residual_exprs, pData_cortex, 'Voineagu', 'Disease.status')

sum(res_sign_genes$Gene.Name %in% cortex_signif_genes$Gene.Name)/nrow(res_sign_genes)

######################################################################################################
# Age and sex counts by disease status

age = pData_cortex[,c('age_group','Disease.status')]
age$Count = 1
age = aggregate(Count ~ age_group + Disease.status, age, sum)
colnames(age) = c('Age','Disease.Status','Count')
ordered_levels = c('Child','10-20','20s','30s','40s','50s')
age = transform(age, Age=factor(Age, levels = ordered_levels))

library(ggplot2)
ggplot(age, aes(x=Age, y=Count, fill=Disease.Status)) + xlab('Age group') + theme_minimal() +
  geom_bar(stat='identity', position=position_dodge()) + 
  scale_fill_manual(values=c('#009999', '#99cc00')) +
  scale_y_continuous(breaks=seq(1,max(age$Count),2)) +
  geom_text(aes(label=Count), vjust=1.6, color='white', position = position_dodge(0.9), size=3.5)

sex = pData_cortex[,c('SEX','Disease.status')]

sex$Count = 1
sex = aggregate(Count ~ SEX + Disease.status, sex, sum)

ggplot(sex, aes(x=SEX, y=Count, fill=Disease.status)) + xlab('Gender') + theme_minimal() +
  geom_bar(stat='identity', position=position_dodge()) + 
  scale_y_continuous(breaks=seq(1,max(sex$Count),5)) +
  geom_text(aes(label=Count), vjust=1.6, color='white', position = position_dodge(0.9), size=3.5)

remove(age, sex)  
######################################################################################################
######################################################################################################
# DE BY BRAIN REGION FOR AUTISM AND CONTROL SAMPLES SEPARATELY

load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed.RData', sep='/'))
exprs = exprs(LumiBatch_v)
pData = pData(LumiBatch_v)

exprs_atsm = exprs[,colnames(exprs) %in% rownames(pData)[pData$Disease.status == 'autism']]
pData_atsm = pData[pData$Disease.status == 'autism',]
exprs_ctrl = exprs[,colnames(exprs) %in% rownames(pData)[pData$Disease.status == 'control']]
pData_ctrl = pData[pData$Disease.status == 'control',]

# 0 DE genes
atsm_signif_genes = identify_significant_genes(exprs_atsm, pData_atsm, 'Voineagu', 'Cortex.area')

# 87 DE genes
ctrl_signif_genes = identify_significant_genes(exprs_ctrl, pData_ctrl, 'Voineagu', 'Cortex.area')

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
map_id_gene = data.frame(map_HEEBO_ILMN_ids()[1])
map_id_gene = unique(map_id_gene[,c('ILMN_ID','Gene_Symbol')])

plot_df = data.frame('ILMN_ID'=names(homocedasticity), 'Homocedasticity'=homocedasticity)
plot_df = merge(plot_df, map_id_gene, by='ILMN_ID', all.x = TRUE)
plot_df$significance = plot_df$Homocedasticity>(0.05/nrow(plot_df)) # Bonferroni correction
plot_df$Gene_Symbol = factor(plot_df$Gene_Symbol, levels=rev(plot_df$Gene_Symbol))

ggplot(plot_df, aes(x=Gene_Symbol, y=Homocedasticity, fill=significance)) + geom_bar(stat='identity') + 
  geom_hline(yintercept=0.05, color='gray') + theme_minimal() + coord_flip() + xlab('Genes') + 
  ylab('Bartlett\'s test p-value') + theme(legend.position="none")

xtable(plot_df[,c('Gene_Symbol','Homocedasticity')])

######################################################################################################
######################################################################################################
# DE BY AGE USING COLANTUONI DATA

load(paste('Data/RDatas', folder, 'LumiBatch_Colantuoni_preprocessed.RData', sep='/'))
exprs = exprs(LumiBatch_c)
pData = pData(LumiBatch_c)

age_group_signif_genes = identify_significant_genes(exprs, pData, 'Colantuoni', 'age_group')

pData_new_age_grps = pData
pData_new_age_grps$age_group[pData$age_group %in% c('Infant','Child')] = 'Infant_Child'
pData_new_age_grps$age_group[pData$age_group %in% c('20s','30s')] = '20s_30s'
pData_new_age_grps$age_group[pData$age_group %in% c('40s','50s')] = '40s_50s'
pData_new_age_grps$age_group[pData$age_group %in% c('60s','70s')] = '60s_70s'

pData_new_age_grps$age_group[pData_new_age_grps$age_group %in% c('20s_30s','40s_50s')] = '20s_50s'
age_group_signif_genes=identify_significant_genes(exprs, pData_new_age_grps, 'Colantuoni',
                                                  'age_group', 0.0005)
gene_idx = substring(age_group_signif_genes$Gene.ID, 2)
age_group_signif_genes$Gene.Name = rownames(exprs[as.numeric(gene_idx),])

pData$age = pData$`age:ch1`
age_signif_genes = identify_significant_genes(exprs, pData, 'Colantuoni', 'age')

sum(age_signif_genes$Gene.Name %in% age_group_signif_genes$Gene.Name)/nrow(age_signif_genes)
nrow(age_group_signif_genes)/ nrow(exprs)

score_relations = merge(data.frame('Gene' = age_group_signif_genes$Gene.Name, 'age_group_rank' = 
  rownames(age_group_signif_genes)), data.frame('Gene' = age_signif_genes$Gene.Name, 
  'age_rank' = rownames(age_signif_genes)), by='Gene')
score_relations$age_group_rank = as.numeric(as.character(score_relations$age_group_rank))
score_relations$age_rank = as.numeric(as.character(score_relations$age_rank))

plot(score_relations$age_group_rank, score_relations$age_rank)

cor(score_relations$age_group_rank, score_relations$age_rank)

ggplot((score_relations), aes(x=age_rank, y=age_group_rank)) + geom_point(alpha=0.25) +
  geom_smooth(method='lm', formula=y~x) + theme_minimal() +
  ggtitle('Gene ranking relation using age and age group') +
  theme(plot.title = element_text(hjust=0.5)) + xlab('Age Rank') +
  ylab('Age Group Rank')

######################################################################################################
# Sex correction

residual_exprs = exprs
residual_exprs[] = 0

for(i in 1:nrow(residual_exprs)){
  df = data.frame('sex' = as.factor(pData$`Sex:ch1`), 'exprs' = exprs[i,])
  lm_fit = lm(exprs ~ sex+1, df)
  residual_exprs[i,] = lm_fit$residuals
}

res_sign_genes = identify_significant_genes(residual_exprs, pData, 'Colantuoni', 'age')

sum(res_sign_genes$Gene.Name %in% age_signif_genes$Gene.Name)/nrow(res_sign_genes)

sum(age_signif_genes$Gene.Name %in% res_sign_genes$Gene.Name)/nrow(age_signif_genes)

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
        show_row_names=F, show_column_names=F, cluster_rows=F, col=rg_palette, show_heatmap_legend=F,
        top_annotation=hm_top_annotation, bottom_annotation=hm_bottom_annotation,
        top_annotation_height=unit(3, 'mm'), bottom_annotation_height=unit(6, 'mm') )

remove(annotation_df, cols_annotation_df, hm_bottom_annotation, hm_top_annotation)