
setwd('~/MSc/Dissertation')
source('utils.R')

library(ComplexHeatmap)
library(circlize)

folder = 'collapsed_probes_max'
load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed_w_cerebellum.RData', sep='/'))
exprs = exprs(LumiBatch_v)
pData = pData(LumiBatch_v)

######################################################################################################
# Statistically significant genes by region

filter_cereb = pData$Cortex.area=='cerebellum'
# 1073 significant genes
cortex_signif_genes = identify_significant_genes(exprs[,!filter_cereb], pData[!filter_cereb,], 'V')

# 48 significant genes
cerebellum_signif_genes = identify_significant_genes(exprs[,filter_cereb], pData[filter_cereb,], 'V')

# 178 significant genes
cerebellum_signif_genes = identify_significant_genes(exprs[,filter_cereb],pData[filter_cereb,],'V',.25)

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
# Linear regression vs age and sex, assess significance of residuals and compare vs above
exprs_cortex = exprs[,!filter_cereb]
pData_cortex = pData[!filter_cereb,]

# By age
library(glmnet)
age_group = ordered(pData_cortex[,'age_group'], levels = c('Child','10-20','20s','30s','40s','50s'))
lambda = cv.glmnet(t(exprs_cortex), age_group, family = 'multinomial')
lm_fit = glmnet(t(exprs_cortex), age_group, family = 'multinomial', lambda = lambda$lambda.1se)


lm_data = t(exprs_cortex)
lm_label = data.frame(pData[,'age_group'])
rownames(lm_label) = rownames(pData)
colnames(lm_label) = 'age_group'
age_group_levels = c('Child','10-20','20s','30s','40s','50s')
lm_label = transform(lm_label, age_group=factor(age_group, levels=age_group_levels))
lm_data = merge(lm_data, lm_label, by='row.names')

library(nnet)
lm_fit = multinom(age_group~., lm_data)




identify_significant_genes(res_age, pData_cortex, 'Voineagu')

common = merge(cortex_signif_genes, age_residuals, by = '?')
paste0(nrow(common),' common significant genes')



