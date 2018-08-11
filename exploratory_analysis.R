setwd('~/MSc/Dissertation')
source('utils.R')

library(GEOquery)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(varhandle)
library(dplyr)

######################################################################################################
######################################################################################################
# LOAD DATA
folder = 'max_var_between' # avg, max_var, max_var_between

load(paste('Data/RDatas', folder, 'LumiBatch_Colantuoni_preprocessed.RData', sep='/'))
exprs_c = exprs(LumiBatch_c)
pData_c = pData(LumiBatch_c)

load(paste('Data/RDatas', folder, 'LumiBatch_Voineagu_preprocessed.RData', sep='/'))
exprs_v = exprs(LumiBatch_v)
pData_v = pData(LumiBatch_v)

######################################################################################################
######################################################################################################
# SCATTER PLOTS BETWEEN SOURCES

# Non fetal samples vs non autism samples
exprs_c_non_fetal = exprs_c[,colnames(exprs_c) %in% rownames(pData_c[pData_c$age_group != 'Fetal',])]
exprs_v_ctrl = exprs_v[,colnames(exprs_v) %in% rownames(pData_v[pData_v$Disease.status=='control',])]
scatterplot_c_vs_v(exprs_c_non_fetal, exprs_v_control, 'Expression Comparisons','sp_all')

# By age
common_age_groups = c('10-20', '20s', '30s', '40s', '50s') # All Child samples in V have autism
age_group_colours = gg_color_hue(10)

for(i in 1:length(common_age_groups)){
  c = exprs_c[,colnames(exprs_c) %in% rownames(pData_c[pData_c$age_group == common_age_groups[i],])]
  v = exprs_v_ctrl[,colnames(exprs_v_ctrl) %in% rownames(pData_v[pData_v$Disease.status == 'control' & 
                                                         pData_v$age_group == common_age_groups[i],])]
  f_name = paste0('sp_', common_age_groups[i],'.png')
  scatterplot_c_vs_v(c, v, paste0('Age group ', common_age_groups[i]), f_name, age_group_colours[i+3])
}

scatterplot_c_vs_v = function(exprs_c, exprs_v, title, file_name, colour = 'black'){
  
  # Calculate mean expression by probe
  mean_exprs_c = data.frame('ID'=rownames(exprs_c), 'mean_c'=apply(exprs_c, 1, mean))
  mean_exprs_v = data.frame('ID'=rownames(exprs_v), 'mean_v'=apply(exprs_v, 1, mean))
  
  # Create dataframe for plot
  mean_exprs = merge(mean_exprs_c, mean_exprs_v, by='ID')
  rownames(mean_exprs) = mean_exprs$ID
  mean_exprs$ID = NULL
  
  # Correlation between datasets
  cor_c_v = cor(mean_exprs$mean_c, mean_exprs$mean_v)
  title = paste0(title, ' (cor = ', round(cor_c_v,4), ')')
  
  # Plot
  if(colour=='black') { colour_line = '#005c99' } else { colour_line = '#595959'}
  scatter_plot = ggplot((mean_exprs), aes(x = mean_c, y = mean_v)) + theme_minimal() + 
    geom_point(alpha = 0.25, color = colour) + geom_smooth(method = 'lm', formula = y ~ x, 
      color = colour_line) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + 
    xlab('Colantuoni\'s data') + ylab('Voineagu\'s data') + coord_fixed()
  
  print(scatter_plot)
  ggsave(paste0('Plots/ExploratoryAnalysis/',file_name), plot = scatter_plot, device = 'png',
         width=3.66, height=3.66, units = 'in')
}

remove(exprs_c_non_fetal, exprs_v_ctrl, common_age_groups, age_group_colours, c, v, i, f_name)
######################################################################################################
######################################################################################################
# COLANTUONI:

conf_vars = data.frame('age'=pData_c$`age:ch1`, 'pmi'=pData_c$`postmortem interval (pmi):ch1`,
                       'gender'=pData_c$`Sex:ch1`,'RIN'=pData_c$`rna integrity number (rin):ch1`,
                       'age_group'=pData_c$age_group)
conf_vars$age = unfactor(conf_vars$age)
ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
conf_vars = transform(conf_vars, age_group=factor(age_group, levels = ordered_levels))

age_to_weeks <- function(age){
  age_weeks = ceiling(age*52-1)
  weeks = 40 + age_weeks
  return(weeks)
}

# Age counts
age_hist = ggplot(data=conf_vars, aes(x=ceiling(age-1), fill=age_group)) +
  geom_bar(stat='count') + guides(fill=FALSE) + xlab('Age') + theme_minimal() +
  ggtitle("Sample's age distribution") + theme(plot.title = element_text(hjust = 0.5))

fetal_hist = ggplot(data=conf_vars[conf_vars$age_group=='Fetal',], aes(x=age_to_weeks(age))) + 
  geom_bar(stat='count', fill=gg_color_hue(10)[1]) + guides(fill=FALSE) + 
  xlab('Gestational weeks') + theme_minimal()# + coord_fixed(ratio=0.3)

infant_hist = ggplot(data=conf_vars[conf_vars$age_group=='Infant',], aes(x=ceiling(age*10)/10)) + 
  geom_bar(stat='count', fill=gg_color_hue(10)[2]) + guides(fill=FALSE) + 
  xlab('Infant age (yrs)') + theme_minimal()# + coord_fixed(ratio=0.05)

group_age_hist = ggplot(data=conf_vars, aes(x=age_group, fill=age_group)) + geom_bar(stat='count') + 
  geom_text(stat='count', aes(label=..count..), color='white', vjust=1.2) + guides(fill=FALSE) + 
  xlab('Age Groups') + theme_minimal()

ggarrange(age_hist + annotation_custom(grob=ggplotGrob(fetal_hist),xmin=30,xmax=55,ymin=15,ymax=35) +
  annotation_custom(grob=ggplotGrob(infant_hist), xmin=60, xmax=Inf, ymin=15, ymax=35), group_age_hist,
  nrow=2, heights=c(3,2))

remove(age_hist, fetal_hist, infant_hist, group_age_hist, age_to_weeks)

# Boxplots of Expression by Age Group
exprs_c_grouped = data.frame(t(aggregate(t(exprs_c), list(pData_c$age_group), mean)))
colnames(exprs_c_grouped) = c('10-20','20s','30s','40s','50s','60s','70s','Child','Fetal','Infant')
exprs_c_grouped = exprs_c_grouped[-1,]
exprs_c_grouped = apply(exprs_c_grouped, 2, function(x) as.numeric(as.character(x)))
exprs_c_grouped = melt(exprs_c_grouped)
exprs_c_grouped = exprs_c_grouped[,-1]
colnames(exprs_c_grouped) = c('age_group','expression')
ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
exprs_c_grouped = transform(exprs_c_grouped, age_group=factor(age_group, levels = ordered_levels))

ggplot(exprs_c_grouped, aes(x=age_group, y=expression, fill=age_group)) + geom_boxplot() +
  xlab('Age Group') + ylab('Expression levels') + theme_minimal() + ggtitle('Expression Levels by Age')+
  theme(legend.position='none', plot.title = element_text(hjust = 0.5))

remove(exprs_c_grouped, ordered_levels)

# Confounding variables
pmi = ggplot(conf_vars, aes(x=age_group, y=unfactor(pmi), fill=age_group)) + geom_boxplot() +
  ggtitle('PMI distribution by Age Group') + theme_minimal() + ylab('Post Mortem Interval') +
  theme(plot.title = element_text(hjust = 0.5), legend.position='none', axis.title.x=element_blank())

gender_df = conf_vars[,c('age_group','gender')]
gender_df = gender_df %>% group_by(age_group, gender) %>% summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% mutate(perc = paste0(round((100*freq)),'%'))

gender = ggplot(gender_df, aes(x=age_group, y=freq, fill=gender)) +  ylab('percentage') + 
  geom_bar(stat='identity') + scale_fill_manual(values=c('#ff8080', '#56B4E9')) +
  geom_text(aes(label = perc, y=freq), position = 'stack', color = 'white', vjust = 1.5) +
  ggtitle('Gender distribution by Age Group') + theme_minimal() + theme(plot.title = 
    element_text(hjust=.5), legend.position='bottom', axis.title.x=element_blank())

grid.arrange(pmi, gender, nrow=2, heights=c(5,4))

remove(pmi, gender, conf_vars, ordered_levels)
######################################################################################################
######################################################################################################
# VOINEAGU PAPER

conf_vars = data.frame('brain_region' = pData_v$Cortex.area, 'disease_status' = pData_v$Disease.status, 
                       'age_group' = pData_v$age_group, 'age' = pData_v$AGE, 'gender' = pData_v$SEX,
                       'pmi' = pData_v$PMI)
ordered_levels = c('Child','10-20','20s','30s','40s','50s')
conf_vars = transform(conf_vars, age_group = factor(age_group, levels = ordered_levels))

# Counts by age group
group_age = ggplot(data=conf_vars, aes(x=age_group, fill=age_group)) + guides(fill=FALSE) + 
  geom_bar(stat='count') + scale_fill_manual(values=gg_color_hue(10)[3:8]) + 
  geom_text(stat='count', aes(label=..count..), color='white', vjust=1.5) + theme_minimal() + 
  theme(plot.title = element_text(hjust=.5), legend.position='bottom', axis.title.x=element_blank()) + 
  ggtitle('Sample by Age')

brain_region_df = conf_vars[,c('age_group','brain_region')]
brain_region_df = brain_region_df %>% group_by(age_group, brain_region) %>% summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% mutate(perc = paste0(round((100*freq)),'%'))

brain_region = ggplot(brain_region_df, aes(x=age_group, y=freq, fill=brain_region)) + 
  ylab('percentage') + geom_bar(stat='identity') + theme_minimal() + 
  geom_text(aes(label = perc, y=freq), position = 'stack', color = 'white', vjust = 1.2) +
  scale_fill_manual(values=c('#7324A6', '#f43e6c', '#FFDB00')) + theme(plot.title = element_text(
    hjust=.5), legend.position='bottom', axis.title.x=element_blank()) + 
  ggtitle('Brain region distribution by Age')

disease_status_df = conf_vars[,c('age_group','disease_status')]
disease_status_df = disease_status_df %>% group_by(age_group, disease_status) %>% summarize(n=n()) %>% 
  mutate(freq = n/sum(n)) %>% mutate(perc = paste0(round((100*freq)),'%'))

disease_status = ggplot(disease_status_df, aes(x=age_group, y=freq, fill=disease_status)) + 
  ylab('percentage') + geom_bar(stat='identity') + theme_minimal() + 
  geom_text(aes(label = perc, y=freq), position = 'stack', color = 'white', vjust = 1.2) +
  scale_fill_manual(values=c('#009999','#99cc00')) + theme(plot.title = element_text(hjust=.5), 
    legend.position='bottom', axis.title.x=element_blank()) + 
  ggtitle('Disease Status distribution by Age')

grid.arrange(group_age, brain_region, disease_status, nrow=3, heights=c(3,4,4))

remove(group_age, brain_region, disease_status, brain_region_df, disease_status_df)

# Boxplots of Expression by Age Group and by Disease Status
exprs_v_grouped = data.frame(t(aggregate(t(exprs_v), list(pData_v$age_group), mean)))
colnames(exprs_v_grouped) = c('10-20','20s','30s','40s','50s','Child')
exprs_v_grouped = exprs_v_grouped[-1,]
exprs_v_grouped = apply(exprs_v_grouped, 2, function(x) as.numeric(as.character(x)))
exprs_v_grouped = melt(exprs_v_grouped)
exprs_v_grouped = exprs_v_grouped[,-1]
colnames(exprs_v_grouped) = c('age_group','expression')
ordered_levels = c('Child','10-20','20s','30s','40s','50s')
exprs_v_grouped = transform(exprs_v_grouped, age_group=factor(age_group, levels = ordered_levels))

age_group = ggplot(exprs_v_grouped, aes(x=age_group, y=expression, fill=age_group)) + geom_boxplot() +
  xlab('Age Group') + ylab('Expression levels') + theme_minimal() +
  theme(legend.position='none', plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values=gg_color_hue(10)[3:8]) 

exprs_v_grouped_ds = data.frame(t(aggregate(t(exprs_v), list(pData_v$Disease.status), mean)))
colnames(exprs_v_grouped_ds) = c('autism','control')
exprs_v_grouped_ds = exprs_v_grouped_ds[-1,]
exprs_v_grouped_ds = apply(exprs_v_grouped_ds, 2, function(x) as.numeric(as.character(x)))
exprs_v_grouped_ds = melt(exprs_v_grouped_ds)
exprs_v_grouped_ds = exprs_v_grouped_ds[,-1]
colnames(exprs_v_grouped_ds) = c('disease_status','expression')

disease_status = ggplot(exprs_v_grouped_ds, aes(x=disease_status, y=expression, fill=disease_status)) + 
  geom_boxplot() + xlab('disease_status') + theme_minimal() +
  theme(legend.position='none', plot.title = element_text(hjust = 0.5), axis.title.y=element_blank()) + 
  scale_fill_manual(values=c('#009999','#99cc00'))

grid.arrange(age_group, disease_status, nrow=1, widths=c(3,1), 
             top='Expression Level by Age Group and Disease Status')

remove(exprs_v_grouped, exprs_v_grouped_ds, ordered_levels)

# Confounding vars
pmi_data = conf_vars[conf_vars$pmi != 'na',]
pmi_data$pmi = as.numeric(gsub(',', '', pmi_data$pmi))
pmi = ggplot(pmi_data, aes(x=disease_status, y=pmi, fill=disease_status)) + geom_boxplot() +
  theme_minimal() + ylab('Post Mortem Interval') +
  theme(plot.title = element_text(hjust = 0.5), legend.position='none', axis.title.x=element_blank()) +
  scale_fill_manual(values=c('#009999','#99cc00'))

gender_df = conf_vars[,c('disease_status','gender')]
gender_df = gender_df %>% group_by(disease_status, gender) %>% summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% mutate(perc = paste0(round((100*freq)),'%'))

gender = ggplot(gender_df, aes(x=disease_status, y=freq, fill=gender)) +  ylab('percentage') + 
  geom_bar(stat='identity') + scale_fill_manual(values=c('#ff8080', '#56B4E9')) +
  geom_text(aes(label = perc, y=freq), position = 'stack', color = 'white', vjust = 1.2) +
  theme_minimal() + theme(plot.title = 
    element_text(hjust=.5), legend.position='right', axis.title.x=element_blank())

grid.arrange(pmi, gender, nrow=1, widths = c(3,4), top='Confounding variables for Disease Status')

remove(pmi_data, gender_df, pmi, gender)

# Counts by Brain region and Disease status
ggplot(data=conf_vars, aes(x=brain_region, fill=disease_status)) +
  geom_bar(stat='count', position='dodge') + xlab('Brain Region') + theme_minimal() + 
  geom_text(stat='count', aes(label=..count..), position = position_dodge(width = 1), vjust=2, 
    color='white') + scale_fill_manual(values=c('#009999','#99cc00')) +
  ggtitle('Sample\'s distribution by Disease Status and Brain Region') + 
  theme(plot.title = element_text(hjust = 0.5))

######################################################################################################
######################################################################################################
# VISUALISATIONS USING DIMENSIONALITY REDUCTION TECHNIQUES

pData = data.frame('ID' = rownames(pData_c), 'age_group' = pData_c$age_group)

# By Age Group
pca = perform_pca(exprs_c, pData, 'PCA plot',1,1)
mds = perform_mds(exprs_c, pData, 'MDS plot')
lda = perform_lda(exprs_c, pData, 'LDA plot')
pls_lda = perform_lda(exprs_c, pData, 'PLS.LDA plot', pls = TRUE)
#tsne = perform_tsne(exprs_c, pData, 'TSNE plot') # Run in terminal, it dies on RStudio!

# Load tsne output ran on terminal and plot
tsne = read.csv('~/MSc/Dissertation/tsne_50.csv')
ordered_levels = c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')
tsne = transform(tsne, age_group=factor(age_group, levels = ordered_levels))

ggplot(tsne, aes(x=x, y=y, colour=age_group)) + geom_point() + ggtitle('TSNE with perplexity=50') +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) + theme_minimal()

# By Disease Status
folder = 'max_var_btwn'
pData = data.frame('ID' = rownames(pData_v), 'disease_status' = pData_v$DISORDER)
pData$disease_status = unfactor(pData$disease_status)
pData$disease_status[pData$disease_status=='Autism'] = 'autism'
pData$disease_status[pData$disease_status=='No Known Disorder'] = 'control'
pData$disease_status[pData$disease_status=='Autism,Chromosome 15q Duplication'] = 'autism (ch 15q dup)'

pca = perform_pca(exprs_v, pData, 'PCA plot', 'DiseaseStatus', folder)
mds = perform_mds(exprs_v, pData, 'MDS plot', 'DiseaseStatus', folder)
lda = perform_lda(exprs_v, pData, 'LDA plot', 'DiseaseStatus', folder)
pls_lda = perform_lda(exprs_v, pData, 'PLS.LDA plot', 'DiseaseStatus', folder, pls = TRUE)

# By Brain Region
pData = data.frame('ID' = rownames(pData_v), 'brain_region' = pData_v$Cortex.area)
pca = perform_pca(exprs_v, pData, 'PCA plot', 'CortexArea', folder)
mds = perform_mds(exprs_v, pData, 'MDS plot', 'CortexArea', folder)
lda = perform_lda(exprs_v, pData, 'LDA plot', 'CortexArea', folder)
pls_lda = perform_lda(exprs_v, pData, 'PLS.LDA plot', 'CortexArea', folder, pls = TRUE)

