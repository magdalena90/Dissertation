setwd('~/MSc/Dissertation')

library(GEOquery)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(wesanderson)
library(varhandle)

######################################################################################################
######################################################################################################
# COLANTUONI PAPER:

gse = getGEO(filename='Data/Colantuoni/GSE30272_series_matrix.txt.gz', getGPL = FALSE)

######################################################################################################
# EXPRESSION DATA
edata = exprs(gse)
#edata = edata[rowSums(is.na(edata)) == 0,]

rand_col = sample(colnames(edata), 1) # 'GSM750017'
rand_col_mean = round(mean(edata[,rand_col]),4)
rand_col_var = round(var(edata[,rand_col]),2)

sample = ggplot(data.frame(edata),aes(edata[,rand_col])) + geom_histogram(binwidth=0.1,fill='#336699') +
  xlab(rand_col) + theme_minimal() +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0('Random sample distribution (mean=',rand_col_mean,', variance=',rand_col_var,')'))

edata_stats = data.frame(cbind(seq(1,ncol(edata)), apply(edata, 2, mean), apply(edata, 2, var)))
colnames(edata_stats) = c('n','mean','var')
mean_mean = formatC(mean(edata_stats$mean), format = 'e', digits = 2)
mean_var = formatC(var(edata_stats$mean), format = 'e', digits = 2)
var_mean = formatC(mean(edata_stats$var), format = 'e', digits = 2)
var_var = formatC(var(edata_stats$var), format = 'e', digits = 2)

mean = ggplot(edata_stats, aes(x=n, y=mean)) + geom_point(color='#6600cc',alpha=0.5) + theme_minimal() +
  ggtitle(paste0('Sample means (mean=',mean_mean,', var=',mean_var,')')) + xlab('') +
  theme(plot.title = element_text(hjust = 0.5))

var = ggplot(edata_stats, aes(x=n, y=var)) + geom_point(color='#ff0066',alpha=0.5) + theme_minimal() +
  ggtitle(paste0('Sample vars (mean=',var_mean,', var=',var_var,')')) + xlab('') +
  theme(plot.title = element_text(hjust = 0.5)) + expand_limits(y = 0)

grid.arrange(sample, mean, var, ncol=1)

######################################################################################################
# METADATA

age_to_weeks <- function(age){
  age_weeks = ceiling(age*52-1)
  weeks = 40 + age_weeks
  return(weeks)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

conf_vars = data.frame('age'=gse$`age:ch1`, 'pmi'=gse$`postmortem interval (pmi):ch1`,
                       'gender'=gse$`Sex:ch1`, 'ethnicity'=gse$`race:ch1`, 
                       'med_ex_office'=gse$`medical examiner office:ch1`, 
                       'brain_bank_src'=gse$`brain bank source:ch1`,
                       'RIN'=gse$`rna integrity number (rin):ch1`)

# Age distribution
conf_vars$age = unfactor(conf_vars$age)
conf_vars$age_group = as.factor(cut(conf_vars$age, c(-0.5, 0, 0.6, 10, 20, 30, 40, 50, 60, 70, 80),
    labels=c('Fetal','Infant','Child','10-20','20s','30s','40s','50s','60s','70s')))

age_hist = ggplot(data=conf_vars, aes(x=ceiling(age-1), fill=age_group)) +
  geom_histogram(stat='count') + guides(fill=FALSE) + xlab('Age') + theme_minimal() +
  ggtitle("Sample's age distribution") + theme(plot.title = element_text(hjust = 0.5))

fetal_hist = ggplot(data=conf_vars[conf_vars$age_group=='Fetal',], aes(x=age_to_weeks(age))) + 
  geom_histogram(stat='count', fill=gg_color_hue(10)[1]) + guides(fill=FALSE) + 
  xlab('Gestational weeks') + theme_minimal()# + coord_fixed(ratio=0.3)

infant_hist = ggplot(data=conf_vars[conf_vars$age_group=='Infant',], aes(x=ceiling(age*10)/10)) + 
  geom_histogram(stat='count', fill=gg_color_hue(10)[2]) + guides(fill=FALSE) + 
  xlab('Infant age (yrs)') + theme_minimal()# + coord_fixed(ratio=0.05)

group_age_hist = ggplot(data=conf_vars, aes(x=age_group, fill=age_group)) + 
  geom_histogram(stat='count') + guides(fill=FALSE) + xlab('Age Groups') + theme_minimal()

ggarrange(age_hist + annotation_custom(grob=ggplotGrob(fetal_hist),xmin=35,xmax=55,ymin=15,ymax=35) +
  annotation_custom(grob=ggplotGrob(infant_hist), xmin=60, xmax=Inf, ymin=15, ymax=35), group_age_hist,
  nrow=2, heights=c(3,2))

remove(age_hist, fetal_hist, infant_hist, group_age_hist, gg_color_hue, age_to_weeks)

# Confounding variables
pmi = ggplot(conf_vars, aes(x=age_group, y=unfactor(pmi), fill=age_group)) + geom_boxplot() +
  xlab('Age Group') + ylab('Post Mortem Interval') + theme_minimal() + theme(legend.position='none')

rin = ggplot(conf_vars, aes(x=age_group, y=unfactor(RIN), fill=age_group)) + geom_boxplot() +
  xlab('Age Group') + ylab('RNA Integrity Number') + theme_minimal() + theme(legend.position='none')

gender = ggplot(conf_vars, aes(x=age_group, fill=gender)) + 
  geom_histogram(stat='count', position='fill') + xlab('Age Group') + ylab('percentage') + 
  theme_minimal() + scale_fill_manual(values=c("#999999", "#ff8080", "#56B4E9")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ethnicity = ggplot(conf_vars, aes(x=age_group, fill=ethnicity)) + 
  geom_histogram(stat='count', position='fill') + xlab('Age Group') + ylab('percentage') + 
  scale_fill_brewer(palette='Dark2') + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

meo = ggplot(conf_vars, aes(x=age_group, fill=med_ex_office)) + theme_minimal() + 
  geom_histogram(stat='count', position='fill') + xlab('Age Group') + ylab('percentage') + 
  scale_fill_manual(values=wes_palette(n=length(unique(conf_vars$med_ex_office)), 
  name='FantasticFox1')) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

brain_bank = ggplot(conf_vars, aes(x=age_group, fill=brain_bank_src)) + theme_minimal() + 
  geom_histogram(stat='count', position='fill') + xlab('Age Group') + ylab('percentage') + 
  scale_fill_manual(values=wes_palette(n=length(unique(conf_vars$med_ex_office)), 
  name='BottleRocket2')) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

grid.arrange(pmi, rin, gender, ethnicity, meo, brain_bank, nrow=3)

remove(pmi, rin, gender, ethnicity, meo, brain_bank)

######################################################################################################
######################################################################################################
# VOINEAGU PAPER

gse = getGEO(filename='Data/Voineagu/GSE28521_series_matrix.txt.gz', getGPL = FALSE)

######################################################################################################
# METADATA EXPLORATORY ANALYSIS

conf_vars = data.frame('brain_region'=gse$`tissue (brain region):ch1`,
                       'disease_status'=gse$`disease status:ch1`, 'row_count'=gse$data_row_count)

ggplot(data=conf_vars, aes(x=brain_region, fill=disease_status)) +
  geom_histogram(stat='count', position='dodge') + xlab('Brain Region') + theme_minimal() + 
  ggtitle("Sample's distribution") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c('#009999','#99cc00'))

######################################################################################################
# GENE EXPRESSION

library(d3heatmap)
library(gplots)

edata = exprs(gse)
#ggplot(data = edata, aes(x = metrics, y = teams)) + geom_tile(aes(fill = performance)) # ???

# Remove rows with NAs (Removed 1,631 rows of illumina points)
clean_edata = edata[rowSums(is.na(edata)) == 0,]

rg_palette <- colorRampPalette(c('green', 'black', 'red'))(n = 100)


color.map = function(disease_status){
  if(disease_status=='autism') '#009999'
  else '#99cc00'
}

sidebarcolors <- unlist(lapply(gse$`disease status:ch1`, color.map))
lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(1.5,4)
lhei = c(1.5,4,1)

edata_sample = clean_edata[1:5000,]

dev.off()
heatmap.2(clean_edata, trace='none', scale='row', col=rg_palette, dendrogram='col',
          labRow = FALSE, labCol = FALSE, margins = c(1, 1), na.rm=TRUE, 
          density.info='none', ColSideColors=sidebarcolors)





