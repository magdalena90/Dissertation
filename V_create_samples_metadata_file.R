setwd('~/MSc/Dissertation')


# BrainBank, chip_array
Voineagu_suppMat = read.csv('Data/Voineagu/aux/Voineagu_suppMat.csv')
Voineagu_suppMat$chip_array = paste(Voineagu_suppMat$Chip, Voineagu_suppMat$Array, sep='_')
Voineagu_suppMat = Voineagu_suppMat[,c('Brain.Bank.Case.Id','Cortex.area','Disease.status','chip_array')]
Voineagu_suppMat$Cortex.area = as.character(Voineagu_suppMat$Cortex.area)

######################################################################################################
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

######################################################################################################
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

######################################################################################################
# JOIN FRONTAL+TEMPORAL INFO WITH CEREBELLUM AND ADD METADATA:

ids_mapping_all = rbind(ids_mapping_f_t, ids_mapping_c)
ids_mapping_all = ids_mapping_all[,c('GSM','chip_array','Brain.Bank.Case.Id','Cortex.area',
                                     'Disease.status')]

samples_full_data = merge(ids_mapping_all, samples_metadata, by='Brain.Bank.Case.Id')
samples_full_data$Cortex.area = as.character(samples_full_data$Cortex.area)

######################################################################################################
# CHECK IDS WITH HEADER OF RAW DATA:

headers = read.delim('Data/Voineagu/aux/headers.csv')
headers = colnames(headers)[2:ncol(headers)]
headers = gsub('\\..*', '', headers)
headers = substring(headers, 2)
headers = data.frame(unique(headers))
colnames(headers) = 'chip_array'

if (setequal(headers$chip_array, samples_full_data$chip_array)){
  write.csv(samples_full_data, file='Data/Voineagu/samples_full_metadata.csv', row.names=FALSE)
}


remove(cortex_area, disease_status, id_split, ids_mapping, ids_mapping_all, ids_mapping_c,
       ids_mapping_f_t, samples_metadata, Voineagu_suppMat, headers)
