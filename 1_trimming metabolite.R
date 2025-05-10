#open library
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)

###########Targeted metabolite (peak area)---------------
#open dataset from Quanbrowser
peak_area <- read_excel("./rawdata/06022025_Peakarea_all_compounds_03022025_For dimentional reduction.xlsx")
peak_area

#clean the cell names and keep in cell_name$cell_name
cell_name_final #date_plate_#_well
peak_area
cell_name <- peak_area$`Component Name`
cell_name <- strsplit(cell_name, split="_")
cell_name <- as.data.frame(cell_name)
cell_name <- t(cell_name)
cell_name <- as.data.frame(cell_name)
cell_name <- cell_name[,c("V2","V7","V8","V10")]
cell_name$V2 <- gsub("20240830","240829",cell_name$V2)
cell_name$V2 <- gsub("20240919","240912",cell_name$V2)
cell_name$V7 <- gsub("P","p",cell_name$V7)
cell_name <- cell_name %>% unite("cell_name", everything(), sep="_", remove = FALSE)

#remove unnecessary columns and convert it to numeric
peak_area$`Cell diameter` <- NULL
peak_area$`Component Name` <- NULL

peak_area <- as.data.frame(lapply(peak_area, as.numeric))
peak_area[is.na(peak_area)] <- 0

#attach cell name
rownames(peak_area) <- cell_name$cell_name

#filter cells based on sequencing result
peak_area$cell_name_final <- rownames(peak_area)
peak_area$cell_name_final

peak_area <- peak_area %>% filter(cell_name_final %in% TPM) #279 cells are retained

#filter the cells without correct pictures (based on CellCelector machine)
peak_area <- peak_area %>% filter(!cell_name_final %in% picture$Cell) #279 cells are retained

#Log-transform the Area (log10(Area+1))
peak_area$cell_name_final <- NULL
peak_area_log10 <- log10(peak_area + 1)

#split the day into two
peak_area_log10_day1 <- peak_area_log10 %>% filter(str_detect(rownames(.), "240829"))
peak_area_log10_day2 <- peak_area_log10 %>% filter(str_detect(rownames(.), "240912"))

#calculate z-score in each day
peak_area_zscore_day1 <- as.data.frame(scale(peak_area_log10_day1))
peak_area_zscore_day2 <- as.data.frame(scale(peak_area_log10_day2))

#merge z-scored data into one dataframe
peak_area_zscore_each <- rbind(peak_area_zscore_day1, peak_area_zscore_day2)

#calculate z-score across all data
peak_area_zscore_whole <- as.data.frame(scale(peak_area_log10))

#store calculated results as csv
write.csv(peak_area_log10, file = "./step1_output/peak_area_log10.csv")
write.csv(peak_area_log10_day1, file = "./step1_output/peak_area_log10_day1.csv")
write.csv(peak_area_log10_day2, file = "./step1_output/peak_area_log10_day2.csv")
write.csv(peak_area_zscore_day1, file = "./step1_output/peak_area_zscore_day1.csv")
write.csv(peak_area_zscore_day2, file = "./step1_output/peak_area_zscore_day2.csv")
write.csv(peak_area_zscore_each, file = "./step1_output/peak_area_zscore_each.csv")
write.csv(peak_area_zscore_whole, file = "./step1_output/peak_area_zscore_whole.csv")

###########Targeted metabolite (concentration)---------------
#open dataset
concentration <- read_excel("./rawdata/06022025_Concentration_all_compounds_03022025.xlsx")
concentration

#clean the cell names
cell_name_final

cell_name <- concentration$`Component Name`
cell_name <- strsplit(cell_name, split="_")
cell_name <- as.data.frame(cell_name)
cell_name <- t(cell_name)
cell_name <- as.data.frame(cell_name)
cell_name <- cell_name[,c("V2","V7","V8","V10")]
cell_name$V2 <- gsub("20240830","240829",cell_name$V2)
cell_name$V2 <- gsub("20240919","240912",cell_name$V2)
cell_name$V7 <- gsub("P","p",cell_name$V7)
cell_name <- cell_name %>% unite("cell_name", everything(), sep="_", remove = FALSE)


#store manual cell type annotation (by targeted metabolite profiles)
cell_type_manual <- concentration[,c("Cell type annotation", "Component Name")]
cell_type_manual <- cbind(cell_type_manual, cell_name$cell_name)


#remove unnecessary columns and convert it to numeric
concentration$...1 <- NULL
concentration$`Cell type annotation` <- NULL
concentration$`Component Name` <- NULL

concentration <- as.data.frame(lapply(concentration, as.numeric))
concentration[is.na(concentration)] <- 0

#attach cell name
rownames(concentration) <- cell_name$cell_name

#filter cells based on sequencing result
concentration$cell_name_final <- rownames(concentration)
concentration <- concentration %>% filter(cell_name_final %in% TPM) #280 cells are retained

#filter the cells without correct pictures (based on CellCelector machine)
concentration <- concentration %>% filter(!cell_name_final %in% picture$Cell) #279 cells are retained

#drop cell name
concentration$cell_name_final <- NULL

#log-transform the concentration (log10(conc+1))
concentration_log10 <- log10(concentration + 1)

#split the day into two
concentration_log10_day1 <- concentration_log10 %>% filter(str_detect(rownames(.), "240829"))
concentration_log10_day2 <- concentration_log10 %>% filter(str_detect(rownames(.), "240912"))

#calculate z-score in each day
concentration_zscore_day1 <- as.data.frame(scale(concentration_log10_day1))
concentration_zscore_day2 <- as.data.frame(scale(concentration_log10_day2))

#merze z-scored data into one dataframe
concentration_zscore_each <- rbind(concentration_zscore_day1, concentration_zscore_day2)

#calculate z-score across all data
concentration_zscore_whole <- as.data.frame(scale(concentration_log10))

#store calculated results as csv
write.csv(cell_type_manual, file = "./step1_output/cell_type_manual_by_Anh.csv")

write.csv(concentration, file = "./step1_output/concentration.csv")
write.csv(concentration_log10, file = "./step1_output/concentration_log10.csv")
write.csv(concentration_log10_day1, file = "./step1_output/concentration_log10_day1.csv")
write.csv(concentration_log10_day2, file = "./step1_output/concentration_log10_day2.csv")
write.csv(concentration_zscore_day1, file = "./step1_output/concentration_zscore_day1.csv")
write.csv(concentration_zscore_day2, file = "./step1_output/concentration_zscore_day2.csv")
write.csv(concentration_zscore_each, file = "./step1_output/concentration_zscore_each.csv")
write.csv(concentration_zscore_whole, file = "./step1_output/concentration_zscore_whole.csv")


###########Make Seurat object from targeted metabolite results (peak area) -------------------
## 1) Transpose the dataframes ----------
peak_area_t <- t(peak_area)
rownames(peak_area_t) <- gsub("_","-",rownames(peak_area_t)) #exchange the feature names (_ to -)

peak_area_log10_t <- t(peak_area_log10)
rownames(peak_area_log10_t) <- gsub("_","-",rownames(peak_area_log10_t)) #exchange the feature names (_ to -)

zscore_t <- t(peak_area_zscore_whole)
rownames(zscore_t) <- gsub("_","-",rownames(zscore_t)) #exchange the feature names (_ to -)

## 2) Make Seurat object (Assay) ---------
#create seurat assay
assay_targeted_metabolite <- CreateAssayObject(counts = peak_area_t)

#make log10(peak_area+1) as a dgCMatrix and append in data slot
peak_area_log10_t <- as(peak_area_log10_t, "sparseMatrix")
assay_targeted_metabolite@data <- peak_area_log10_t

#append zscore in scale.data slot
assay_targeted_metabolite@scale.data <- zscore_t

#save as rds
saveRDS(assay_targeted_metabolite, file = "./step1_output/targeted_metabolite_peakarea_Mar2025.rds")

###########Make Seurat object from targeted metabolite results (concentration) -------------------
## 1) Transpose the dataframes ----------
concentration_t <- t(concentration)
rownames(concentration_t) <- gsub("_","-",rownames(concentration_t)) #exchange the feature names (_ to -)

concentration_log10_t <- t(concentration_log10)
rownames(concentration_log10_t) <- gsub("_","-",rownames(concentration_log10_t)) #exchange the feature names (_ to -)

zscore_t <- t(concentration_zscore_whole)
rownames(zscore_t) <- gsub("_","-",rownames(zscore_t)) #exchange the feature names (_ to -)

## 2) Make Seurat object (Assay) ---------
#create seurat assay
assay_targeted_metabolite <- CreateAssayObject(counts = concentration_t)

#make log10(concentration+1) as a dgCMatrix and append in data slot
concentration_log10_t <- as(concentration_log10_t, "sparseMatrix")
assay_targeted_metabolite@data <- concentration_log10_t

#append zscore in scale.data slot
assay_targeted_metabolite@scale.data <- zscore_t

#save as rds
saveRDS(assay_targeted_metabolite, file = "./step1_output/targeted_metabolite_concentration_Mar2025.rds")
