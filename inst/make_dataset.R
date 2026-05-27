library(scpdata)
library(DaparToolshed)

# IMPORT DATASET -----
dta <- dou2019_mouse()


# CREATE QFEATURES WITH DAPARTOOLSHED -----

## PREPARE DATA -----
#keep only wanted conditions
dtacond <- c("C10", "SVEC", "RAW")
colselect <- which(SummarizedExperiment::colData(dta)$SampleType %in% dtacond)

#create info for coldata
nomcoldta <- colnames(SummarizedExperiment::assay(dta[[13]][, colselect]))
nomcond <- SummarizedExperiment::colData(dta)[which(SummarizedExperiment::colData(dta)$SampleType %in% dtacond), "SampleType"]
#met tt ensemble
sampledta <- data.frame(Sample.name = nomcoldta, 
                        Condition = nomcond, 
                        Bio.Rep = seq_along(nomcoldta))
sampledta$Sample.name <- gsub(".", "_", sampledta$Sample.name, fixed = TRUE)
sampledta <- sampledta[order(sampledta$Condition), ]


#select assay
assaydta <- SummarizedExperiment::assay(dta[[13]][, colselect])
colnames(assaydta) <- gsub(".", "_", colnames(assaydta), fixed = TRUE)
rownames(assaydta) <- seq_len(nrow(assaydta))

#get rowdata
rowdta <- SummarizedExperiment::rowData(dta[[13]])
rownames(rowdta) <- seq_len(nrow(rowdta))

#create columns for metacell
metadta <- as.data.frame(replicate(ncol(assaydta), 
                                   rep("By MS/MS", nrow(assaydta)), 
                                   simplify = FALSE))
metadta <- setNames(metadta, paste0("Id_", colnames(assaydta)))

#merge together
quantdta <- cbind(assaydta, rowdta, metadta)
quantdta <- as.data.frame(quantdta)
colnames(quantdta) <- gsub(".", "_", colnames(quantdta), fixed = TRUE)


## CREATE QFEATURES -----
objt <- createQFeatures(data = quantdta,
                        file = 'dou2019_mouse',
                        sample = sampledta,
                        indQData = 1:72,
                        keyId = "sequence",
                        indexForMetacell = 87:158,
                        logData = TRUE,
                        force.na = TRUE,
                        typeDataset = "peptide",
                        parentProtId = "DatabaseAccess",
                        analysis = "Pept_Data",
                        software = "maxquant")



# PROCESS DATA THROUGH WORKFLOW -----
## SELECT SUBSET -----
obj <- objt[1:1000,]


## FILTERING -----
### Create filters -----
#create filter to remove empty lines
filter_emptyline <- FunctionFilter("qMetacellWholeLine",
                                   cmd = 'delete',
                                   pattern = 'Missing MEC')
#filter lines with too many NA is every condition
filter_manyNA <- FunctionFilter('qMetacellOnConditions',
                                cmd = 'delete',
                                mode = 'AllCond',
                                pattern = c('Missing MEC', 'Missing POV'),
                                conds = SummarizedExperiment::colData(obj)$Condition,
                                percent = TRUE,
                                th = 0.8,
                                operator = '>')
#create filter to remove contaminant
filter_contaminant <- QFeatures::VariableFilter(field = "isContaminant",
                                                value = "TRUE",
                                                condition = "==",
                                                not = TRUE)

### Apply filters -----
obj <- filterFeaturesOneSE(object = obj,
                           i = length(obj),
                           name = "Filtered",
                           filters = list(filter_emptyline, filter_manyNA, filter_contaminant))

#remove proteins with no associated peptides
X <- QFeatures::adjacencyMatrix(obj[[length(obj)]])
SummarizedExperiment::rowData(obj[[length(obj)]])$adjacencyMatrix <- X[, -which(Matrix::colSums(X) == 0)]


## NORMALIZATION -----
obj <- normalizeFunction(obj,
                         method = "MeanCentering",
                         scaling = TRUE,
                         type = "overall")


## IMPUTATION -----
obj <- wrapperPirat(data = obj,
                    adjmat = SummarizedExperiment::rowData(obj[[length(obj)]])$adjacencyMatrix,
                    extension = "base")


## AGGREGATION -----
obj <- RunAggregation(obj,
                      adjMatrix = 'adjacencyMatrix',
                      includeSharedPeptides = 'Yes_Iterative_Redistribution',
                      operator = 'Mean',
                      considerPeptides = 'allPeptides',
                      ponderation = "Global",
                      max_iter = 500
)

#apply the filter to remove empty lines created
obj <- filterFeaturesOneSE(object = obj,
                           i = length(obj),
                           name = "FilterProt",
                           filters = list(filter_emptyline))



# SAVE DATA -----
saveRDS(obj, file = "data3c24e_FNIAmini1000.qf")

