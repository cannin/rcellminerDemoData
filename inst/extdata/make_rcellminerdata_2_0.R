library(rcellminer)
library(stringr)

#--------------------------------------------------------------------------------------------------
# HELPER FUNCTIONS
#--------------------------------------------------------------------------------------------------

validateDataTab <- function(dataTab, keyCol = "Gene name", featureDataColNums = 1:9){
  stopifnot(keyCol %in% colnames(dataTab))
  
  # Remove any whitespace in feature data columns.
  for(j in featureDataColNums){
    if(is.character(dataTab[, j])){
      dataTab[, j] <- stringr::str_trim(dataTab[, j])
    }
  }
  
  # Make sure expected numeric data is numeric (symbols to be read in as NAs are
  # properly handled).
  stopifnot(all(c(lapply(dataTab[, -featureDataColNums], is.numeric), recursive = TRUE)))
  
  # Set row names of data table to names in key column after checks.
  stopifnot(all(!is.na(dataTab[, keyCol])))
  stopifnot(all(dataTab[, keyCol] != ""))
  stopifnot(all(dataTab[, keyCol] != "1-Mar"))    # No Excel conversion to dates.
  stopifnot(all(!duplicated(dataTab[, keyCol])))
  rownames(dataTab) <- dataTab[, keyCol]
  
  return(dataTab)
}

#--------------------------------------------------------------------------------------------------
# LOAD DATA: MRNA EXPRESSION.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [RNA: 5 Platform Gene Transcript, select: Average z score]
# Processing of data file (RNA__5_Platform_Gene_Transcript_Average_z_scores.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

#----[z score data]------------------------------------------------------------
filePath <- "inst/extdata/cellminer_2_0/RNA__5_Platform_Gene_Transcript_Average_z_scores.txt"
expTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
expTabOrig <- validateDataTab(expTabOrig, keyCol = "Gene name",
                              featureDataColNums = featureDataCols)

expData <- ExpressionSet(as.matrix(expTabOrig[, -featureDataCols]))
featureData(expData) <- new("AnnotatedDataFrame", data=expTabOrig[, featureDataCols])

###################################################################################################
# Set CellMiner NCI-60 Cell Line Names
cmNci60Names <- colnames(exprs(expData))
stopifnot(identical(cmNci60Names, stringr::str_trim(cmNci60Names)))

expectedCellLineNames <- readLines("inst/extdata/cellminer_2_0/CellMinerNci60Line_2_0_Names.txt")
stopifnot(identical(cmNci60Names, expectedCellLineNames))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: GENE COPY.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [DNA: Combined aCGH, select: Gene summary]
# Processing of data file (DNA__Combined_aCGH_Gene_summary.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

filePath <- "inst/extdata/cellminer_2_0/DNA__Combined_aCGH_Gene_summary.txt"
copTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
copTabOrig <- validateDataTab(copTabOrig, keyCol = "Probe id",
                              featureDataColNums = featureDataCols)

copData <- ExpressionSet(as.matrix(copTabOrig[, -featureDataCols]))
featureData(copData) <- new("AnnotatedDataFrame", data=copTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(copData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: EXOME/MUTATION.
#--------------------------------------------------------------------------------------------------

#----[gene level function altering mutations]-------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [DNA: Exome Seq, select: protein function affecting]
# Processing of data file (DNA__Exome_Seq_Protein_function_affecting.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

filePath <- "inst/extdata/cellminer_2_0/DNA__Exome_Seq_Protein_function_affecting.txt"
mutTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
mutTabOrig <- validateDataTab(mutTabOrig, keyCol = "Gene name",
                              featureDataColNums = featureDataCols)

# NAs indicate that not enough reads were available to determine variant allele percent
# conversion (see CellMiner spreadsheet footnotes); treated as zeros for analyses.
for(cLine in colnames(mutTabOrig[, -featureDataCols])){
  naIndexSet <- which(is.na(mutTabOrig[, cLine]))
  mutTabOrig[naIndexSet, cLine] <- 0
}

mutData <- ExpressionSet(as.matrix(mutTabOrig[, -featureDataCols]))
featureData(mutData) <- new("AnnotatedDataFrame", data=mutTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(mutData)), cmNci60Names))

#----[variant level exome sequencing data]--------------------------------------------------
filePath <- "inst/extdata/cellminer_2_0/DNA__Exome_Seq_none.txt"
exoTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:18
exoTabOrig <- validateDataTab(exoTabOrig, keyCol = "Probe id",
                              featureDataColNums = featureDataCols)

exoData <- ExpressionSet(as.matrix(exoTabOrig[, -featureDataCols]))
featureData(exoData) <- new("AnnotatedDataFrame", data=exoTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(exoData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: CELL LINE METADATA.
#--------------------------------------------------------------------------------------------------

filePath <- "inst/extdata/cellminer_2_0/CELLMINER_CELL_LINE_METADATA.txt"
mdaTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings=c("", "NA", "?"))

# Fix cell line names
mdaTabOrig$`Cell Line Name` <- expectedCellLineNames
rownames(mdaTabOrig) <- mdaTabOrig$`Cell Line Name`

# Keep a limited features, for the example
quantFeatures <- c("p53", "doubling time")
mdaQuantTab <- mdaTabOrig[, quantFeatures]
colnames(mdaQuantTab) <- c("IS_P53_MUT", "DOUBLING_TIME")

# Numeric features tend to be preferred, but may not be necessary dependent on the task
mdaQuantTab$IS_P53_MUT[stringr::str_trim(mdaQuantTab$IS_P53_MUT) == "MT"] <- 1
mdaQuantTab$IS_P53_MUT[stringr::str_trim(mdaQuantTab$IS_P53_MUT) == "WT"] <- 0
mdaQuantTab$IS_P53_MUT <- as.integer(mdaQuantTab$IS_P53_MUT)

mdaQuantTab$DOUBLING_TIME <- as.numeric(mdaQuantTab$DOUBLING_TIME)

mdaTabSampleInfo <- mdaTabOrig[, setdiff(colnames(mdaTabOrig), quantFeatures)]

stopifnot(all(c(lapply(mdaQuantTab, is.numeric), recursive = TRUE)))
mdaData <- ExpressionSet(t(mdaQuantTab))
stopifnot(is.numeric(exprs(mdaData)))

mdaAnnot <- data.frame(Name = rownames(exprs(mdaData)), Footnote = NA, stringsAsFactors = FALSE)
rownames(mdaAnnot) <- mdaAnnot$Name
mdaAnnot["IS_P53_MUT", "Footnote"] <- "p53 status as determined by yeast growth functional assay: PM O'Conner, et al. (Cancer Res. 1997 Oct 1;57(19):4285-300)."
mdaAnnot["DOUBLING_TIME", "Footnote"] <- "Doubling times described at NCI/DTP site."

featureData(mdaData) <- new("AnnotatedDataFrame", data=mdaAnnot)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(mdaData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: DRUG ACTIVITY.
#--------------------------------------------------------------------------------------------------

# activity data -------------------------------------------------------------------------
filePath <- "inst/extdata/cellminer_2_0/DTP_NCI60_ZSCORE.txt"
actTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="",
                         na.strings=c("", "na", "-"))

# This re-orders the columns to put the drug information together
actTabOrig <- actTabOrig[, c(1:6, 67, 68, 7:66)]

# Keep just the drug activity data for actData
featureDataCols <- 1:8
actTabOrig <- validateDataTab(actTabOrig, keyCol = "NSC #",
                              featureDataColNums = featureDataCols)
actData <- ExpressionSet(as.matrix(actTabOrig[, -featureDataCols]))

# Keep a limited features, for the example
featureDataCols <- 1:3
drugInfoTab <- actTabOrig[, featureDataCols]
colnames(drugInfoTab) <- c("NSC", "NAME", "FDA_STATUS")
drugInfoTab$NSC <- as.character(drugInfoTab$NSC)

stopifnot(identical(rownames(exprs(actData)), rownames(drugInfoTab)))
featureData(actData) <- new("AnnotatedDataFrame", data=drugInfoTab)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(actData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# Make NCI-60 sample info (shared by molData and drugData objects to be constructed).
#--------------------------------------------------------------------------------------------------
stopifnot(identical(mdaTabSampleInfo$`Cell Line Name`, cmNci60Names))

nci60Miame <- new("MIAME", name="CellMiner 2.0", lab="NCI/DTB",
                  samples=list(Name = cmNci60Names,
                               Gender = mdaTabSampleInfo[, "sex"],
                               PriorTreatment = mdaTabSampleInfo[, "prior treatment"],
                               Histology = mdaTabSampleInfo[, "histology"],
                               Source = mdaTabSampleInfo[, "source"],
                               Ploidy = mdaTabSampleInfo[, "ploidy"],
                               Institute = mdaTabSampleInfo[, "Institute"],
                               Contributor = mdaTabSampleInfo[, "Contributor"],
                               Reference = mdaTabSampleInfo[, "Reference"]))

#--------------------------------------------------------------------------------------------------
# Make NCI-60 MolData object.
#--------------------------------------------------------------------------------------------------

nci60ESetList <- list()

nci60ESetList[["exp"]] <- expData
nci60ESetList[["cop"]] <- copData
nci60ESetList[["mut"]] <- mutData

nci60ESetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")

#--------------------------------------------------------------------------------------------------
# Make NCI-60 DrugData object.
#--------------------------------------------------------------------------------------------------

# NOTE: For simplicity of the example, the actData variable is used again for the repeat data
drugData <- new("DrugData", act = actData, repeatAct = actData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")

#--------------------------------------------------------------------------------------------------
# UPDATE SAMPLE DATA FOR DRUG DATA OBJECT
#--------------------------------------------------------------------------------------------------
library(rcellminer)
rcmMolData  <- rcellminerData::molData  # rcellminer MolData object
rcmDrugData <- rcellminerData::drugData # rcellminer DrugData object

stopifnot(identical(rcellminer::getSampleData(rcmMolData)[["Name"]],
                    rcellminer::getSampleData(rcmDrugData)[["Name"]]))

drugData <- new("DrugData", 
                act = rcmDrugData@act,     # keep existing object
                repeatAct = rcmDrugData@repeatAct, # keep existing object
                sampleData = rcmMolData@sampleData)    # take correct molData object version

save(drugData, file = "data/drugData.RData")
#--------------------------------------------------------------------------------------------------
