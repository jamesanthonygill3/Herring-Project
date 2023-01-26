
#this file creates a DGElist for the Marrowstone herring embryo exposure in 2018 for
#all 3 populaitons (Sitka, Prince William Sound, and Cherry Point) and 3 timepoints D4, D6, & D10.


#Step 1.  Iport the tximport data, txiv3a. 
#tx2geneV3a was used in the txiport and is based off of a list of transcripts and gene names

#upload DGElist
load("/Box/ph18_V2_TC/results/tx2geneV3/txiV3a.csv")


#-------------import Sample names -----------------------------------------------------------------------------------------------
txi.samples <- colnames(txiV3a$counts)
Sample.Key <- as.data.frame(colnames(txiV3a$counts))
colnames(Sample.Key) <- c("Sequence.name")
Sample.Key$order <- as.integer(rownames(Sample.Key))

#extract experimental conditions from names
Sample.Key$pop <- gsub(x = Sample.Key$Sequence.name, pattern = "(\\S+)-(\\S+)-(\\S+)", replacement = "\\1")
Sample.Key$dpf <- gsub(x = Sample.Key$Sequence.name, pattern = "(\\S+)-(\\S+)-(\\S+)", replacement = "\\2")

Sample.Key$trt.tank.rep <- gsub(x = Sample.Key$Sequence.name, pattern = "(\\S+)-(\\S+)-(\\S+)", replacement = "\\3")
Sample.Key$trt <- gsub(x = Sample.Key$trt.tank.rep, pattern = "(\\S+)(\\d+)(\\S+)", replacement = "\\1")
Sample.Key$tank <- gsub(x = Sample.Key$trt.tank.rep, pattern = "(\\S+)(\\d+)(\\S+)", replacement = "\\2")
Sample.Key$rep <- gsub(x = Sample.Key$trt.tank.rep, pattern = "(\\S+)(\\d+)(\\S+)", replacement = "\\3")

#set names for grouping later
Sample.Key$trt.tank          <- paste0(Sample.Key$trt, Sample.Key$tank)
Sample.Key$pop.trt.tank      <- paste0(Sample.Key$pop, Sample.Key$trt.tank)
Sample.Key$dpf.pop.trt       <- paste0(Sample.Key$dpf, Sample.Key$pop, Sample.Key$trt)
Sample.Key$dpf.pop.trt.tank  <- paste0(Sample.Key$dpf, Sample.Key$pop, Sample.Key$trt.tank)
Sample.Key$pop.trt.tank.rep  <- paste0(Sample.Key$pop, Sample.Key$trt.tank.rep)

#assign groups in key for DGElist
Sample.Key$Group <- factor(paste(Sample.Key$dpf, Sample.Key$pop, Sample.Key$trt, sep = ":" ))

#--------------------------------------------------------------------------------------------------------
#add body burden from key

bb.key <- read.csv("/Box/ph18_V2_TC/results/tx2geneV3/BodyBurden_Key.csv")
bb.key <- bb.key [ , -1]
colnames(bb.key) <- c("Group", "body.burden")
bb.key <- bb.key [!is.na(bb.key[,1]),]

Sample.Key$pop.trt <- paste0(Sample.Key$pop, "-", Sample.Key$trt)
Sample.Key <- merge(x = Sample.Key, y = bb.key, by.x = "pop.trt", by.y = "Group")

rm(bb.key)

Sample.Key <- Sample.Key[order(Sample.Key$order),]
Sample.Key <- Sample.Key[,-1]

#--------------------------------------------------------------------------------------------------------
#add plate from key

plate.rep <- read.csv("/Box/ph18_V2_TC/results/tx2geneV3/reps_per_plate.csv", header = T, na.strings = c("", NA))
plate.rep <- as.data.frame (plate.rep[!is.na(plate.rep[,1]), ])

rep.key <- data.frame (matrix(ncol = 2, nrow = 0))
colnames(rep.key) <- c("sample_id", "plate")

for (i in 1:5) {
  df <- as.data.frame (plate.rep [,i])
  df [,2] <- colnames (plate.rep)[i]
  colnames(df) <- c("sample_id", "plate")
  rep.key <- rbind (rep.key, df)
  rm(df)
}

rep.key <- rep.key[!is.na(rep.key[,1]), ]
rm(plate.rep)

Sample.Key <- merge(x = Sample.Key, y = rep.key, by.x = "Sequence.name", by.y = "sample_id")
Sample.Key <- Sample.Key[order(Sample.Key$order),]

rm(rep.key)

Sample.Key$plate <- gsub(pattern = "plate", replacement = "", Sample.Key$plate)

#------------------Graphing Color assignments --------------------------------------------------------------
#assign graphing colors by populations
Sample.Key$Color.pop <- Sample.Key$pop
Sample.Key$Color.pop <- gsub(pattern = "CP", replacement = "blue", Sample.Key$Color.pop)
Sample.Key$Color.pop <- gsub(pattern = "PW", replacement = "orange", Sample.Key$Color.pop)
Sample.Key$Color.pop <- gsub(pattern = "S", replacement = "red", Sample.Key$Color.pop)

#assign graphing colors by Treatment
Sample.Key$Color.trt <- Sample.Key$trt
Sample.Key$Color.trt <- gsub("C", "black", Sample.Key$Color.trt)
Sample.Key$Color.trt <- gsub("ML", "blue", Sample.Key$Color.trt)
Sample.Key$Color.trt <- gsub("MH", "red", Sample.Key$Color.trt)
Sample.Key$Color.trt <- gsub("M", "orange", Sample.Key$Color.trt)
Sample.Key$Color.trt <- gsub("L", "purple", Sample.Key$Color.trt)
Sample.Key$Color.trt <- gsub("H", "dark red", Sample.Key$Color.trt)

#assign graphing colors by dpf
Sample.Key$Color.dpf <- Sample.Key$dpf
Sample.Key$Color.dpf <- gsub("D2", "black", Sample.Key$Color.dpf)
Sample.Key$Color.dpf <- gsub("D4", "purple", Sample.Key$Color.dpf)
Sample.Key$Color.dpf <- gsub("D6", "blue", Sample.Key$Color.dpf)
Sample.Key$Color.dpf <- gsub("D8", "orange", Sample.Key$Color.dpf)
Sample.Key$Color.dpf <- gsub("D10", "red", Sample.Key$Color.dpf)

Sample.Key$Color.plate <- Sample.Key$plate
Sample.Key$Color.plate <- gsub("1", "blue", Sample.Key$Color.plate)
Sample.Key$Color.plate <- gsub("2", "red", Sample.Key$Color.plate)
Sample.Key$Color.plate <- gsub("3", "black", Sample.Key$Color.plate)
Sample.Key$Color.plate <- gsub("4", "purple", Sample.Key$Color.plate)
Sample.Key$Color.plate <- gsub("5", "orange", Sample.Key$Color.plate)

#=============================================================================================================================
#set Factors in Sample.Key
Sample.Key$trt <- factor(Sample.Key$trt, levels = c("C", "L", "ML", "M", "MH", "H"))
Sample.Key$dpf <- factor(Sample.Key$dpf, levels = c("D2", "D4", "D6", "D8", "D10"))
Sample.Key$pop <- factor(Sample.Key$pop, levels = c("CP", "S", "PW"))
Sample.Key$plate <- factor(Sample.Key$plate, levels = c("1", "2", "3", "4", "5"))

row.names(Sample.Key) <- Sample.Key$Sequence.name
Sample.Key$Group <- factor(paste(Sample.Key$dpf.pop.trt , sep = ":" ))

#=============================================================================================================================
#==== tximport from Mike Love
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR

#=============================================================================================================================
# counts were countsFromAbundance = "lengthScaledTPM" in tximport for txiV2.RData
library(edgeR)
MHEE.DGEList <- DGEList(counts = txiv3a$counts, group = Sample.Key$Group)
MHEE.DGEList$samples <- cbind(MHEE.DGEList$samples, Sample.Key)
dim(MHEE.DGEList)

save(MHEE.DGEList, file = "/Box/ph18_V2_TC/results/tx2geneV3/MHEE_DGElist3b.RDS")


