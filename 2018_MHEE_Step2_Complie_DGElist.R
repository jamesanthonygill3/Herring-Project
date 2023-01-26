#=================== Upload DGElist and remove bad data and low count genes =====================
#see previous versions for details on how DGElist was created with Sample.key.
library(edgeR)

load("/Box/ph18_V2_TC/results/tx2geneV3/MHEE_DGElist3a.RDS" )

Days <- c("D4", "D6", "D10", "All")
################ work on annotatin this better##########################

if (!dir.exists( "/Box/ph18_V2_TC/results/tx2geneV3" )) {
  dir.create ("/Box/ph18_V2_TC/results/tx2geneV3")
  print ("Made Folder 1")
}

for (i in 1:4) {

  Data <- MHEE.DGEList
  Day <- Days[i]
  print(Day)

  for (j in 1:3) {
    if (j == 1) {
      if (!Day == "All") {
        population <- c("PW", "S", "CP")
        folder <- paste0( "/Box/ph18_V2_TC/results/tx2geneV3/2018_MHEE_", Day, "/")
        Filename <- paste0("2018_MHEE_", Day, "_lcpm_V2.0")
      }}
    if (j == 2) {
      if (!Day == "All") {
        population <- c("PW", "S")
        folder <- paste0( "/Box/ph18_V2_TC/results/tx2geneV3/2018_MHEE_", Day, "_AK_pops/")
        Filename <- paste0("2018_MHEE_", Day, "_lcpm_V2.0_AK_only" )
      }}
    if (j == 3) {
      if (Day == "All") {
        Day         <- levels(Data$samples$dpf)
        population  <- c("PW")
        folder      <- paste0( "/Box/ph18_V2_TC/results/tx2geneV3/2018_MHEE_PW_timeseries/")
        Filename    <- paste0("2018_MHEE_PW_timeseries_lcpm_V2.0")
        Data <- Data [ ,!grepl( pattern = "-D2-", rownames(Data$samples))]
      } else { #only use this code when Day == "All"
        break ()
      }}

    print(population)
    if (!dir.exists(folder)){
      dir.create(folder)
      print ("Made Folder 2")}

    #all treatments
    treatmnets <- paste0 (paste0 ("-", levels (Data$samples$trt)))
    treatmnets

    Data <- Data [,grepl( pattern = paste0(Day, collapse = "|"), Data$samples$Sequence.name)]
    Data <- Data [,grepl( pattern = paste0(population, collapse = "|"), Data$samples$Sequence.name)]
    Data <- Data [,grepl( pattern = paste0(treatmnets, collapse = "|"), Data$samples$Sequence.name)]

    dim(Data)

    #remove samples with low lib size.
    Data <- Data[,Data$samples$lib.size > 7.5 * 10^5]
    dim(Data)

    #----- remove outliers ----------------------
    Data <- Data [,!grepl( pattern = "-M3", rownames(Data$samples))]

    #Day 10
    Data <- Data [,!grepl( pattern = "S-D10-ML4", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "S-D10-M2A", rownames(Data$samples))]

    Data <- Data [,!grepl( pattern = "CP-D10-L1A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "CP-D10-MH3A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "CP-D10-ML4A", rownames(Data$samples))]  #outlier on MDS Dim 2
    Data <- Data [,!grepl( pattern = "CP-D10-C2B", rownames(Data$samples))]  #outlier on MDS Dim 1

    Data <- Data [,!grepl( pattern = "CP-D10-L3A", rownames(Data$samples))]  #outlier on MDS Dim 1

    Data <- Data [,!grepl( pattern = "PW-D10-H2B", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D10-ML2B", rownames(Data$samples))]

    #Day 6
    Data <- Data [,!grepl( pattern = "PW-D6-H4A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D6-M1B", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D6-M1B", rownames(Data$samples))]

    Data <- Data [,!grepl( pattern = "CP-D6-C2A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "CP-D6-H4A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "CP-D6-MH4A", rownames(Data$samples))]

    Data <- Data [,!grepl( pattern = "S-D6-H3A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "S-D6-ML4A", rownames(Data$samples))]

    #Day 8
    Data <- Data [,!grepl( pattern = "PW-D8-H4B", rownames(Data$samples))]

    #Day 4
    Data <- Data [,!grepl( pattern = "PW-D4-L3A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D4-ML4A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D4-H3A", rownames(Data$samples))]

    Data <- Data [,!grepl( pattern = "S-D4-C3A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "S-D4-MH1A", rownames(Data$samples))]

    #PWS other times
    Data <- Data [,!grepl( pattern = "PW-D2-ML4A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D2-C2A", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D2-C4B", rownames(Data$samples))]
    Data <- Data [,!grepl( pattern = "PW-D2-ML2A", rownames(Data$samples))]

    #------------------ filter & normalize -------------------

    par(mfrow=c(1,2))
    hist(aveLogCPM(Data))

    keep <- filterByExpr (y = Data, group = Data$samples$group, min.count = 35, min.total.count = 500, large.n = 6, min.prop = 0.7)
    table(keep)
    Data <- Data [keep,,keep.lib.sizes=FALSE]
    rm (keep)
    dim(Data)
    hist(aveLogCPM(Data))

    output <- paste0(folder, "/", Filename, "_hist.png")
    dev.print(png, filename = output,  width = 2500, height = 1250, res = 150)
    dev.off()


    library(RColorBrewer)
    nsamples <- ncol(Data)
    col <- brewer.pal(nsamples, "Paired")

    par(mfrow=c(1,1))
    lcpm <- cpm(Data, log=TRUE)
    boxplot(lcpm, las=2, col=col, main="")
    title(main="A. Unnormalised data", ylab="Log-cpm")

    output <- paste0(folder, "/", Filename, "_Unnormalized.png")
    dev.print(png, filename = output,  width = 2500, height = 1250)
    dev.off()


    Data <- calcNormFactors(Data, method = "TMM") # "TMM", "RLE", or "upperquartile" or "none".

    par(mfrow=c(1,1))
    lcpm <- cpm(Data, log=TRUE)
    boxplot(lcpm, las=2, col=col, main="")
    title(main="B. Normalised data", ylab="Log-cpm")

    output <- paste0(folder, "/", Filename, "_Normalized.png")
    dev.print(png, filename = output,  width = 2500, height = 1250)
    dev.off()


    Data.cpm <- cpm(Data, log=FALSE)
    Data.cpmC <- Data.cpm + 2.5
    Data.lcpmC <- log(Data.cpmC, 2)

    rm(nsamples)
    rm(col)

    #============ MDS plots =================================================================================================================

    if (j != 3) {
      Header <- paste0("2018 Marrowstone ", Day, ", pops = ", paste(population, collapse = ", "), ", by plate ")
      par(mfrow=c(2,3))
      for (h in 1:3) {
        for (k in 1:3) {
          if (k != h) {
            plotMDS (as.data.frame(Data.cpm), col = Data$samples$Color.plate, labels = Data$samples$pop.trt.tank.rep, dim.plot = c(h,k), gene.selection = "pairwise")
            title(main = Header)
          }
        }}

      output <- paste0(folder, "/", Filename, "_MDSbyPlate.png")
      dev.print(png, filename = output, width = 2500, height = 1250, res = 150)
      dev.off()

      Header <- paste0("2018 Marrowstone ", Day, ", pops = ", paste(population, collapse = ", "), ", by trt ")
      par(mfrow=c(2,3))
      for (h in 1:3) {
        for (k in 1:3) {
          if (k != h) {
            plotMDS (as.data.frame(Data.cpm), col = Data$samples$Color.trt, labels = Data$samples$pop.trt.tank.rep, dim.plot = c(h,k), gene.selection = "pairwise")
            title(main = Header)
          }
        }}
      output <- paste0(folder, "/", Filename, "_MDSbyTrt.png")
      dev.print(png, filename = output, width = 2500, height = 1250, res = 150)
      dev.off()


    } else { #(j = 3)
      Header <- paste0("2018 Marrowstone timecourse: pops = ", paste(population, collapse = ", "), ", by plate ")
      par(mfrow=c(2,3))
      for (h in 1:3) {
        for (k in 1:3) {
          if (k != h) {
            if ( j == 3 ) {
              plotMDS (as.data.frame(Data.cpm), col = Data$samples$Color.plate, labels = Data$samples$dpf.pop.trt.tank, dim.plot = c(h,k), gene.selection = "pairwise")
              title(main = Header)
            }}}}
      output <- paste0(folder, "/", Filename, "_MDSbyPlate.png")
      dev.print(png, filename = output, width = 2500, height = 1250, res = 150)
      dev.off()


      Header <- paste0("2018 Marrowstone timecourse: pops = ", paste(population, collapse = ", "), ", by trt ")
      par(mfrow=c(2,3))
      for (h in 1:3) {
        for (k in 1:3) {
          if (k != h) {
            if ( j == 3 ) {
              plotMDS (as.data.frame(Data.cpm), col = Data$samples$Color.trt, labels = Data$samples$dpf.pop.trt.tank, dim.plot = c(h,k), gene.selection = "pairwise")
              title(main = Header)
            }}}}
      output <- paste0(folder, "/", Filename, "_MDSbyPlate.png")
      dev.print(png, filename = output, width = 2500, height = 1250, res = 150)
      dev.off()


      Header <- paste0("2018 Marrowstone timecourse: pops = ", paste(population, collapse = ", "), ", by dpf ")
      par(mfrow=c(2,3))
      for (h in 1:3) {
        for (k in 1:3) {
          if (k != h) {
            if ( j == 3 ) {
              plotMDS (as.data.frame(Data.cpm), col = Data$samples$Color.dpf, labels = Data$samples$dpf.pop.trt.tank, dim.plot = c(h,k), gene.selection = "pairwise")
              title(main = Header)
            }}}}
      output <- paste0(folder, "/", Filename, "_MDSbyDPF.png")
      dev.print(png, filename = output, width = 2500, height = 1250, res = 150)
      dev.off()

    }

    #=====  prepare data and run linear model ======================================================================================

    Data$samples$trt <- factor(Data$samples$trt)
    Data$samples$pop <- factor(Data$samples$pop)
    Data$samples$dpf <- factor(Data$samples$dpf)
    Data$samples$plate <- factor (Data$samples$plate)

    levels(Data$samples$trt)
    levels(Data$samples$pop)
    levels(Data$samples$dpf)
    levels(Data$samples$plate)

    #Set up data for lm():
    #ttrans <- t(Data.cpm)  #Normalized counts
    ttrans <- t(Data.lcpmC)  #Normalized count with the 2.5 constant to each count
    ttrans <- as.data.frame(ttrans)

    ttrans$trt <- Data$samples$trt
    ttrans$pop <- Data$samples$pop
    ttrans$dpf <- Data$samples$dpf
    ttrans$plate <- as.numeric (Data$samples$plate)

    #ttrans$bb <- Data$samples$body.burden
    ttrans$bb <- log (Data$samples$body.burden, 2)

    dim(ttrans)

    #Set up data frame to store p-values:
    d <- dim(ttrans)[2]- 5 # remove the number of factors from d
    glist.y <- colnames(ttrans[,1:d])

    if (j == 1) {
      cofs <- c("Model.EFFECT", "Main.bb.aov", "Main.pop.aov", "Main.popXbb.aov", "r2", "Main.(Intercept)", "Main.bb.lm", "Main.S", "Main.PW", "Main.bb:S", "Main.bb:PW", "Estimate.Intercept", "Estimate.bb", "Estimate.S", "Estimate.PW", "Estimate.bb:S", "Estimate.bb:PW", "Tukey.CP", "Tukey.PW", "Tukey.S")
    } else {
      if (j == 2) {
        cofs <- c("Model.EFFECT", "Main.bb.aov", "Main.pop.aov", "Main.popXbb.aov", "r2", "Main.(Intercept)", "Main.bb.lm", "Main.PW",  "Main.bb:PW", "Estimate.Intercept", "Estimate.bb", "Estimate.PW",  "Estimate.bb:PW", "Tukey.PW", "Tukey.S")
      } else {
        if (j == 3) {
          Mains     <- c( "Main.(D2.Intercept)", "Main.bb.lm", "Main.D4",  "Main.D6", "Main.D8", "Main.D10",
                          "Main.bbXD2", "Main.bbXD4, Main.bbXD6", "Main.bbXD8", "MainbbXD10" )
          Estimates <- c( "Estimate.(D2.Intercept)", "Estimate.bb.lm", "Estimate.D4",  "Estimate.D6", "Estimate.D8", "Estimate.D10",
                          "Estimate.bbXD4", "Estimate.bbXD6", "Estimate.bbXD8", "Estimate.bbXD10" )
          Tukeys    <- c( "Tukey.D2", "Tukey.D4", "Tukey.D6", "Tukey.D8", "Tukey.D10")
          cofs <- c("Model.EFFECT", "Main.bb.aov", "Main.dpf.aov", "Main.dpfXbb.aov", "r2", Mains, Estimates, Tukeys)
          rm(Mains)
          rm(Estimates)
          rm(Tukeys)
        }}}

    pdat.y <- data.frame (matrix(nrow=length(glist.y), ncol=length(cofs)))
    colnames(pdat.y) <- cofs
    rownames(pdat.y) <- glist.y

    rm(cofs)
    rm(glist.y)

    library(agricolae)

    #Run lm() on each gene:
    if (j == 3) {
      print(paste0(Day, "_", population))
      for (h in 1:ncol(ttrans[,1:d])) {
        gene <- which(rownames(pdat.y) == colnames(ttrans)[h])
        #type III sum of squares for the model effect
        Model <- lm(formula=ttrans[,h] ~ bb * dpf  + (1 | plate), data=ttrans)

        #ANOVA for model effects
        f = summary(Model)$fstatistic
        Model.p <- pf(f[1],f[2],f[3],lower.tail=F)
        Model.p

        a <- aov(Model)
        a.sum <- summary (a)[[1]]
        Main.p <- a.sum$`Pr(>F)`[1:3]
        r2 <- summary(Model)$adj.r.squared
        a.sum

        # regression analysis
        M.sum <- summary(Model)
        M.sum

        Estimates <- M.sum$coefficients[,1]
        Main.M <- M.sum$coefficients[,4]

        Tukey.dpf <- HSD.test (a, "dpf")
        Tukey.dpf <- Tukey.dpf$groups
        colnames(pdat.y)
        pdat.y[gene,] <- c (Model.p, Main.p, r2, Main.M, Estimates, Tukey.dpf$groups)
      }
    } else {
      print(paste0(Day, "_", population))
      for (h in 1:ncol(ttrans[,1:d])) {
        gene <- which(rownames(pdat.y) == colnames(ttrans)[h])
        #type III sum of squares for the model effect
        Model <- lm(formula=ttrans[,h] ~ bb * pop + (1 | plate), data=ttrans)
        M.sum <- summary(Model)
        M.sum
        #ANOVA for model effects
        f = summary(Model)$fstatistic
        Model.p <- pf(f[1],f[2],f[3],lower.tail=F)

        # regression analysis
        Estimates <- M.sum$coefficients[,1]
        Main.M <- M.sum$coefficients[,4]
        r2 <- summary(Model)$adj.r.squared

        a <- aov(Model)
        a.sum <- summary (a)[[1]]
        a.sum

        Main.p <- a.sum$`Pr(>F)`[1:3]

        Tukey.pop <- HSD.test (a, "pop")
        Tukey.pop <- Tukey.pop$groups
        Tukey.pop <- Tukey.pop [ order (row.names(Tukey.pop)), ]
        colnames(pdat.y)
        pdat.y[gene,] <- c (Model.p, Main.p, r2, Main.M, Estimates, Tukey.pop$groups)
      }}

    rm(d)
    rm(gene)
    rm(Main.p)
    rm(Model.p)
    rm(a)
    rm(Model)
    rm(ttrans)
    rm(f)
    rm(r2)
    rm(Tukey.pop)
    rm(Tukey.dpf)
    rm(a.sum)
    rm(M.sum)
    rm(Main.M)

    output <- paste0(folder, "/", Filename, "_pdat.csv")

    pdat.y <- as.data.frame(pdat.y)
    write.csv(pdat.y, output)

    #===========================================================================================
    #Start new run with upload

    pdat.y <- read.csv(output)
    colnames(pdat.y)[1] <- "Gene"

    #### Adjust p-values for FDR
    #Vectorize and apply p.adjust to pdat.y dataframe

    adj_pval_all <- pdat.y

    if (j == 3) {

      p.adj.Model <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Model.EFFECT, meth))
      p.adj.bb  <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Main.bb.aov, meth))
      p.adj.dpf  <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Main.dpf.aov, meth))
      p.adj.bbXdpf  <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Main.dpfXbb.aov, meth))

      p.adj.Model <- as.data.frame(p.adj.Model)
      p.adj.bb  <- as.data.frame(p.adj.bb)
      p.adj.dpf  <- as.data.frame(p.adj.dpf)
      p.adj.bbXdpf  <- as.data.frame(p.adj.bbXdpf)

      pResults <- cbind( p.adj.Model$fdr, p.adj.dpf$fdr, p.adj.bb$fdr, p.adj.bbXdpf$fdr, pdat.y )

      rm(p.adj.Model)
      rm(p.adj.bb)
      rm(p.adj.dpf)
      rm(p.adj.bbXdpf)

    } else {

      p.adj.Model <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Model.EFFECT, meth))
      p.adj.bb  <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Main.bb.aov, meth))
      p.adj.pop  <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Main.pop.aov, meth))
      p.adj.bbXpop  <- sapply(p.adjust.methods, function(meth) p.adjust(adj_pval_all$Main.popXbb.aov, meth))

      p.adj.Model <- as.data.frame(p.adj.Model)
      p.adj.bb  <- as.data.frame(p.adj.bb)
      p.adj.pop  <- as.data.frame(p.adj.pop)
      p.adj.bbXpop  <- as.data.frame(p.adj.bbXpop)

      pResults <- cbind( p.adj.Model$fdr, p.adj.pop$fdr, p.adj.bb$fdr, p.adj.bbXpop$fdr, pdat.y )
      Results <- pResults [ order(pResults$r2), ]
      dim (Results)

      rm(p.adj.Model)
      rm(p.adj.bb)
      rm(p.adj.pop)
      rm(p.adj.bbXpop)
    }


    Results.sort <- pResults [ order(pResults$r2), ]
    dim (Results.sort)
    rm(pResults)
    rm(adj_pval_all)

    #================ Generate data for graphing with removeBatchEffect ========================

    if (j == 3) {
      MHEE.design <- model.matrix(~ 0 + dpf + dpf:log(body.burden, 2), data = Data$samples)
    } else {
      MHEE.design <- model.matrix(~ 0 + pop + pop:log(body.burden, 2), data = Data$samples)
    }

    colnames (MHEE.design)
    v <- voom(Data, MHEE.design, plot=FALSE)
    v$E <- removeBatchEffect(v$E, batch = Data$samples$plate)
    rm(MHEE.design)

    #==========================================================================================================================
    library(pheatmap)

    Group <- paste (Data$samples$dpf, Data$samples$pop, Data$samples$trt, sep = "_")

    #make dataframe with mean data by factor (each trt in each pop)
    Data.mean.by.factor <- t(v$E)
    Data.mean.by.factor <- aggregate(x = Data.mean.by.factor, by = list(Group), FUN = mean)
    row.names(Data.mean.by.factor) <- Data.mean.by.factor$Group.1
    Data.mean.by.factor <- Data.mean.by.factor[,-1]
    Data.mean.by.factor <- as.data.frame (t(Data.mean.by.factor))

    Data.mean.by.factor.DDct <- as.data.frame (rownames(Data$counts))
    row.names(Data.mean.by.factor.DDct) <- Data.mean.by.factor.DDct[,1]

    colnames(Data.mean.by.factor)

    Group <- as.data.frame (Group [ !duplicated (Group)])
    colnames (Group) <- "Group"
    Group$dpf <- gsub ( Group, pattern = "(\\S+)_(\\S+)_(\\S+)", replacement = "\\1", x = Group$Group )
    Group$pop <- gsub ( Group, pattern = "(\\S+)_(\\S+)_(\\S+)", replacement = "\\2", x = Group$Group )
    Group$trt <- gsub ( Group, pattern = "(\\S+)_(\\S+)_(\\S+)", replacement = "\\3", x = Group$Group )
    Group$trt <- factor (Group$trt, levels = c("C","L", "ML", "M", "MH", "H"))

    if (j == 3) {
      #all samples Normalized by control mean for each pop
      Group$dpf <- factor (Group$dpf, levels = levels(Data$samples$dpf))
      for (h in 1:length(levels(Group$dpf)))  {
        subset.df <- Data.mean.by.factor [ , Group$dpf == levels(Group$dpf)[h]]
        subset.trt <- Group [ Group$dpf == levels(Group$dpf)[h], ]
        subset.ddt <- as.data.frame (subset.df - subset.df[, subset.trt$trt == "C" ])
        subset.ddt <- subset.ddt [ , order (subset.trt$trt) ]
        Data.mean.by.factor.DDct <- cbind(Data.mean.by.factor.DDct, subset.ddt)
      }
      Data.mean.by.factor <- Data.mean.by.factor [ ,order (Group$dpf, Group$trt) ]

    } else {
      #all samples Normalized by control mean for each pop
      for (h in 1:length(unique(Group$pop)))  {
        subset.trt <- Group [ Group$pop == unique(Group$pop)[h], ]
        subset.df <- Data.mean.by.factor [ , Group$pop == unique(Group$pop)[h]]
        subset.ddt <- as.data.frame (subset.df - subset.df[, subset.trt$trt == "C" ])
        subset.ddt <- subset.ddt [ , order (subset.trt$trt) ]
        Data.mean.by.factor.DDct <- cbind(Data.mean.by.factor.DDct, subset.ddt)
      }}

    Data.mean.by.factor.DDct <- Data.mean.by.factor.DDct[, -1]

    rm(subset.trt)
    rm(subset.df)
    rm(subset.ddt)
    rm(v)

    #=============================================================================================================================

    #Make result file with lcpm_ddct data
    annotation <- read.csv("/Box/ph18_V2_TC/results/tx2geneV3/annotationV3.csv", na.strings = c("", "NA"))

    MHFT.results <- annotation[!duplicated(annotation$gene.name), ]
    dim(MHFT.results)
    dim(Results.sort)
    Results.sort$Gene

    MHFT.results <- merge(MHFT.results, Results.sort, by.x = "gene.name", by.y = "Gene")
    dim(MHFT.results)

    Data.mean.by.factor.DDct$Tx <- row.names(Data.mean.by.factor.DDct)
    Data.mean.by.factor$Tx <- row.names(Data.mean.by.factor)
    Data.means <- merge(Data.mean.by.factor, Data.mean.by.factor.DDct, by.x = "Tx", by.y = "Tx")
    dim (Data.means)

    MHFT.results <- merge(MHFT.results, Data.means, by.x = "gene.name", by.y = "Tx")
    dim(MHFT.results)

    MHFT.results <- MHFT.results [ order (MHFT.results$`p.adj.Model$fdr`),]
    MHFT.results <- MHFT.results[, -c(2,3)]
    dim(MHFT.results)

    output <- paste0(folder, Filename, "_Results.csv")

    write.csv(MHFT.results, output)

  }
}



