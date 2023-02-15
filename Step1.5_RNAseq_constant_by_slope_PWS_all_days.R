

Data.cpm <- cpm(Data, log=FALSE)
Data.cpmC <- Data.cpm + 2.5
Data.lcpm <- log(Data.cpm, 2)
Data.lcpmC <- log(Data.cpmC, 2)


#################################
### Test which constant produced a slope approachig zero.

js <- c(0, 18.6220 ) #pws timeseries w/o D2
js <- c(0, 17.80592) #PWS timeseries all days

js <- c(0.001,0.01,0.1,0.5)

js <- c(2.1,2.5,3,3.5)

js <- c(2.8,2.9,3.1,3.2)

js <- seq(from = 0.399921, to = 0.4001, by = 0.00001)
js <- seq(from = 2.899921, to = 2.9001, by = 0.00001)
js <- seq(from = 2.859, to = 2.88, by = 0.001)
js <- seq(from = 0, to = 1, by = 0.1)

x <- 17.8059
y <- 17.8060
z <- .00001

js <- seq(from = x, to = y, by = z)

js

slopes <- as.data.frame (c())

for (j in 1:length(js))  {
  #Transform count + constant:
  print (js[j])
  print(j)
  print(js)

  constant <- js[j]
  Data.cpmC <- Data$counts + constant
  trans <- cpm(Data.cpmC, log=TRUE)

  mean_counts <- colMeans(trans)
  norm_counts <- sweep(trans, 2, colMeans(trans), "-")
  grandmean <- mean(mean_counts)
  norm_counts <- norm_counts + grandmean
  #head(norm_counts)
  range(norm_counts) #[1] -0.0311301 19.8623468

  ttrans <- t(head(norm_counts))
  ttrans <- t(norm_counts)
  ttrans <- as.data.frame(ttrans)

  genelist <- colnames(ttrans)
  ttrans$trt <- Data$samples$dpf
  treat <- unique(ttrans$trt)

  meandat <- matrix(nrow=length(treat), ncol=length(genelist))
  colnames(meandat) <- genelist
  rownames(meandat) <- treat
  for (i in 1:length(treat)) {
    a <- which(ttrans$trt == treat[i])
    n <- apply(ttrans[a,1:(dim(ttrans)[2]-1)], 2, mean)
    b <- which(rownames(meandat)==treat[i])
    meandat[b,] <- n
  }

  #swap this for sd instead of variance

  vardat <- matrix(nrow=length(treat), ncol=length(genelist))
  colnames(vardat) <- genelist
  rownames(vardat) <- treat
  for (i in 1:length(treat)) {
    if (length(which(ttrans$trt==treat[i])) > 1) {
      a <- which(ttrans$trt==treat[i])
      #print(a)
      n <- apply(ttrans[a,1:(dim(ttrans)[2]-1)], 2, var, na.rm=TRUE)

      test <- ttrans[a,1:(dim(ttrans)[2]-1)]

      #print(n)
      b <- which(rownames(vardat)==treat[i])
      #print(b)
      vardat[b,] <- n
    }}

  stdev <- sqrt(vardat)
  mv <- list(mean=meandat, stdev=stdev)

  #mv <- MeanVar(ttrans)
  mea <- as.vector(mv$mean)
  sd <- as.vector(mv$stdev)

  slope <- (coef(lm(sd~mea))[2])
  print(slope)
  slopes <- rbind (slopes, c (js[j], slope[[1]]))

  plot(mea,sd, pch=16, cex=0.5, main= paste0("Mean-SD relationship of Normalized Log(counts + ", js[j], " slope = ", slope))
  abline(coef(lm(sd~mea)),col="red",lwd=3)
  lines(lowess(mea,sd,f=0.2),col="blue",lwd=3)

  }

colnames(slopes) <- c("constant", "slope")
slopes
