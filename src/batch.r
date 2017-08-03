
#require(devtools)
#install_github("theislab/kBET")
require(kBET)
suppressPackageStartupMessages(library(Rtsne))
require(ggplot2)
require(gridExtra)

data <- read.table(file = '/home/jovyan/jamboree/kolodziejczyk_2015/counts//kolodziejczk_counttable_es.tsv', sep='\t')
annotation <- read.table(file='/home/jovyan/jamboree/kolodziejczyk_2015//counts//kolodziejczk_annotations.tsv', 
                         sep='\t', header=TRUE)
mat <- as.matrix(data)

fidx <- apply(mat, 1, function(x){ any(x > 50) })
tmp <- mat[fidx, ]
tmp <- log10(tmp + 1)
control <- kBET(tmp[,annotation$medium=='2i'], batch = annotation$batch[annotation$medium=='2i'])$average.pval
message(control)


control <- kBET(pcRegressed[,annotation$medium=='2i'], batch = annotation$batch[annotation$medium=='2i'])$average.pval
message(control)
# transformation

rmPCbet <- function(tmp, keep = 0, anno){
  # return kBET results for the first transformed PCs  
  pcCor <- rep(NA, times = 5) 
    for(i in 1:5){
      keep <- i
     pca_tmp <- prcomp(tmp)
     pcrm <- pca_tmp$x[,keep, drop = FALSE] %*% t( pca_tmp$rotation[,keep, drop = FALSE]) 
     pca <- prcomp(pcrm)
    
     pcCor[i] <- kBET(pcrm, batch = anno$batch)$average.pval
     
    }
    return(pcCor)
}

# run the regression for each biological cluster
bioParts <- list()
for(i in levels(annotation$medium)){
  level <- i
  res <- rmPCbet(tmp[, annotation$medium==i], anno = annotation[annotation$medium==i, ])
  keep <- which.max(res)
  pca_tmp <- prcomp(tmp[, annotation$medium==i])
  pcrm <- pca_tmp$x[,keep, drop = FALSE] %*% t( pca_tmp$rotation[,keep, drop = FALSE]) 
  bioParts[[length(bioParts) + 1]] <- pcrm
}

# combine the regressed parts
pcRegressed <- do.call(cbind, bioParts )


pca <- prcomp(pcRegressed)

da <- data.frame(dim1 = pca$rotation[,1 ], dim2 = pca$rotation[,2 ],
                 Batch = as.factor(annotation$batch),
                 Real = annotation$medium)
ggb <- ggplot(da, aes(x = dim1, y = dim2)) + geom_point(alpha = 0.75,
                 aes(colour = Batch)) + ggtitle('Mouse 3 Media (Colour = Batch)')
ggr <- ggplot(da, aes(x = dim1, y = dim2)) + geom_point(alpha = 0.75,
                                                       aes(colour = Real)) +
                 ggtitle('Mouse 3 Media (Colour = Media)')

grid.arrange(ggb, ggr, ncol = 2, top = 'Removal of Batch Effect Correlating PCs')

#pseu <- 1 + abs(min(pcrm))
#tSNE.2i <- Rtsne(t(log10(pseu+pcrm)))
#tSNE.df <- data.frame(tSNE1 = tSNE.2i$Y[,1], 
#                      tSNE2=tSNE.2i$Y[,2], 
#                      medium = annotation$medium,
#                      batch= as.factor(annotation$batch))
#gg <- ggplot(tSNE.df, aes(tSNE1, tSNE2, color=batch)) + geom_point(alpha=0.5) + theme_bw()
#gg <- ggplot(tSNE.df, aes(tSNE1, tSNE2, color=medium)) + geom_point(alpha=0.5) + theme_bw()

#keep <- c(which(res > 0.8))
#keep <- c(which(res > 0))
#pca_tmp <- prcomp(tmp)
#pcrm <- pca_tmp$x[,keep, drop = FALSE] %*% t( pca_tmp$rotation[,keep, drop = FALSE]) 

#plot(gg)
#res <- kBET(pcrm, batch = annotation$batch)
#return(res$average.pval)
#}

#res <- rmPCbet(tmp, annotation = annotation)
#mean(entPro)

#p <- kBET(pca$rotation[, 1:2], batch = annotation$batch)$average.pval

#da <- data.frame(dim1 = pca$rotation[,1 ], dim2 = pca$rotation[,2 ],
#                 Batch = as.factor(annotation$batch),
#                 Real = annotation$medium)
#gg <- ggplot(da, aes(x = dim1, y = dim2)) + geom_point(alpha = 0.75,
#                 aes(colour = Batch, shape = Real)) + ggtitle(keep)
#  ggtitle(paste('rm::', n  ))

#pseu <- 1 + abs(min(pcrm))
#tSNE.2i <- Rtsne(t(log10(pseu+pcrm[, annotation$medium=='2i'])))
#tSNE.df <- data.frame(tSNE1 = tSNE.2i$Y[,1], 
#                      tSNE2=tSNE.2i$Y[,2], 
#                      batch= as.factor(annotation$batch[annotation$medium=='2i']))
#gg <- ggplot(tSNE.df, aes(tSNE1, tSNE2, color=batch)) + geom_point(alpha=0.5) + theme_bw()
#gg <- gg + ggtitle(keep)
#plot(gg)
#res <- kBET(pcrm, batch = annotation$batch)
#return(res$average.pval)



