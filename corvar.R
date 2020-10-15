#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("4 argument must be supplied, \nUsage: Rscript --vanilla corvar.R <../yourpath/filename_of_genotype.gz> <../yourpath/filename_of_local_gebv.gz> <../yourpath/filename_of_segment_information.gz> <../yourpath_of_output>", call.=FALSE)
}

library(data.table)
library(WGCNA)

#---define variables
genofn <- args[1] #(e.g., test.geno.gz)
lgebvfn <- args[2] #(e.g., test.lgebv.gz)
segfn <- args[3] #(e.g., test.lgseg.gz)
outputpathfn <- args[4] #(e.g., test)

#---read in data---
cat(paste0('data reading started at ',Sys.time()),sep='\n')
geno <- fread(cmd=paste0('zcat ',genofn))
lgebv <- fread(cmd=paste0('zcat ',lgebvfn))
seg <- fread(cmd=paste0('zcat ',segfn))
#---define output log----
cat(paste0('Analysis started at ',Sys.time()),sep='\n')
write(paste('data reading finished at',Sys.time()),file=paste0(outputpathfn,'.log'),append=TRUE)
#--determine total number of traits for the analysis
NTrTot <- ncol(lgebv)-2
#---loop through segment
seglist <- list()
for (j in unique((unlist(seg[,which(colnames(seg)=='seg'),with=F])))){
lgebvdum <- lgebv[seg==j]
dumgeno <- geno[,c(1,unlist(seg[seg==j]$order+1)),with=F]
write(paste(paste0('seg',j,' started at'),Sys.time()),file=paste0(outputpathfn,'.log'),append=TRUE)
corlist <- list()
for (i in sprintf('%0.2d',1:NTrTot)){
lgebvdum1 <- lgebv[,c(1,2,grep(i,colnames(lgebv))),with=F]
lgebvdum2 <- merge(lgebvdum1,dumgeno,by='ID',sort=F)
corlist[[paste0('tr',i)]] <- var(lgebvdum2[[paste0('tr',i)]])*(WGCNA::cor(lgebvdum2[,-c(1:2)],use='pairwise.complete.obs')[1,]^2)
}
cordat <- as.data.frame(rowSums(do.call(cbind,corlist),na.rm=T))[-1,,drop=F]
colnames(cordat)[1] <- 'corvar'
rownames(cordat) <- as.numeric(gsub('V','',rownames(cordat)))-1
cordat1 <- merge(cordat,seg,by.x='row.names',by.y='order')
colnames(cordat1)[1] <- 'order'
write(paste(paste0('loop within seg',j,' finished at'),Sys.time()),file=paste0(outputpathfn,'.log'),append=TRUE)
seglist[[j]] <- cordat1
}
segdat <- do.call(rbind,seglist)
#---last update on log file
write(paste(paste0('all loops',' finished at'),Sys.time()),file=paste0(outputpathfn,'.log'),append=TRUE)
#---write out results
write.table(segdat,gzfile(paste0(outputpathfn,'.segcorvar.txt.gz')),row.names=F,quote=F,sep='\t')
cat(paste0('Analysis finished at ',Sys.time()),sep='\n')
cat(paste0('Results saved to ',outputpathfn,'.segcorvar.txt.gz'),sep='\n')
cat(paste0('Running log can be found in ',outputpathfn,'.log'),sep='\n')
