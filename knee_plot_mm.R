#!/usr/bin/env Rscript

##Mike Mariani, Frietze Lab, UVM, 2020

##This is based off the dropseq cookbook from 
##BROAD/McCarroll lab

##This should produce a "cumulative distibution plot"
##Which I believe is the same thing as a "knee plot"
##Like we have looked at in cellranger qc htmls from
##the 10x platform. 

##Original R code from the manual:
##a=read.table("100cells_numReads_perCell_XC_mq_10.txt.gz", header=F, stringsAsFactors=F) 
##x=cumsum(a$V1) 
##x=x/max(x) 
##plot(1:length(x), 
##     x, 
##     type='l', 
##     col="blue", 
##     xlab="cell barcodes sorted by number of reads [descending]", 
##     ylab="cumulative fraction of reads", xlim=c(1,500)) 

library(cowplot)
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)

knee.plot <- function(working.dir, height, width){
  library(cowplot)
  library(dplyr)
  library(magrittr)
  library(data.table)
  library(ggplot2)
  files.in <- list.files(path=working.dir,
                         pattern="_readcounts.txt",
                         full.names=TRUE)
  ##print(files.in)
  ##length(files.in)
  for(z in 1:length(files.in)){
        title <- paste0(gsub(".star_cell_readcounts.txt",".kneeplot",basename(files.in[z])),".pdf")
        frame.now <- read.table(files.in[z], header=F, stringsAsFactors = F)
        ##print(nrow(frame.now))
        ##print(colnames(frame.now))
        x=cumsum(frame.now$V1) 
        x=x/max(x)
        plot.out <- ggplot() +
          geom_line(aes(x=1:length(x), y=x %>% sort(decreasing = TRUE)), color="blue") +
          xlab("cell barcodes sorted by number of reads [descending]") + 
          ylab("cumulative fraction of reads") +
          theme_bw() +
          ggtitle(title) +
          theme(plot.title=element_text(hjust=0.5))
        ggsave(filename=paste0(working.dir,"/",title),
               height = height, 
               width = width,
               device="pdf",
               plot=plot.out)
  }
}

working.dir.1 <- "/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock"
working.dir.2 <- "/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/wt"

plots.out.1 <- knee.plot(working.dir=working.dir.1, 
                         height=8, 
                         width=8)
plots.out.2 <- knee.plot(working.dir=working.dir.2, 
                         height=8, 
                         width=8)
