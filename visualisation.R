#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
bedFileName <- args[1]
outputFolder <- args[2]
nameOfOutputFile <- args[3]
#setwd("/Users/german/Desktop/Parseq/DELETIONS/pipeline")
#bedFileName <- "4.bed"
distances <- read.table("distance.xls", header=T, stringsAsFactors=F, row.names=1)
nameOfOutputFile = paste(nameOfOutputFile, ".pdf")

bed <- read.table(bedFileName, header=T, stringsAsFactors=F, sep="\t", skip=1)
first_chr <- 0
dir.create(file.path(outputFolder), showWarnings = FALSE)
setwd(outputFolder)
pdf(file=nameOfOutputFile, width=0.1 * nrow(distances))
for (j in seq(from=1, to=ncol(distances), by = 3)) {
  prev=1
  vertical_lines_chromosomes <- c()
  chr_names <- c()
  dist_to_norm <- c()
  dist_to_del <- c()
  dist_to_dup <- c()
  counter_of_existing_amplicons <- 0
  exons <- c()
  current_exon <- 0
  number_of_exons_per_chromosome <- c()
  number_of_exons_current_chromosome <- 1
  for (i in seq(1:nrow(bed))) {
    if ((bed[i,4] %in% row.names(distances))) {
      ampl_name <- bed[i,4]
      counter_of_existing_amplicons = counter_of_existing_amplicons + 1
      if (first_chr != bed[i,1]) {
        first_chr = bed[i,1]
        chr_names <- c(chr_names, bed[i,1])
        vertical_lines_chromosomes <- c(vertical_lines_chromosomes, counter_of_existing_amplicons)
        number_of_exons_per_chromosome <- c(number_of_exons_per_chromosome, number_of_exons_current_chromosome)
      }
      if (current_exon != bed[i,5]) {
        exons <- c(exons, counter_of_existing_amplicons)
        number_of_exons_current_chromosome = number_of_exons_current_chromosome + 1
        current_exon = bed[i,5]
      }
      
      dist_to_norm <- c(dist_to_norm, distances[ampl_name, j])
      dist_to_del <- c(dist_to_del, distances[ampl_name, j + 1])
      dist_to_dup <- c(dist_to_dup, distances[ampl_name, j + 2])      
    }
  }
  number_of_exons_per_chromosome <- c(number_of_exons_per_chromosome, number_of_exons_current_chromosome)
  
  plot(dist_to_norm, ylim=c(min(dist_to_norm,dist_to_del,dist_to_dup) - 0.1, max(dist_to_norm,dist_to_del,dist_to_dup) + 0.1), col="blue",cex=0.7,
       xlab=paste(chr_names, collapse = ', '), ylab="Green:dup, Blue:normal, Red:deletion",pch=19,cex.lab=1.5, cex.axis=1.5)
  title(sub=colnames(distances)[j],)
  points(dist_to_del, col="red",cex=0.7,pch=19)
  points(dist_to_dup, col="green4",cex=0.7,pch=19)
  vertical_lines_chromosomes <- c(vertical_lines_chromosomes, length(dist_to_norm) + 1)
  for (k in seq(1:length(chr_names) - 1)) {
    lines(c(rep("NA",prev - 1), dist_to_norm[seq(prev, vertical_lines_chromosomes[k + 1] - 1)]))
    lines(c(rep("NA",prev - 1), dist_to_del[seq(prev, vertical_lines_chromosomes[k + 1] - 1)]))
    lines(c(rep("NA",prev - 1), dist_to_dup[seq(prev, vertical_lines_chromosomes[k + 1] - 1)]))
    prev=vertical_lines_chromosomes[k + 1]
  }
  abline(h=c(3,-3))
  abline(v=vertical_lines_chromosomes[-c(1, length(vertical_lines_chromosomes))])
  abline(v=exons,lty=3)
  
  prev=1
  for (k in seq(1:length(chr_names) - 1)) {
    if (vertical_lines_chromosomes[k + 1] - prev > 10 ) {
    plot(dist_to_norm[seq(prev, vertical_lines_chromosomes[k + 1] - 1)], ylim=c(min(dist_to_norm,dist_to_del,dist_to_dup) - 0.1, max(dist_to_norm,dist_to_del,dist_to_dup) + 0.1), col="blue",cex=0.7,
         xlab=paste(chr_names[k]), ylab="Green:dup, Blue:normal, Red:deletion",pch=19)
    title(sub=colnames(distances)[j],)
    points(dist_to_del[seq(prev, vertical_lines_chromosomes[k + 1] - 1)], col="red",cex=0.7,pch=19)
    points(dist_to_dup[seq(prev, vertical_lines_chromosomes[k + 1] - 1)], col="green4",cex=0.7,pch=19)
    lines(c(dist_to_norm[seq(prev, vertical_lines_chromosomes[k + 1] - 1)]))
    lines(c(dist_to_del[seq(prev, vertical_lines_chromosomes[k + 1] - 1)]))
    lines(c(dist_to_dup[seq(prev, vertical_lines_chromosomes[k + 1] - 1)]))
    abline(v=exons[seq(number_of_exons_per_chromosome[k], number_of_exons_per_chromosome[k + 1])] - prev + 1, lty=3)
    prev=vertical_lines_chromosomes[k + 1]
    }
    abline(h=c(3,-3),lwd=2,lty=3)
  }
}
