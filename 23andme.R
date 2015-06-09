#23andME raw data for funsies 
#http://www.vincebuffalo.com/2012/03/12/23andme-gwascat.html

#snpStats

source("http://bioconductor.org/biocLite.R")
biocLite("gwascat")
biocLite("ggplot")
library(gwascat)


# Unfortunately they do so without SNP call confidence data, but in a personal correspondence with a 
# 23andme representative they stated:
# Data reproducibility of our genotyping platforms is estimated at about 99.9%. 
# Average call rate is about 99%. When samples do not meet sufficient call rate 
# thresholds, we repeat the analysis, and/or request a new sample. We do not 
# return data to customers that does not meet our quality thresholds.


dat <- read.table("/Volumes/Macintosh HD/Users/Abbie/Dropbox/23andme/data/genome_Abigail_Groff_Full_20140406135813.txt", sep="\t", header=FALSE,
                colClasses=c("character", "character", "numeric", "character"),
                col.names=c("rsid", "chrom", "position", "genotype"))

tmp <- dat$chrom
dat$chrom = ordered(dat$chrom, levels=c(seq(1, 22), "X", "Y", "MT"))

## It's never a bad idea to check your work
stopifnot(all(as.character(tmp) == as.character(dat$chrom)))

# SNPs by chrom:
ggplot(dat) + geom_bar(aes(chrom))

