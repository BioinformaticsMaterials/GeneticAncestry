#23andME raw data for funsies 
#http://www.vincebuffalo.com/2012/03/12/23andme-gwascat.html

#snpStats

source("http://bioconductor.org/biocLite.R")
biocLite("gwascat")
biocLite("ggplot2")
biocLite("TxDb.Hsapiens.UCSC.hg18.knownGene") #Hg18 is what 23andme references! 
biocLite("org.Hs.eg.db")
biocLite("AnnotationDbi")
biocLite("GenomicRanges")
biocLite("ggbio")
biocLite("biovizBase")

library(gwascat)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(IRanges)
library(GenomicRanges)
library(ggbio)
library(biovizBase)
#biomart
#rtracklayer 
#SQLite




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

txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
class(txdb)
transcripts(txdb)

tx.by.gene <- transcriptsBy(txdb, "gene")
tx.by.gene



columns(org.Hs.eg.db)
select(org.Hs.eg.db, keys="APOE", columns=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="SYMBOL")
# SYMBOL ENTREZID         GENENAME
# 1   APOE      348 apolipoprotein E
# A word of caution: Entrez Gene IDs are names and thus they need to be quoted 
# when working with GRangesList objects from transcript databases.

tx.by.gene["348"]
# tx.by.gene[348]
# If I had used tx.by.gene[348] the 348th element of the list would 
# have been returned, not the transcript data for the APOE gene 
# (which has Entrez Gene ID “348”).

#Now, do any SNPs fall in this region? 
#Let’s build a GRanges object from my genotyping data, and look for overlaps. 
#Before I do, it’s worth mentioning another gotcha about working with 
#bioinformatics data: chromosome naming schemes. Different databases use 
#all sorts of schemes, and you should always check them. 23andme returns 
#just numbers, X, Y, and MT. Let’s change it to use the same as the 
#Bioconductor annotation.


levels(dat$chrom) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")
my.snps <- with(dat, GRanges(seqnames=chrom,IRanges(start=position, width=1),rsid=rsid, genotype=genotype)) # this goes into metadata

apoe.i <- findOverlaps(tx.by.gene["348"], my.snps)


#apoe.i is an object of class RangesMatching. Note that had we not matched 
#chromosome names, Bioconductor gives us a nice warning that sequence names 
#don’t match. We could look at the slots of apoe.i but output can be seen 
#with matchMatrix:
  
hits <- matchMatrix(apoe.i)[, "subject"] ## DOESNT WORK function defunct 
hits

my.snps[hits]

### RISK VARIANTS 
#We can use the metadata provided by gwascat to further look for interesting 
#variants in our 23andme data. I would recommend interpreting this data with 
#caution, as summarizing these findings in a single element metadata data 
#frame is hard: there’s definitely lost information.

data(gwrngs19)
gwrngs.emd <- as.data.frame(elementMetadata(gwrngs19))
datm <- merge(dat, gwrngs.emd, by.x="rsid", by.y="SNPs")


risk.alleles <- gsub("[^\\-]*-([ATCG?])", "\\1", datm$Strongest.SNP.Risk.Allele)
i.have.risk <- mapply(function(risk, mine) {
  risk %in% unlist(strsplit(mine, ""))
}, risk.alleles, datm$genotype)
datm$i.have.risk <- i.have.risk


my.risk <- datm[datm$i.have.risk, ]
rel.cols <- c(colnames(dat), "Disease.Trait", "Risk.Allele.Frequency",
              "p.Value", "i.have.risk", "X95..CI..text.")

head(my.risk[order(my.risk$Risk.Allele.Frequency), rel.cols], 1)
datm[which(datm$rsid == "rs2377339"), "Initial.Sample.Size"]
head(my.risk[grep("European", my.risk$Initial.Sample.Size), rel.cols], 30)
### EXPLORE AND THOROUGHLY UNDERSTAND THESE!!!

## make function to easily summarize everyones! 
## update this table? 



#VIS!

#single chrom
#p <- Ideogram(genome="hg18",cytoband=FALSE)

#mult chroms
data(ideoCyto, package = "biovizBase")
p<-autoplot(seqinfo(ideoCyto$hg18), layout = "karyogram")

#with cytoband info 
biovizBase::isIdeogram(ideoCyto$hg18)
autoplot(ideoCyto$hg19, layout = "karyogram", cytoband = TRUE)

(elementMetadata(gwrngs19)$my.genotype <- 
   dat$genotype[(match(elementMetadata(gwrngs19)$SNPs, dat$rsid))])

elementMetadata(gwrngs19)$my.risk <- with(elementMetadata(gwrngs19), 
                                        mapply(function(risk, mine) {
                                          risk %in% unlist(strsplit(mine, ""))
                                        }, gsub("[^\\-]*-([ATCG?])", "\\1", Strongest.SNP.Risk.Allele), my.genotype))

p + geom_hotregion(gwrngs19, aes(color=my.risk))
