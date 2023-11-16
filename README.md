# DADA2

```{r}
#Load required libraries 
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
```

```{r}
#Set the working directory to where your sequence files are.
setwd("/Users/larawarriner/Desktop/Seqs_Assignment3")
path<-"/Users/larawarriner/Desktop/Seqs_Assignment3"
list.files(path)
```

```{r}
#If you have a different name format, change the script but keep the object names the same 
# For these samples the forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq.
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names. Again check the format. This function assumes: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

```{r}
#Inspect the quality of sequences (forward) 
plotQualityProfile(fnFs[1:2])
#The cycle is equivalent of a base pair. Quality is best close to the primer.  
#Mean quality score is green (focus on this)
#If you trim too much there might not be enough overlap for ends to match.
#Quality score of 30 is a good cut off point. 
#Ignore the negative control in terms of trimming the sequence. 
#cut at 250 here for instance.  
```

```{r} 
#Do the same for reverse reads.
plotQualityProfile(fnRs[1:2])
#1:2 means look at the first two samples for both the forward and reverse
#Cut at 190 here for instance.
```

```{r}
#Making a filtered folder for forward and reverse.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

#For forward reads cutting at 250, for reverse reads cutting at 190. This will need to be changed depending on your data. 
#We are throwing out anything with an unknown nucleotide (N)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,190),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
#Read out shows you how many are lost 
#These are stringent filtering parameters

```{r} 
errF <- learnErrors(filtFs, multithread=TRUE)
# 38724750 total bases in 154899 reads from 3 samples will be used for learning error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
# 29330810 total bases in 154899 reads from 3 samples will be used for learning error rates.
plotErrors(errF, nominalQ=TRUE)
#Estimating number of times there was error in sequencing 
```

```{r}
#Assessing unique sequences
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
```

```{r}
#Merge paired reads 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#Use head to check if it did it right.
#abundance is how many copies of the sequence are in the sample 
#Get rid of copies before this step so it's computationally easier
```

```{r}
#Construct an ASV table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Some pipelines remove outliers
```

```{r}
#Removing chimeras. In the sequencing reaction sometimes sequences get fused together that aren't suppose to.
#DADA2 has built in chimera detection.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```

```{r}
#Compare to a training set of known species data
taxa <- assignTaxonomy(seqtab.nochim, "/Users/larawarriner/Desktop/Seqs_Assignment3/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
```

```{r}
#Optional Extention to species level
taxa <- addSpecies(taxa,"/Users/larawarriner/Desktop/Seqs_Assignment3/silva_species_assignment_v138.1.fa.gz")
```

```{r}
taxa.print <- taxa # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

#NOT in protocol. We want to save the data to computer
write.csv(taxa, file="/Users/larawarriner/Desktop/seqs_Assignment3/taxa")
write.csv(seqtab.nochim, file= "/Users/larawarriner/Desktop/seqs_Assignment3/seqtab.nochim.csv")

#Continue from where we left off 
taxa<- read.csv (file="/Users/larawarriner/Desktop/seqs_Assignment3/taxa", sep =',', row.names =1)

seqtab.nochim<- read.csv(file="/Users/larawarriner/Desktop/seqs_Assignment3/seqtab.nochim.csv", sep =',' , row.names=1)

#First convert data frame into a matrix. If we didn't use first row as a header, would be difficult to make it into a matrix. Bold is not data it's a header. 
seqtab.nochim<-as.matrix(seqtab.nochim)
taxa<-as.matrix(taxa)

```{r}
#It is useful to have taxonomic identity and the abundance. First step transpose (ie., flip) the seqtab.nochim.
flipped_seqtab.nochim<- as.data.frame(t(seqtab.nochim))
dim(flipped_seqtab.nochim)
```

```{r}
#Now we will merge the files
OTUabund<- cbind(flipped_seqtab.nochim, taxa)
#Look at sheet to find the shared ASVs. 
```

write.csv(OTUabund, file="/Users/larawarriner/Desktop/seqs_Assignment3/OTUabund.csv")

```{r}
#Now lets make some graphs, First lets make a data frame
#Take row names from seqtab (names of samples)
samples.out<- rownames(seqtab.nochim)
#Make data frame of sample names
samdf<-data.frame(samples.out)
rownames(samdf) <- samples.out
```

```{r}
#Hand off to phyloseq. Perform phyloseq, abundance table (otu_table) is in this file. We don't have taxa. 
#Then saying sample data is in an object called sample data frame. Where we actually have the information is in taxa. It is just called OTU because that is what used to be used. 
ps<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf),tax_table(taxa))
ps
```

```{r}
dna<-Biostrings::DNAStringSet(taxa_names(ps))
names(dna)<-taxa_names(ps)
ps<- merge_phyloseq(ps, dna)
taxa_names(ps)<-paste0("ASV", seq(ntaxa(ps)))
ps
```

```{r}
#If we want to use ggplot alone we need to get data out of phyloseq into a dataframe.Use the psmelt function.
ps.table<-psmelt(ps)
#Need to factor the column. We are focusing on phylym so we will factor it.
ps.table$Phylum<-factor(ps.table$Phylum)
ggplot(ps.table,mapping =  aes(x=Sample, y=Abundance))+ geom_bar(stat = "identity", position="stack") + aes(fill=Phylum)
#If making geom bar need stat= identity if making the x and y.
```  

```{r}
#Often presented as relative abudance (Example: what percentage is Baceriodata). Easy way to compare between samples. Percentage of the total basically. 
ps.prop<-transform_sample_counts(ps,function(OTU) OTU/sum(OTU))
#Transform into dataframe 
ps.proptable<-psmelt(ps.prop)
ps.proptable$Phylum<-factor(ps.proptable$Phylum)
```

```{r}
#now graph it
ggplot(ps.proptable,mapping =  aes(x=Sample, y=Abundance))+ geom_bar(stat = "identity", position="stack") + aes(fill=Phylum) +
  xlab("Sample") +
  ylab("Relative Abundance")                 +
  ggtitle("Relative Abundance of Phyla in Permafrost Samples from the Canadian Arctic")+ scale_x_discrete(labels=c("Hummock", "Negative Control", "Trough"))
```

```{r}
#now graph it for order
ps.proptable$Order<-factor(ps.proptable$Order)
ggplot(ps.proptable,mapping =  aes(x=Sample, y=Abundance))+ geom_bar(stat = "identity", position="stack") + aes(fill=Order) +
  xlab("Sample") +
  ylab("Relative Abundance")                 +
  ggtitle("Relative Abundance of Order in Permafrost Samples from the Canadian Arctic")+ scale_x_discrete(labels=c("Hummock", "Negative Control", "Trough"))
```

```{r}
ggplot(ps.proptable,mapping =  aes(x=Sample, y=Phylum))+ geom_point (aes(size=Abundance*100)) +aes(color=Phylum) + guides(color="none") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+   xlab("Sample") +
  ylab("Phylum") + ggtitle("Relative Abundance of Phyla in Permafrost Samples from the Canadian Arctic")+ scale_x_discrete(labels=c("Hummock", "Negative Control", "Trough"))+ labs(size='Relative Abundance (%)') 
```

