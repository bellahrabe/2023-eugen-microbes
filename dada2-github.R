############ DADA2 - Eugen Microbiome workflow ###############################
# # Q1: Do your sequences derive from different runs? If yes -> Create sequence tables separately and merge afterwards.
# # Q2: Do your sequences contain Adapters/Primer? If yes -> truncate sequences accordingly -> change trimLeft in function out()
# # Q3: How much do your runs overlap? How much can you truncate at the end of each run? -> consider to change truncLen in function out()
# # Q4: Do you have single- or paired end reads?

# # This workflow uses dada2; please cite 
# # Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). "DADA2: High-resolution sample inference from Illumina amplicon data." _Nature Methods_, *13*, 581-583. doi: 10.1038/nmeth.3869 (URL: https://doi.org/10.1038/nmeth.3869)
# # DADA2 workflow for small data. This workflow uses unzipped and zipped forward and reverse reads (or only forward or only reverse).
# # For help check out https://benjjneb.github.io/dada2/index.html

# # # calculate length of overlap:
# #  Eugen 18S
FrwdPrmrStart <- 1422  # start of your forward primer, e.g. 515 for 515F
RvrsPrmrStart <- 1797  # start of your reverse primer, e.g. 909 for 909R

# # Eugen 16S
FrwdPrmrStart <- 515  # start of your forward primer, e.g. 515 for 515F
RvrsPrmrStart <- 909  # start of your reverse primer, e.g. 909 for 909R

ReadLength <- 301     # length of your reads
print(paste("Your reads have an overlap of",(ReadLength-((RvrsPrmrStart-FrwdPrmrStart+1)/2))*2,"bases."))

# # dada2 installation instructions: https://benjjneb.github.io/dada2/dada-installation.html
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.17")
# BiocManager::install("DECIPHER")

library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")

# # ==== User Input ==== 
# # !!! In addition to the input provided here, you have to set your filtering parameters in filterandTrim depending on your data.
PairedReads='paired'                          # set to 'single' or 'paired'
ProjectName='2022_SYES_18S_dada2'                     
ProjectName='2022_SYES_16S_e36_260-190_dada2'

PathToDadaInputFiles=file.path("~/SeqData/BaseCalls_Angelis_Run2new_18S")                                      # Select the folder where dada can find your sequences (fastq files)
PathToDadaOutputFiles=file.path("~/dada2",paste(ProjectName,"_dada2",sep=''))  # Select a folder where dada should save your analysis. Inside this folder a folder with the name of your analysis and _dada2 attached to it will be created
TaxonomyFolder=file.path('~/DatabasesForDada2')                                # Select the folder where you have saved the SILVA taxonomy files 

classmeth <- "dada2" # "dada2" use dada2 algorithm (naive Bayesian classifier) or "dec" for the decipher algorithm. Be aware that depending on your choice you'll need different databases
database="silva1381" 

# which database you want to use 
# for dada2 classifying: 
# database = "silva1381" # for silva_nr_v138.1
# database = "silva138" # for silva_nr_v138
# database = "silva132" # for silva_nr_v132
# database = "rdp18" # for RDP_18
# database = "pr2" # for PR2
# database = "unite" # for unite

# for decipher classifiying: 
# database = "silva138" # for SILVA_SSU_r138_2019
# database = "rdp18" # for RDP_v18
# database = "pr2v4" # for PR2_v4
# database = "unite2021" # for UNITE_v2021_May2021

# # If you only have forward or reverse reads leave the other field empty
PatternForwardReads="_R1.fastq" # Change to the Pattern at the end of your forward reads, if e.g. SAMPLENAME_PATTERN.fastq, provide '_PATTERN.fastq'
PatternReverseReads="_R2.fastq" # Change to the Pattern at the end of your reverse reads

MultiThread=F               # Use MultiThread=F (Windows) or MultiThread=TRUE

reversecomplement=T         # for tryRC in taxonomy assignment. Set T if you want to try if the reverse complement better fits. Might be necessary for reverse reads.

Pool=F                      # Set T, F, or "pseudo". Pooling can increase the sensitivity to rare per-sample variants. Needs more memory and computation time. Pseudo-pooling approximates pooling with less time requirement.

# # Eugen 18S
maxerror <- c(2,2) # for maxEE in dada2::filterAndTrim(). Sets the maximum number of "expected errors" allowed in a read. Default c(2,2) for merged reads and c(2) for single reads. " If you want to speed up downstream computation, consider tightening maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails. Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later." (dada2 homepage) 
qualityscore <- 2 # for truncQ in dada2::filterAndTrim(). Sets the quality score. truncQ=2 will truncate the read at the first nucleotide with a quality score of 2 (if there is one). Read will only  thrown away if it is then shorter than truncLen
trim.lefty <- c(21,20)  # If you need to cut off primers or adapters from your reads. the first number is what will be cut off at the start of the forward reads, the second number is what will be cut off at the start of the reverse reads. if you have single reads (e.g. only forward or only reverse) then only use one number

# # Eugen 16S
maxerror <- c(3,6) # for maxEE in dada2::filterAndTrim(). Sets the maximum number of "expected errors" allowed in a read. Default c(2,2) for merged reads and c(2) for single reads. " If you want to speed up downstream computation, consider tightening maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails. Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later." (dada2 homepage) 
qualityscore <- 2 # for truncQ in dada2::filterAndTrim(). Sets the quality score. truncQ=2 will truncate the read at the first nucleotide with a quality score of 2 (if there is one). Read will only  thrown away if it is then shorter than truncLen
trim.lefty <- c(10,10)  # If you need to cut off primers or adapters from your reads. the first number is what will be cut off at the start of the forward reads, the second number is what will be cut off at the start of the reverse reads. if you have single reads (e.g. only forward or only reverse) then only use one number

if(PairedReads=="paired"){
  print(paste("Your merged reads will have a length of",(RvrsPrmrStart-FrwdPrmrStart-sum(trim.lefty)+1) ,"bases."))
}
# # ==== Create analysis directories ====

dir.create(PathToDadaOutputFiles)
dir.create(file.path(PathToDadaOutputFiles,"plots"))
dir.create(file.path(PathToDadaOutputFiles,"filtered"))

workspaceName=paste(ProjectName,"_dada2",".RData",sep="")

# # ==== Reads fastq's ==== 

list.files(PathToDadaInputFiles) # here, all your fastq files should be listed

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
if (PairedReads=='paired') {
  fnFs <- sort(list.files(PathToDadaInputFiles, pattern=PatternForwardReads, full.names = TRUE))  # creates a list containing all your forward reads
  fnRs <- sort(list.files(PathToDadaInputFiles, pattern=PatternReverseReads, full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), PatternForwardReads), `[`, 1)                   # Extracts sample names
  sample.names <- sapply(strsplit(basename(fnRs), PatternReverseReads), `[`, 1)
} else if (PairedReads=='single') {
  fnFs <- sort(list.files(PathToDadaInputFiles, pattern=PatternForwardReads, full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), PatternForwardReads), `[`, 1)
}

# # ==== Plots and Saves Quality profiles ==== 

if (PairedReads=='paired') {
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Quality_Forward_16.pdf",sep="")),width=10,height=8,useDingbats=FALSE)
  print(plotQualityProfile(fnFs[1:16]))
  dev.off()
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Quality_Reverse_16.pdf",sep="")),width=10,height=8,useDingbats=FALSE)
  print(plotQualityProfile(fnRs[1:16]))
  dev.off()
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Quality_Forward_all.pdf",sep="")),width=20,height=20,useDingbats=FALSE)
  print(plotQualityProfile(fnFs))
  dev.off()
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Quality_Reverse_all.pdf",sep="")),width=20,height=20,useDingbats=FALSE)
  print(plotQualityProfile(fnRs))
  dev.off()
} else if (PairedReads=='single') {
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Quality_Forward_16.pdf",sep="")),width=10,height=8,useDingbats=FALSE)
  print(plotQualityProfile(fnFs[1:16]))
  dev.off()
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Quality_Forward_all.pdf",sep="")),width=20,height=20,useDingbats=FALSE)
  print(plotQualityProfile(fnFs))
  dev.off()
}

# # INPUT NEEDED!!!! # # 
# # Your reads should at least have an overlap of 20 + biological length variation nucleotides.
# # the v3v4 region is bimodal with a peak at ~440 and one of ~460 nucleotides. (https://github.com/benjjneb/dada2/issues/857)

# # Eugen 18S
trim.righty <- c(215,200)   # First number gives the length of the forward reads you want to keep (after which base the sequence will be discarded), the second number is the same for the reverse reads. Keep in mind that you will still need to have overlap! For single reads only use one number.

# # Eugen 16S
trim.righty <- c(260,190)   # First number gives the length of the forward reads you want to keep (after which base the sequence will be discarded), the second number is the same for the reverse reads. Keep in mind that you will still need to have overlap! For single reads only use one number.

print(paste("If you cut your reads like that you'll have",as.character((FrwdPrmrStart+trim.righty[1]-1)-(RvrsPrmrStart-trim.righty[2]+1)+1),"bases overlap."))

# # ==== Filters and trims - Your input needed! ==== 

# removePrimers(fnRs,filtRs,'GGACTACHVGGGTWTCTAAT',orient = T)

# Considerations for your own data: The standard filtering parameters are starting points, not set in stone. 
# If you want to speed up downstream computation, consider tightening maxEE. 
# If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails. 
# Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later. 


if (PairedReads=='paired') {
  filtFs <- file.path(PathToDadaOutputFiles, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # Assigns the filenames for the filtered fastq.gz files
  filtRs <- file.path(PathToDadaOutputFiles, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trim.righty,                            # change the parameters according to your sequencees 
                       maxN=0, maxEE=maxerror, truncQ=qualityscore, rm.phix=TRUE,
                       compress=TRUE, multithread=MultiThread,trimLeft = trim.lefty)
} else if (PairedReads=='single') {
  filtFs <- file.path(PathToDadaOutputFiles, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # Assigns the filenames for the filtered fastq.gz files
  names(filtFs) <- sample.names
  out <- filterAndTrim(fnFs, filtFs,  truncLen=trim.righty,
                       maxN=0, maxEE=maxerror, truncQ=qualityscore, rm.phix=TRUE,
                       compress=TRUE, multithread=MultiThread,trimLeft = trim.lefty)
}

out<-as.data.frame(out)
out$ratio <- out$reads.out/out$reads.in

head(out) # shows reads in and reads out for the first 6 samples

message('ratio filtered/input')
summary(out$ratio)

# # ==== Learns error rates ==== 

if (PairedReads=='paired') {
  errF <- learnErrors(filtFs, multithread=MultiThread)
  errR <- learnErrors(filtRs, multithread=MultiThread)
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Error_Forward.pdf",sep="")),width=8,height=8,useDingbats=FALSE)
  print(plotErrors(errF, nominalQ=TRUE))
  dev.off()
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Error_Reverse.pdf",sep="")),width=8,height=8,useDingbats=FALSE)
  print(plotErrors(errR, nominalQ=TRUE))
  dev.off()
} else if (PairedReads=='single') {
  errF <- learnErrors(filtFs, multithread=MultiThread)
  pdf(file.path(PathToDadaOutputFiles,"plots",paste(ProjectName,"_Error_Forward.pdf",sep="")),width=8,height=8,useDingbats=FALSE)
  print(plotErrors(errF, nominalQ=TRUE))
  dev.off()
}

filtFs <- filtFs[file.exists(filtFs)] # in case some samples were excluded previously they will be filtered out here from the list
filtRs <- filtRs[file.exists(filtRs)]
# # ==== Applies Sample Inference algorithm ==== 

getN <- function(x) sum(getUniques(x))

if (PairedReads=='paired') {
  dadaFs <- dada(filtFs, err=errF, multithread=MultiThread,pool=Pool)
  dadaRs <- dada(filtRs, err=errR, multithread=MultiThread,pool=Pool)
  sapply(dadaFs, getN)/out$reads.out
  sapply(dadaRs, getN)/out$reads.out
  message('ratio denoisedF/filtered')
  summary(sapply(dadaFs, getN)/out$reads.out)
  message('ratio denoisedR/filtered')
  summary(sapply(dadaRs, getN)/out$reads.out)
} else if (PairedReads=='single') {
  dadaFs <- dada(filtFs, err=errF, multithread=MultiThread,pool=Pool)
  sapply(dadaFs, getN)/out$reads.out
  message('ratio denoisedF/filtered')
  summary(sapply(dadaFs, getN)/out$reads.out)
}

# # ==== Merges paired reads ==== 

if (PairedReads=='paired') {
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
  head(mergers[[1]]) # Inspect the merger data.frame from the first sample
  message('ratio merged/denoised')
  summary(sapply(mergers, getN)/((sapply(dadaFs, getN)+sapply(dadaRs, getN))/2))
} 

# # ==== Construct sequence table ==== 

if (PairedReads=='paired') {
  seqtab <- makeSequenceTable(mergers)
} else if (PairedReads=='single') {
  seqtab <- makeSequenceTable(dadaFs)
}

dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

saveRDS(seqtab, file.path(PathToDadaOutputFiles,paste(ProjectName,"_seqtab.rds",sep="")))

# # ==== Optional: Filter out non-target-length sequences ==== 
# seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] 

# # ==== Removes chimeras ==== 

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=MultiThread, verbose=TRUE)
dim(seqtab.nochim)
dim(seqtab.nochim)[2]/dim(seqtab)[2]
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim,file.path(PathToDadaOutputFiles,paste(ProjectName,"_seqtab_nochim.rds",sep=""))) 

message('ratio nochim/merged')
summary(rowSums(seqtab.nochim)/sapply(mergers, getN))

# # ==== Track your reads through the pipeline ==== 

if (PairedReads=='paired') {
    track <- out; rownames(track) <- sample.names
    track <- merge(track,as.data.frame(sapply(dadaFs, getN)) ,by="row.names", all=T); rownames(track) <- track[,1]; track <- track[,-1]
    track[,5] <- track[,4]/track[,2]
    track <- merge(track,as.data.frame(sapply(dadaRs, getN)) ,by="row.names", all=T); rownames(track) <- track[,1]; track <- track[,-1]
    track[,7] <- track[,6]/track[,2]
    track <- merge(track,as.data.frame(sapply(mergers, getN)) ,by="row.names", all=T); rownames(track) <- track[,1]; track <- track[,-1]
    track[,9] <- track[,8]/((track[,6]+track[,4])/2)
    track <- merge(track,as.data.frame(rowSums(seqtab.nochim)) ,by="row.names", all=T); rownames(track) <- track[,1]; track <- track[,-1]
    track[,11] <- track[,10]/track[,8]
    track[,12] <- track[,10]/track[,1]
    colnames(track) <- c("input", "filtered","ratio-filtered-input", "denoisedF","ratio-denoisedF-filtered", "denoisedR","ratio-denoisedR-filtered", "merged","ratio-merged-denoised", "nonchim","ratio-nochim-merged",'ratio-out-in')
} else if (PairedReads=='single') {
    track <- out; rownames(track) <- sample.names
    track <- merge(track,as.data.frame(sapply(dadaFs, getN)) ,by="row.names", all=T); rownames(track) <- track[,1]; track <- track[,-1]
    track[,5] <- track[,4]/track[,2]
    track <- merge(track,as.data.frame(rowSums(seqtab.nochim)) ,by="row.names", all=T); rownames(track) <- track[,1]; track <- track[,-1]
    track[,7] <- track[,6]/track[,4]
    track[,8] <- track[,6]/track[,1]
    colnames(track) <- c("input", "filtered","ratio-filtered-input", "denoisedF","ratio-denoisedF-filtered","nonchim","ratio-nochim-denoisedF","ratio-out-in")
} 

message("For the following samples no reads remained:")
rownames(track)[rowSums(is.na(track))>0]

head(track)
message('ratio filtered/input')
summary(track$`ratio-filtered-input`)
if (PairedReads=='paired') {
  message('ratio denoisedF/filtered')
  print(summary(track$`ratio-denoisedF-filtered`))
  message('ratio denoisedR/filtered')
  print(summary(track$`ratio-denoisedR-filtered`))
  message('ratio merged/denoised')
  print(summary(track$`ratio-merged-denoised`))
  message('ratio nochim/merged')
  print(summary(track$`ratio-nochim-merged`))
} else if (PairedReads=='single') {
  message('ratio denoisedF/filtered')
  print(summary(track$`ratio-denoisedF-filtered`))
  message('ratio nochim/denoisedF')
  print(summary(track$`ratio-nochim-denoisedF`))
}
message('ratio out/in')
print(summary(track$`ratio-out-in`))

write.table(track,file.path(PathToDadaOutputFiles,paste(ProjectName,"_track.txt",sep="")))

track.overview <- do.call(cbind, lapply(track, summary))
write.table(track.overview,file.path(PathToDadaOutputFiles,paste(ProjectName,"_track-overview.txt",sep="")))

save.image(file.path(PathToDadaOutputFiles,workspaceName))
gc()

# # ==== Assigns taxonomy ==== 

if (classmeth == "dada2") {
  if (database=='silva132'){
    taxa <- assignTaxonomy(seqtab.nochim, file.path(TaxonomyFolder,"dada","silva_nr_v132_train_set.fa.gz"), multithread=MultiThread,tryRC =reversecomplement) 
    taxa.species <- addSpecies(taxa,file.path(TaxonomyFolder,"dada","silva_species_assignment_v132.fa.gz"),tryRC =reversecomplement)
  } else if (database == 'silva138') {
    taxa <- assignTaxonomy(seqtab.nochim, file.path(TaxonomyFolder,"dada","silva_nr_v138_train_set.fa.gz"), multithread=MultiThread,tryRC =reversecomplement) 
    taxa.species <- addSpecies(taxa,file.path(TaxonomyFolder,"dada","silva_species_assignment_v138.fa.gz"),tryRC =reversecomplement)
  } else if (database == "silva1381") {
    taxa <- assignTaxonomy(seqtab.nochim, file.path(TaxonomyFolder,"dada","silva_nr99_v138.1_train_set.fa.gz"), multithread=MultiThread,tryRC =reversecomplement) 
    taxa.species <- addSpecies(taxa,file.path(TaxonomyFolder,"dada","silva_species_assignment_v138.1.fa.gz"),tryRC =reversecomplement)
  } else if (database == "rdp18") {
    taxa <- assignTaxonomy(seqtab.nochim, file.path(TaxonomyFolder,"dada","rdp_train_set_18.fa.gz"), multithread=MultiThread,tryRC =reversecomplement) 
    taxa.species <- addSpecies(taxa,file.path(TaxonomyFolder,"dada","rdp_species_assignment_18.fa.gz"),tryRC =reversecomplement)
  } else if (database == "unite") {
    taxa.species <- assignTaxonomy(seqtab.nochim, file.path(TaxonomyFolder,"dada","sh_general_release_dynamic_all_16.10.2022_dev.fasta"), multithread=MultiThread,tryRC =reversecomplement)
  } else if (database == "pr2") {
    taxa.species <- assignTaxonomy(seqtab.nochim, file.path(TaxonomyFolder,"dada","pr2_version_5.0.0_SSU_dada2.fasta.gz"), multithread=MultiThread,tryRC =reversecomplement, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))
  }
} else if (classmeth == "dec") {
  dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
  if (database == "silva138") {
    load(file.path(TaxonomyFolder,"decipher","SILVA_SSU_r138_2019.RData")) 
  } else if (database == "pr2v4") {
    load(file.path(TaxonomyFolder,"decipher","PR2_v4_13_March2021.RData")) 
  } else if (database == "unite2021") {
    load(file.path(TaxonomyFolder,"decipher","UNITE_v2021_May2021.RData")) 
  } else if (database == "rdp18") {
    load(file.path(TaxonomyFolder,"decipher","RDP_v18-mod_July2020.RData")) 
  }
  ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
  taxa.species <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa.species <- x$taxon[m]
    taxa.species[startsWith(taxa.species, "unclassified_")] <- NA
    taxa.species
  }))
  colnames(taxa.species) <- ranks; rownames(taxa.species) <- getSequences(seqtab.nochim)
}
gc()

taxa.print <- taxa.species # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

## save taxonomy files
if(classmeth == "dada2") {
  if(database == "unite" | database == "pr2") {
    saveRDS(taxa.species, file.path(PathToDadaOutputFiles,paste(ProjectName,"_",classmeth,"_",database,"_taxa_species.rds",sep="")))
  } else {
    saveRDS(taxa, file.path(PathToDadaOutputFiles,paste(ProjectName,"_",classmeth,"_",database,"_taxa.rds",sep="")))
    saveRDS(taxa.species, file.path(PathToDadaOutputFiles,paste(ProjectName,"_",classmeth,"_",database,"_taxa_species.rds",sep="")))
  }
} else if (classmeth == "dec") {
  saveRDS(taxa.species, file.path(PathToDadaOutputFiles,paste(ProjectName,"_",classmeth,"_",database,"_taxa_species.rds",sep="")))
}

save.image(file.path(PathToDadaOutputFiles,workspaceName))

##############################################################
# ##Alternative 2: Download Silva files to your project folder##
# download.file(url="https://zenodo.org/record/1172783/files/silva_nr_v138_train_set.fa.gz?download=1",destfile=file.path(PathToDadaOutputFiles,"silva_nr_v138_train_set.fa.gz"))
# download.file(url="https://zenodo.org/record/1172783/files/silva_species_assignment_v138.fa.gz?download=1",destfile=file.path(PathToDadaOutputFiles,"silva_species_assignment_v138.fa.gz"))