# High throughput data processing

library(dada2)
library(R.utils)
library(data.table)
library(ggplot2)
library(Biostrings)
library(ShortRead)

# 1. Data import and quality controls----
path = "/home/.../fastq" # Lists of reverse sequences and forward sequences
fnFs = sort(list.files(path, pattern="_R1.fq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_R2.fq.gz", full.names = TRUE))

sample.names = sapply(strsplit(basename(fnFs), "_R1|R2.fq.gz"), `[`, 1) # Extract sample names

plotQualityProfile(fnFs[1:5]) # quality profiles of fwd reads
plotQualityProfile(fnRs[1:5]) # quality profiles of rev reads

# 2. Filter and Trim----
       
    ## 2.1. Are the primer sequences in the reads?----

primer_F = "CCTACGGGNGGCWGCAG"
primer_R = "GACTACHVGGGTATCTAATCC"

pattern = ".fq"

primer_search = function(path, pattern, primer_F, primer_R, size, max_mismatch) {
  sub_samples = sample(list.files(path, pattern = pattern, full.names = T), size)
  seqs = sread(readFastq(sub_samples))
  match_F = vmatchPattern(primer_F, seqs, max.mismatch = 2)
  match_R = vmatchPattern(primer_R, seqs, max.mismatch = 2)
  nb_match = sum(lengths(match_F) > 0) + sum(lengths(match_R) > 0)
  print(paste("number of reads including primers = ", nb_match))
}

primer_search1 = primer_search(path = path, pattern = pattern, primer_F, primer_R, size = 5, max_mismatch = 2)

    ## 2.2. Preparation Filter and Trim----

filtFs = file.path(path = getwd(), "filtered", paste0(sample.names, "_R1_filt.fq"))
filtRs = file.path(path = getwd(), "filtered", paste0(sample.names, "_R2_filt.fq"))

names(filtFs) = sample.names
names(filtRs) = sample.names

out = filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  trimLeft = c(25, 25),
  truncLen = c(290, 245),
  maxN = 0,
  maxEE = c(2,2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

    ## 2.3. Checking primer clipping----

pattern2 = "fastq.gz"
primer_search2 = primer_search(path = path2, pattern = pattern2, primer_F, primer_R, 
                               size = 5, max_mismatch = 2) # primers clipping verification

primer_search2 = primer_search(path = "sub_samples/filtered/", primer_F, primer_R, size = 5, max_mismatch = 2) # primers clipping verification

    ## 2.4. Tracking reads filtering----

out_F = data.table(sample = rownames(out_F), out_F)
out_R = data.table(sample = rownames(out_R), out_R)

out = rbind(out_F, out_R)
out[, diff := reads.in - reads.out]

out[grep("R1", sample), group := "R1"]
out[grep("R2", sample), group := "R2"]
out[, name := sub("_R[12]\\.fq$", "", out[,sample])]

ggplot(out, aes(x = name, y = diff, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Diff", fill = "Group")

# 3. Learn the Error Rates----

errF = learnErrors(filtFs, multithread=F, verbose = T)
errR = learnErrors(filtRs, multithread=F, verbose = T)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# 4. Sample Inference----

dadaFs = dada(filtFs, err=errF, multithread=F)
dadaRs = dada(filtRs, err=errF, multithread=F)
dadaFs[[1]] # dada returns a "dada-class" object including several values, explore "dadaFs[[1,2,3,4...]]"
dadaRs[[1]]

# 5. Merge Paired Reads----

mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# 6. Construct sequence table----

seqtab = makeSequenceTable(mergers)

# 7. Remove chimeras----

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
asv_table = data.table(ASV = colnames(seqtab.nochim), Counts = colSums(seqtab.nochim))

# 8. Track reads through the pipeline----

getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample.names

# 9. Assign taxonomy----

taxa = assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", tryRC = TRUE, multithread=TRUE)
taxa.print = taxa # Removing sequence rownames for display only
rownames(taxa.print) = NULL

unqs.mock = seqtab.nochim[grepl("MOCK2", rownames(seqtab.nochim)), ]