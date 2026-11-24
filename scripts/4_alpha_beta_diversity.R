# Alpha and beta diversity analysis

library(data.table)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(ggbeeswarm)
library(ggsci)
library(dplyr)
library(vegan)
library(MeanRarity)
library(hillR)
library(mgcv)
library(DHARMa)

  # 1. Overview of sequencing depth----

source("/.../2_seq_depth.R") # print a figure "sequencing depth by sample"

  # 2. Import output data from dada2 and production of a table----

setwd("/.../datasets/")

seqtab.nochim = fread("4_1_seqtab.nochim.AVUL.txt") # ASV file from dada2
taxa = fread("4_2_taxa.txt") # taxonomy file from dada2

seqs = colnames(seqtab.nochim[, -"names"]) # assign an ID to each sequence

tax_dt = as.data.table(taxa)
tax_dt[, DNA_SEQ := seqs]
tax_dt[, ASV := seq_len(.N)]

counts_dt = as.data.table(seqtab.nochim, keep.rownames = "names")

counts_dt = melt(
  counts_dt,
  id.vars = "names",
  variable.name = "DNA_SEQ",
  value.name = "counts",
  variable.factor = FALSE
) # Transition from a wide format to a long format (one line = combinaison sample x ASV)
setnames(counts_dt, "names", "samples")

counts_dt[, DNA_SEQ := as.character(DNA_SEQ)]

dt = merge(
  counts_dt,
  tax_dt,
  by = "DNA_SEQ",
  all.x = TRUE,
  sort = FALSE
) # merging counts and taxonomy

    ## 2.1. Data set overview----
    
uniqueN(dt$ASV) # total number of unique ASV across all samples
dt[, length(unique(Phylum))] # number of phylums
dt[, length(unique(Class))] # number of classes
dt[, length(unique(Order))] # number of orders
dt[, length(unique(Family))] # number of families
dt[, length(unique(Genus))] # number of genera

#If needed for exploration: 
# dt2 = dt[grep("AVUL", samples),] # to calculate the proportion of reads removed due to the filter with blanks
# dt2[Genus %in% filter, sum(counts)]/dt2[, sum(counts)] * 100 # to calculate the proportion of reads removed due to the filter with blanks

# dt_filt = dt[!Genus %in% filter,] # number of ASV after filtering, see below for filters
# dt_filt = dt_filt[!Genus == "", ]
# dt_filt = dt_filt[grep("AVUL", samples), ]
# dt_filt[counts > 0, length(samples)]

# dt_genus[!Genus == "" & total_reads > 0, length(samples)] # As there are multiple lines per sample (one line = one ASV), it is necessary to calculate the number of unassigned samples using this method

    ## 2.2. Mocks samples : DNA extraction, PCR and sequencing biases----

mocks = dt[grep("MOCK", samples)]
mocks[, plate := sub(".*(P[1-2]).*", "\\1", samples)] # column plate : P1, P2

mocks = mocks[, .(total_reads = sum(counts)), by = .(Genus, samples)]
mocks = mocks[!Genus == "", ] # remove unassigned ASVs
mocks = mocks[total_reads > 0 ,]
mocks = mocks[order(total_reads, decreasing = T)] # order by number of reads

mocks[, plate := sub(".*(P[1-2]).*", "\\1", samples)]

mocks = mocks[, .(reads = sum(total_reads)), by = .(samples, Genus, plate)]
mocks = mocks[, prop := round(reads/sum(reads), 3), by = .(samples)]
mocks[, sum(prop), by = Genus] # check that prop = 1 for all mocks

ggplot(mocks[order(-prop), head(.SD, 8), by = samples], aes(x = plate, y = prop, fill = Genus)) +
  geom_bar(stat = "identity") + 
  scale_fill_rickandmorty() +
  theme_bw() + xlab("Plate") +
  ylab("Proportions of reads assigned")

mocks_expect = fread("4_3_mocks_expected_prop.txt")
mocks_expect[, sum(Prop_expected)]
mocks_expect = mocks_expect[!Genus %in% c("Saccharomyces", "Cryptococcus"), ] # among the 10 microorganisms genera in the mock, 2 are yeast genera -> exclusion
mocks_expect[, Prop_expected := (Prop_expected*100)/96] # Recalculation of proportions
mocks_expect[, sum(Prop_expected)]

mocks_obs = mocks[order(-prop), head(.SD, 8), by = samples] # selection of eight most bacterial genera (eight are expected, see above)
mocks_obs = mocks_obs[, .(prop = sum(prop)), by = .(Genus)]
mocks_obs = mocks_obs[, Prop_obs := (prop * 100)/sum(prop)] # normalization
mocks_obs[, plate := "All plates"] # combine the proportions for the two plates, as there does not appear to be any noticeable difference between them
mocks_obs[, sum(Prop_obs)]
mocks_obs[Genus == "Escherichia/Shigella", Genus := "Escherichia"]

mocks_obs[, prop := NULL]
all_mocks = mocks_obs[, Prop_expected := mocks_expect$Prop_expected]
melt_mocks = melt(all_mocks, id.vars = "Genus", measure.vars = c("Prop_expected", "Prop_obs"), variable.name = "Type", value.name = "Prop")
melt_mocks[, Prop := Prop/100]

ggplot(melt_mocks, aes(x = Type, y = Prop, fill = Genus)) + 
  geom_bar(stat = "identity") + 
  scale_fill_rickandmorty(labels = c(
    "<i>Bacillus</i>",
    "<i>Enterococcus</i>",
    "<i>Escherichia</i>",
    "<i>Lactobacillus</i>",
    "<i>Listeria</i>",
    "<i>Pseudomonas</i>",
    "<i>Salmonella</i>",
    "<i>Staphylococcus</i>")) + 
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, 0.25),
                     labels = function(x) gsub("\\.00$", "", sprintf("%.2f", x))) +
  scale_x_discrete(breaks = c("Prop_obs", "Prop_expected"), labels = c("Observed", "Expected"), ) +
  xlab("") +
  ylab("Relative abundance") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.text = ggtext::element_markdown())

  # 3. Table transformations, filters and matrix----

dt_genus = dt[, .(total_reads = sum(counts)), by = .(Genus, samples)] # extraction of a dataset with total reads by genus x samples combinaison

    ## 3.1. Removing ASV assigned to "Chloroplast" or "Mitochondria"----

dt = dt[Order != "Chloroplast",]
dt = dt[Family != "Mitochondria",]

dt_genus = dt[, .(total_reads = sum(counts)), by = .(Genus, samples)]
# dt_genus[Genus == "" & total_reads > 0, length(samples)] # number of unassinged ASV at genus level

    ## 3.2. Filtering based on bacteria detected in blanks----

blanks = dt_genus[grep("Ctrl_ADN|Ctrl_PCR", samples), ] # Extract blanks
blanks[, plate := sub(".*(P[1-2]).*", "\\1", samples)] # column plate : P1, P2
blanks[, type := sub(".*(ADN|PCR).*", "\\1", samples)] # type of controls : DNA/PCR
blanks[type == "ADN", type := "DNA"]

# blanks[Genus == "" & total_reads > 0, length(samples)] # number of unassigned ASV at genus level
blanks = blanks[!Genus == "", ] # remove unassigned ASVs
blanks = blanks[total_reads > 0 ,]
blanks = blanks[order(total_reads, decreasing = T)] # order by number of reads
blanks_sum = blanks[ , .(total_reads = sum(total_reads, na.rm = TRUE)), by = .(samples, Genus)] # Sum the reads counts for each sample–Genus pair
tab = data.table::dcast(blanks_sum, samples ~ Genus, value.var = "total_reads", fill = 0)

mat = as.matrix(tab[, -1, with = FALSE])
rownames(mat) = tab$samples
storage.mode(mat) = "numeric"

rarecurve(
  mat,
  step = 10,
  label = FALSE,
  col = "steelblue",
  lwd = 1,
  xlab = "Number of reads",
  ylab = "Observed richness (number of genera)",
  main = "Rarefaction curve by sample"
) # rarefaction curves
abline(v = 9000, col = "red", lty = 2, lwd = 2)

depths = rowSums(mat) # total number of reads per samples
summary(depths)

thresholds = c(9000, 10000, 20000)

data.frame(
  depth_target = thresholds,
  n_samples = sapply(thresholds, function(d) sum(depths >= d)),
  perc_conserved = round(
    sapply(thresholds, function(d) mean(depths >= d) * 100),
    1
  )
) # function estimating the percentage of conserved samples depending on depth_target

depth_target = 9000 # Select a depth of x reads based on summary(depths)

mat_filt = mat[depths >= depth_target, , drop = FALSE] # Keep samples with enough reads

set.seed(123)
mat_rare = rrarefy(mat_filt, sample = depth_target) # rarefaction

summary(rowSums(mat_rare)) # number of reads counts after rarefaction
nrow(mat_rare) # number of conserved samples
all(rowSums(mat_rare) == depth_target) # verify that each sample has the target depth

blanks_rare = data.table(Genus = colnames(mat_rare), total_reads = colSums(mat_rare)) # total abundance in each genus in blanks
blanks_rare[grep("Burkhol", Genus), Genus := "B-C-P"] # abbreviation for "Burkholderia-Caballeronia-Paraburkholderia"
setorder(blanks_rare, -total_reads) # order by decreasing abundance
blanks_rare = blanks_rare[, prop := round(total_reads/sum(total_reads, na.rm = TRUE),3)]

ggplot(blanks_rare[prop > 0.01,], aes(x = reorder(Genus, prop), y = prop)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.1, linetype = "dashed",
             colour = "black",
             linewidth = 1) +
  coord_flip() +
  labs(
    x = "Genus",
    y = "Relative abundance"
  ) +
  theme_bw()

filter = blanks_rare[prop > 0.1, Genus] # filter based on blanks : remove genus with a frequency > 10%

dt_genus2 = dt_genus[grep("AVUL", samples), ] # keep only samples and exclude negative controls
dt_genus2[, freq := total_reads / sum(total_reads), by = samples] # verification that the genera included in this filter are not predominant in the samples

# dt_genus_filt = dt_genus2[!Genus %in% filter, ] # creation of a filtered data set but with Wolbachia for analysis of reads assigned to Wolbachia (see wolbachia_load.R)
# fwrite(dt_genus_filt, "dt_genus_filt_with_wolbachia.txt")

top20_genus = dt_genus2[ , .(freq_total = sum(freq, na.rm = TRUE)), by = Genus][order(-freq_total)][1:20]
dt_genus2 = dt_genus[!Genus %in% c("", "Wolbachia"), ]
dt_genus2[, freq := total_reads / sum(total_reads), by = samples] # same but without unassigned ASV at genus level and Wolbachia
top20_genus = dt_genus2[ , .(freq_total = sum(freq, na.rm = TRUE)), by = Genus][order(-freq_total)][1:20]

dt_genus_filt = dt_genus2[!Genus %in% filter, ]
dt_genus_filt[, uniqueN(Genus)] # number of genera after curation, 1246 + 1 if Wolbachia is considered
dt_genus_filt[, freq := total_reads / sum(total_reads), by = samples] 
top20_genus = dt_genus_filt[ , .(freq_total = sum(freq, na.rm = TRUE)), by = Genus][order(-freq_total)][1:20]

    ## 3.3. Preparation of a dataset for statistical analyses----

dt_genus_filt = dt_genus_filt[grep("AVUL", samples), ]

dt_genus_filt[grepl("AVUL", samples), timepoint := sub(".*(T[1-5]).*", "\\1", samples)] # add a column "timepoint" (T1, T2, T3, T4, T5)
dt_genus_filt[grepl("AVUL", samples), condition := sub(".*(ctrl|test).*", "\\1", samples)] # add a column "condition" (ctrl/test)
dt_genus_filt[, samples := gsub(pattern = "*_F_filt.fastq.gz", replacement = "", x = samples)]
dt_genus_filt = dt_genus_filt[Genus != "",] # Removes the unassigned ASV because they are not informative for alpha diversity
dt_genus_filt = dt_genus_filt[total_reads > 0] # 0 are useless for Shannon index estimation

      ### 3.3.1. Production of a matrix (rows : samples, columns : genera, cells : reads counts)----

dt_sum = dt_genus_filt[ , .(total_reads = sum(total_reads, na.rm = TRUE)), by = .(samples, Genus)] # Sum the reads counts for each sample–Genus pair

tab = data.table::dcast(dt_sum, samples ~ Genus, value.var = "total_reads", fill = 0)

mat = as.matrix(tab[, -1, with = FALSE])
rownames(mat) = tab$samples
storage.mode(mat) = "numeric"
    
      ### 3.3.2. Dataset rarefaction----

rarecurve(
  mat,
  step = 10,
  label = FALSE,
  col = "steelblue",
  lwd = 1,
  xlab = "Number of reads",
  ylab = "Observed richness (number of genera)",
  main = "Rarefaction curve by sample"
)
abline(v = 35000, col = "black", lty = 2, lwd = 2)

depths = rowSums(mat) # total number of reads per samples
summary(depths)

thresholds = c(10000, 25000, 30000, 35000, 40000, 50000)

data.frame(
  depth_target = thresholds,
  n_samples = sapply(thresholds, function(d) sum(depths >= d)),
  perc_conserved = round(
    sapply(thresholds, function(d) mean(depths >= d) * 100),
    1
  )
)

depth_target = 35000 # Select a depth of x reads based on summary(depths) and rarefaction curves

mat_filt = mat[depths >= depth_target, , drop = FALSE] # Keep samples with enough reads

set.seed(123)
mat_rare = rrarefy(mat_filt, sample = depth_target)

summary(rowSums(mat_rare)) # number of reads counts after rarefaction
nrow(mat_rare) # number of conserved samples
ncol(mat_rare) # number of conserved genus
all(rowSums(mat_rare) == depth_target) # Verify that each sample has the target depth

  # 4. Alpha diversity : Hill numbers (q0: Richness, q1: Exponentiated Shannon, q2: Inverse Simpson)----

hill_rich = hill_taxa(mat_rare, q = 0) # calculate richness
hill_sh = hill_taxa(mat_rare, q = 1) # calculate exponentiaed Shannon
hill_simps = hill_taxa(mat_rare, q = 2) # calculate inverse Simpson

hill_dt = function(res, value_name) {
  dt = data.table(
    sample = names(res),
    value = as.numeric(res)
  )
  setnames(dt, "value", value_name)
  dt[, timepoint := sub(".*(T[1-5]).*", "\\1", sample)]
  # dt[, treatment := sub(".*(ctrl|test).*", "\\1", sample)]
  dt[, treatment := ifelse(grepl("ctrl", sample), "Control", "Exposed")]
  dt[, sample := gsub("_F_filt\\.fastq\\.gz$", "", sample)]
  dt[, `:=`(
    timepoint = factor(timepoint, levels = paste0("T", 1:5)),
    treatment = factor(treatment)
  )]
  return(dt)
} # function to transform matrices into dt with metadata (treatments, timepoint etc.)

hill_rich = hill_dt(hill_rich, "hill_rich")
hill_sh = hill_dt(hill_sh, "hill_sh")
hill_simps = hill_dt(hill_simps, "hill_simps")

all0 = left_join(hill_rich, hill_sh, by = c("sample", "timepoint", "treatment"))
all = left_join(all0, hill_simps, by = c("sample", "timepoint", "treatment"))
setcolorder(all, c("sample", "timepoint", "treatment", "hill_rich", "hill_sh", "hill_simps"))
all[, timepoint_num := as.numeric(gsub(".*([1-5]).*", "\\1", timepoint))]
data.table::set(all, j = "timepoint_num", value = as.integer(sub(".*([1-5]).*", "\\1", all$timepoint)))
str(all)

levels_x = sort(unique(all$timepoint_num))

all_melt = melt(all, id.vars = c("sample", "timepoint", "timepoint_num", "treatment"), measure.vars = c("hill_rich", "hill_sh", "hill_simps"), variable.name = "index", value.name = "value")

hill_rich1 = ggplot() +
  geom_boxplot(
    data = all,
    aes(x = as.factor(timepoint_num), y = hill_rich, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data =  all,
    aes(x = timepoint_num, y = hill_rich, color = treatment),
    # alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "#cd4f39",
             linewidth = 1) +
  scale_x_discrete(
    breaks = levels_x,
    labels = c(
      "Before HW",
      "End of the HW",
      "16 days post HW",
      "32 days post HW",
      "64 days post HW"
    )
  ) +
  scale_color_manual(values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  scale_fill_manual(values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  labs(
    x = "Sampling point",
    y = "Richness") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("A")

hill_sh2 = ggplot() +
  geom_boxplot(
    data = all,
    aes(x = as.factor(timepoint_num), y = hill_sh, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = all,
    aes(x = timepoint_num, y = hill_sh, color = treatment),
    # alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "#cd4f39",
             linewidth = 1) +
  scale_x_discrete(
    breaks = levels_x,
    labels = c(
      "Before HW",
      "End of the HW",
      "16 days post HW",
      "32 days post HW",
      "64 days post HW"
    )
  ) +
  scale_color_manual(values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  scale_fill_manual(values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  labs(
    x = "Sampling point",
    y = "Exponentiated Shannon") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("B")

hill_simps3 = ggplot() +
  geom_boxplot(
    data =  all,
    aes(x = as.factor(timepoint_num), y = hill_simps, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data =  all,
    aes(x = timepoint_num, y = hill_simps, color = treatment),
    # alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "#cd4f39",
             linewidth = 1) +
  scale_x_discrete(
    breaks = levels_x,
    labels = c(
      "Before HW",
      "End of the HW",
      "16 days post HW",
      "32 days post HW",
      "64 days post HW"
    )
  ) +
  scale_color_manual(values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  scale_fill_manual(values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  labs(
    x = "Sampling point",
    y = "Inverse Simpson") + 
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("C")

hill_rich1|hill_sh2|hill_simps3

    ## 4.1. Analyse 1: Has the heatwave exposure affected gut bacterial community?----

hist(all$hill_rich)
hist(all$hill_sh)
hist(all$hill_simps)

gam_mod_hill_rich1 = gam(hill_rich ~ treatment + treatment : timepoint_num + s(timepoint_num, k = 5),
                         data = all, method = "REML", family = gaussian) # GAM with interaction treatment x experiment time

summary(gam_mod_hill_rich1)
anova(gam_mod_hill_rich1) # check whether the interaction is significant
gam.check(gam_mod_hill_rich1)
qq.gam(gam_mod_hill_rich1)

gam_mod_hill_rich2 = gam(hill_rich ~ treatment + s(timepoint_num, k = 5), data = all, method = "REML", family = gaussian) # Test the effect of timepoint_num (non-linear smooth effect) on hill_rich, while controlling for treatment

summary(gam_mod_hill_rich2)
anova(gam_mod_hill_rich2) # no significant effect of s(timepoint_num) -> linear model LM
gam.check(gam_mod_hill_rich2)
qq.gam(gam_mod_hill_rich2)

lm1 = lm(hill_rich ~ treatment + treatment * timepoint_num, data = all)  # LM with interaction treatment x experiment time
summary(lm1)
anova(lm1) # no significant interaction -> LM without interaction (lm2)
plot(simulateResiduals(lm1))

lm2 = lm(hill_rich ~ treatment + timepoint_num, data = all)
summary(lm2)
anova(lm2)
plot(simulateResiduals(lm2))

lm3 = lm(hill_rich ~ treatment, data = all[timepoint_num<2,]) # test whether richness was significantly different prior to the heatwave
summary(lm3)
anova(lm3)

gam_mod_sh1 = gam(hill_sh ~ treatment + treatment : timepoint_num + s(timepoint_num, k = 5),
                  data = all, method = "REML", family = gaussian)  # GAM with interaction treatment x experiment time

summary(gam_mod_sh1)
anova(gam_mod_sh1) # check whether the interaction is significant
gam.check(gam_mod_sh1)
qq.gam(gam_mod_sh1)

gam_mod_sh2 = gam(hill_sh ~ treatment + s(timepoint_num, k = 5), data = all, method = "REML", family = gaussian) # Test the effect of timepoint_num (non-linear smooth effect) on hill_sh, while controlling for treatment

summary(gam_mod_sh2)
anova(gam_mod_sh2) # significant non-linear smooth effect -> different temporal patterns by treatment? (gam_mod_sh3)
gam.check(gam_mod_sh2)
qq.gam(gam_mod_sh2)

anova(gam_mod_sh2, gam_mod_sh3, test = "F") # allowing treatment-specific smooth terms significantly improved the model

gam_mod_sh3 = gam(hill_sh ~ treatment + s(timepoint_num, by = treatment, k = 5), data = all, method = "REML", family = gaussian) # different temporal patterns by treatment?

summary(gam_mod_sh3)
anova(gam_mod_sh3) # different temporal patterns by treatment
gam.check(gam_mod_sh3)
qq.gam(gam_mod_sh3)

lm1 = lm(hill_sh ~ treatment, data = all[timepoint_num<2,]) # test whether exponentiated Shannon was significantly different between treatment prior to the heatwave
summary(lm1)
anova(lm1)

gam_mod_simps1 = gam(hill_simps ~ treatment + treatment : timepoint_num + s(timepoint_num, k = 5),
                    data = all, method = "REML", family = gaussian) # GAM with interaction treatment x experiment time

summary(gam_mod_simps1)
anova(gam_mod_simps1) # check whether the interaction is significant
gam.check(gam_mod_simps1)
qq.gam(gam_mod_simps1)

gam_mod_simps2 = gam(hill_simps ~ treatment + s(timepoint_num, k = 5), data = all, method = "REML", family = gaussian) # Test the effect of timepoint_num (non-linear smooth effect) on hill_simps, while controlling for treatment

summary(gam_mod_simps2)
anova(gam_mod_simps2)
gam.check(gam_mod_simps2)
qq.gam(gam_mod_simps2)

gam_mod_simps3 = gam(hill_simps ~ treatment + s(timepoint_num, by = treatment, k = 5), data = all, method = "REML", family = gaussian) # different temporal patterns between treatments?
summary(gam_mod_simps3)
anova(gam_mod_simps3) # different temporal patterns between treatments
gam.check(gam_mod_simps3)
qq.gam(gam_mod_simps3)

anova(gam_mod_simps2, gam_mod_simps3, test = "F") # allowing treatment-specific smooth terms significantly improved the model

    ## 4.2. Plotting GAM model predictions----

levels_x = sort(unique(all$timepoint_num))

hill_rich1 = ggplot() +
  geom_boxplot(
    data = all,
    aes(x = as.factor(timepoint_num), y = hill_rich, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = all,
    aes(x = timepoint_num, y = hill_rich, color = treatment),
    # alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "#cd4f39",
             linewidth = 1) +
  stat_summary(data = all, 
               mapping = aes(x = timepoint_num, y = hill_rich, group = treatment, fill = treatment),
               fun = mean,
               geom = "point",
               shape = 23,
               size = 3,
               position = position_dodge(width = 0.8)
  ) +
  scale_x_discrete(
    breaks = levels_x,
    labels = c(
      "Before HW",
      "End of the HW",
      "16 days post HW",
      "32 days post HW",
      "64 days post HW"
    )
  ) +
  scale_color_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  scale_fill_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  # scale_y_log10() +
  labs(
    x = "Sampling point",
    y = expression(paste("Richness"))) +
  theme_bw() +
  theme(legend.position = "none",
        #       legend.text = element_text(size = 14),
        #       legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("A")

newdata = expand.grid(
  timepoint_num = seq(min(all$timepoint_num ), max(all$timepoint_num), length.out = 200),
  treatment = unique(all$treatment)
)
newdata$pred = predict(gam_mod_sh3, newdata, se.fit = TRUE)$fit
newdata$se = predict(gam_mod_sh3, newdata, se.fit = TRUE)$se.fit

pred_link = predict(gam_mod_sh3, newdata = newdata, type = "link", se.fit = TRUE)

newdata$fit = gam_mod_sh3$family$linkinv(pred_link$fit)
newdata$lower = gam_mod_sh3$family$linkinv(pred_link$fit - 1.96 * pred_link$se.fit)
newdata$upper = gam_mod_sh3$family$linkinv(pred_link$fit + 1.96 * pred_link$se.fit)

levels_x = sort(unique(all$timepoint_num))

hill_sh2 = ggplot() +
  geom_boxplot(
    data = all,
    aes(x = as.factor(timepoint_num), y = hill_sh, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = all,
    aes(x = timepoint_num, y = hill_sh, color = treatment),
    # alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_line(
    data = newdata,
    aes(x = timepoint_num, y = fit,
        colour = treatment,
        group = treatment),
    linewidth = 1.2
  ) +
  geom_ribbon(
    data = newdata,
    aes(x = timepoint_num,
        ymin = lower,
        ymax = upper,
        fill = treatment,
        group = treatment),
    alpha = 0.15,
    colour = NA
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "#cd4f39",
             linewidth = 1) +
  stat_summary(data = all, 
               mapping = aes(x = timepoint_num, y = hill_sh, group = treatment, fill = treatment),
               fun = mean,
               geom = "point",
               shape = 23,
               size = 3,
               position = position_dodge(width = 0.8)
  ) +
  scale_x_discrete(
    breaks = levels_x,
    labels = c(
      "Before HW",
      "End of the HW",
      "16 days post HW",
      "32 days post HW",
      "64 days post HW"
    )
  ) +
  scale_color_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  scale_fill_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  # scale_y_log10() +
  labs(
    x = "Sampling point",
    y = expression(paste("Exponentiated Shannon"))) +
  theme_bw() +
  theme(legend.position = "none",
        #       legend.text = element_text(size = 14),
        #       legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("B")

newdata = expand.grid(
  timepoint_num = seq(min(all$timepoint_num ), max(all$timepoint_num), length.out = 200),
  treatment = unique(all$treatment)
)
newdata$pred = predict(gam_mod_simps3, newdata, se.fit = TRUE)$fit
newdata$se = predict(gam_mod_simps3, newdata, se.fit = TRUE)$se.fit

pred_link = predict(gam_mod_simps3, newdata = newdata, type = "link", se.fit = TRUE)

newdata$fit = gam_mod_simps3$family$linkinv(pred_link$fit)
newdata$lower = gam_mod_simps3$family$linkinv(pred_link$fit - 1.96 * pred_link$se.fit)
newdata$upper = gam_mod_simps3$family$linkinv(pred_link$fit + 1.96 * pred_link$se.fit)

levels_x = sort(unique(all$timepoint_num))

hill_simps3 = ggplot() +
  geom_boxplot(
    data = all,
    aes(x = as.factor(timepoint_num), y = hill_simps, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = all,
    aes(x = timepoint_num, y = hill_simps, color = treatment),
    # alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_line(
    data = newdata,
    aes(x = timepoint_num, y = fit,
        colour = treatment,
        group = treatment),
    linewidth = 1.2
  ) +
  geom_ribbon(
    data = newdata,
    aes(x = timepoint_num,
        ymin = lower,
        ymax = upper,
        fill = treatment,
        group = treatment),
    alpha = 0.15,
    colour = NA
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "#cd4f39",
             linewidth = 1) +
  scale_x_discrete(
    breaks = levels_x,
    labels = c(
      "Before HW",
      "End of the HW",
      "16 days post HW",
      "32 days post HW",
      "64 days post HW"
    )
  ) +
  stat_summary(data = all, 
               mapping = aes(x = timepoint_num, y = hill_simps, group = treatment, fill = treatment),
               fun = mean,
               geom = "point",
               shape = 23,
               size = 3,
               position = position_dodge(width = 0.8)
  ) +
  scale_color_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  scale_fill_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  # scale_y_log10() +
  labs(
    x = "Sampling point",
    y = expression(paste("Inverse Simpson"))) +
  theme_bw() +
  theme(legend.position = "none",
        #       legend.text = element_text(size = 14),
        #       legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("C")

hill_rich1 | hill_sh2 | hill_simps3 # Figure 2 in paper

    ## 4.3. Mean and standard errors of hill numbers----

hill_summary = all[, .(
  n = .N,
  
  mean_rich = mean(hill_rich, na.rm = TRUE),
  se_rich = sd(hill_rich, na.rm = TRUE) / sqrt(sum(!is.na(hill_rich))),
  
  mean_sh = mean(hill_sh, na.rm = TRUE),
  se_sh = sd(hill_sh, na.rm = TRUE) / sqrt(sum(!is.na(hill_sh))),
  
  mean_simps = mean(hill_simps, na.rm = TRUE),
  se_simps = sd(hill_simps, na.rm = TRUE) / sqrt(sum(!is.na(hill_simps)))
  
), by = .(timepoint, timepoint_num, treatment)]

setorder(hill_summary, timepoint_num)

  # 5. Beta diversity : Bray-Curtis dissimilarities----

    ## 5.1. Bacterial community overview: relative abundances of the ten most frequent phylum and genera----

      ### 5.1.1. Dataset preparation----

dt_rare = as.data.table(mat_rare, keep.rownames = "samples")
dt_rare = melt(dt_rare, id.vars = "samples", variable.name = "genus", value.name = "total_reads")
dt_rare = dt_rare[order(-total_reads)]
dt[, samples := sub("_F_filt\\.fastq\\.gz$", "", samples)]
tax_phylum = unique(dt[, .(samples = samples, genus = Genus, Phylum)]) # add a phylum column
dt_rare = merge(dt_rare, tax_phylum, by = c("samples", "genus"), all.x = TRUE) # add a phylum column to dt_rare
dt_rare[, condition := sub(".*_(ctrl|test)_.*", "\\1", samples)]
dt_rare[, timepoint := sub(".*_(T[1-5])_.*", "\\1", samples)]
dt_rare[, timepoint_num := as.integer(sub(".*_T([1-5])_.*", "\\1", samples))]
dt_rare[, unique(Phylum)]

dt_rare[, total_sample := sum(total_reads, na.rm = TRUE), by = samples]
dt_rare[, rel_abundance := total_reads / total_sample]
phylum_sample = dt_rare[ , .(rel_abundance = sum(rel_abundance, na.rm = TRUE)), by = .(samples, Phylum)]

summary = phylum_sample[
  ,
  .(
    mean = round(mean(rel_abundance) * 100, 3),
    sd = round(sd(rel_abundance) * 100, 3),
    n = .N,
    se = round(sd(rel_abundance) / sqrt(.N) * 100, 3),
    lower = round(mean(rel_abundance) - qt(0.975, .N - 1) * sd(rel_abundance) / sqrt(.N), 5) * 100,
    upper = round(mean(rel_abundance) + qt(0.975, .N - 1) * sd(rel_abundance) / sqrt(.N), 5) * 100
  ),
  by = Phylum
] # mean and se of relative abundances of the phylums

summary[order(mean, decreasing = T)]

      ### 5.1.2. The five most frequent phyla----

top10_phylum = dt_rare[, .(total_reads = sum(total_reads, na.rm = TRUE)), by = Phylum][order(-total_reads)][1:10]
top10_phylum[, relative_abundance := 100 * total_reads / sum(dt_rare$total_reads, na.rm = TRUE),]

ggplot(top10_phylum,
       aes(x = reorder(Phylum, total_reads),
           y = relative_abundance)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Phylum",
    y = "Relative abundance (%)",
    title = "Ten more abundant phyla"
  ) +
  theme_minimal()

filter = top10_phylum[order(relative_abundance, decreasing = T), Phylum][1:5]

dt_phyl = dt_rare[Phylum %in% filter,]

dt_phyl[, freq := total_reads/sum(total_reads), by = .(condition, timepoint_num)]

dt_phyl[Phylum == "Proteobacteria", Phylum := "Pseudomonadota"]
dt_phyl[Phylum == "Actinobacteria", Phylum := "Actinobacteriota"]
dt_phyl[Phylum == "Tenericutes", Phylum := "Mycoplasmatota"]
dt_phyl[Phylum == "Bacteroidetes", Phylum := "Bacteroidota"]
dt_phyl[Phylum == "Firmicutes", Phylum := "Bacillota"]

levels_x = unique(dt_phyl[, timepoint])

pal = c(
  "Actinobacteriota"   = "#FFD700",
  "Bacteroidota"       = "#458B00",
  "Bacillota"          = "#CD6090",
  "Pseudomonadota"     = "darkorange",
  "Mycoplasmatota"     = "#CD6839")


phylum = ggplot(dt_phyl, aes(y = freq, x = as.factor(timepoint_num), fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid2(~condition, labeller = labeller(condition = c("ctrl" = "Control", "test" = "Exposed")), 
              strip = strip_themed(
                background_x = list(
                  element_rect(fill = "#7f7f7f"),
                  element_rect(fill = "#cd4f39")
                )
              )) +
  scale_fill_manual(values = pal) +
  # scale_fill_rickandmorty() +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, 0.25),
                     labels = function(x) gsub("\\.00$", "", sprintf("%.2f", x))) +
  ylab("Relative abundance") +
  xlab("Sampling timepoint") +
  scale_x_discrete(breaks = levels_x, labels = c("Before HW", "End of the HW", "16 days post HW", "32 days post HW", "64 days post HW"), )  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "left", strip.text = element_text(
    face = "bold",
    colour = "white"),
    text = element_text(size = 16),
    legend.text = ggtext::element_markdown(size = 14), 
    plot.title = element_text(size=18, face="bold")) + ggtitle("A") # Figure 3A

      ### 5.1.3. The ten most frequent genera----

top10_gen = dt_rare[, .(total_reads = sum(total_reads, na.rm = TRUE)), by = genus][order(-total_reads)][1:10]
top10_gen[, relative_abundance := 100 * total_reads / sum(dt_rare$total_reads, na.rm = TRUE)]

ggplot(top10_gen,
       aes(x = reorder(genus, total_reads),
           y = relative_abundance)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Phylum",
    y = "Relative abundance (%)",
    title = "Ten more abundant phyla"
  ) +
  theme_minimal()

filter = top10_gen[order(relative_abundance, decreasing = T), genus][1:10]

dt_gen = dt_rare[genus %in% filter,]

dt_gen[, freq := total_reads/sum(total_reads), by = .(condition, timepoint_num)]
dt_gen[genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", genus := "Rhizobium spp."]
dt_gen[genus == "Burkholderia-Caballeronia-Paraburkholderia", genus := "B-C-P"]

pal = c(
  "Vibrio"                  = "#F2D06B",
  "Phyllobacterium"         = "#E28B2D",
  "Candidatus_Hepatoplasma" = "#F8F66E",
  "Streptomyces"            = "#8B5422",
  "Rhizobium spp."          = "#8ED1C6",
  "Pseudomonas"             = "#4F7332",
  "Promicromonospora"       = "palegreen3",
  "Paracoccus"              = "#8DC3E6",
  "Halomonas"               = "#2F3F73",
  "B-C-P"                   = "lightsalmon3"
)

genus = ggplot(dt_gen, aes(y = freq, x = as.factor(timepoint_num), fill = genus)) +
  geom_bar(stat = "identity") +
  facet_grid2(~condition, labeller = labeller(condition = c("ctrl" = "Control", "test" = "Exposed")), 
              strip = strip_themed(
                background_x = list(
                  element_rect(fill = "#7f7f7f"),
                  element_rect(fill = "#cd4f39")
                )
              )) +
  scale_fill_manual("Genus",values = pal,
                    labels = c("B-C-P",
                               "Candidatus <i>Hepatoplasma</i>",
                               "<i>Halomonas</i>",
                               "<i>Paracoccus</i>",
                               "<i>Phyllobacterium</i>",
                               "<i>Promicromonospora</i>",
                               "<i>Pseudomonas</i>",
                               "<i>Rhizobium</i> spp.",
                               "<i>Streptomyces</i>",
                               "<i>Vibrio</i>")) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, 0.25),
                     labels = function(x) gsub("\\.00$", "", sprintf("%.2f", x))) +
  ylab("Relative abundance") +
  xlab("Sampling timepoint") +
  scale_x_discrete(breaks = levels_x, labels = c("Before HW", "End of the HW", "16 days post HW", "32 days post HW", "64 days post HW"), )  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", strip.text = element_text(
    face = "bold",
    colour = "white"),
    text = element_text(size = 16),
    legend.text = ggtext::element_markdown(size = 14), 
    plot.title = element_text(size=18, face="bold")) + ggtitle("B") # Figure 3B

phylum | genus # Figure 3 in paper

    ## 5.2. Bray-Curtis dissimilarities matrix----

dt_rare[, prop := total_reads/sum(total_reads), by = samples]
dt_rare[, sum(prop), by = samples]
dt_rare[condition == "ctrl", condition := "Control"]
dt_rare[condition == "test", condition := "Exposed"]

d = dcast(dt_rare, samples ~ genus, value.var = "prop", fill = 0)
rownames(d) = d$samples
samples = d$samples
d = d[, -1]

dist_mat = vegdist(d, method = "bray")

nmds = metaMDS(d, distance = "bray", k = 2, trymax = 100)

meta = unique(dt_rare[, .(samples, timepoint, condition)])

rownames(d) = samples
meta = meta[match(rownames(d), samples)]
meta[, timepoint := as.factor(timepoint)]
meta[, condition := as.factor(condition)]

nmds$stress
stressplot(nmds)
plot(nmds)

scores = as.data.table(scores(nmds))
scores = cbind(scores, meta[, .(condition, timepoint)])

ggplot(scores, aes(x = sites.NMDS1, y = sites.NMDS2,colour = condition, fill = condition)) +
  stat_ellipse(
    aes(fill = condition, colour = condition),
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_point(size = 3) +
  scale_colour_manual(
    name = "Treatment",
    values = c("#7f7f7f", "#cd4f39")
  ) +
  scale_fill_manual(
    values = c("#7f7f7f", "#cd4f39"),
    guide = "none"
  ) +
  ylab("NMDS2") +
  xlab("NMDS1") +
  theme_minimal() +
  theme(text = element_text(size = 14))

ggplot(scores,
       aes(sites.NMDS1, sites.NMDS2,
           colour = condition)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  facet_wrap(~ timepoint) +
  scale_colour_manual(values = c("#7f7f7f", "#cd4f39"),
                      name = "Treatment") +
  theme_minimal() +
  theme(text = element_text(size = 14)) # NMDS by timepoint

    ## 5.3. Analyse 2: PERMANOVA----

# H0 : Groups do not differ in spread or position multivariate space.

disp1 = betadisper(dist_mat, meta$condition)
anova(disp1) # p > 0.05 -> similar variability between groups "condition" & "test"
permutest(disp1)
plot(disp1)

disp2 = betadisper(vegdist(dist_mat, method = "bray"), meta$timepoint)
anova(disp2)
permutest(disp2)
plot(disp2)
TukeyHSD(disp2)

set.seed(123)
terms = adonis2(dist_mat ~ condition * timepoint, data = meta, method = "bray", permutations = 9999, by = "terms") # global permanova, results not allow us to conclude if "condition" or "timepoint" have an effect on the bacterial community composition. But the interaction between this 2 factors is significant. So we have to test whether each factor has an effect independant from each other -> by = "margin"
terms

set.seed(123)
margin = adonis2(dist_mat ~ condition * timepoint, data = meta, method = "bray", permutations = 9999, by = "margin") # only the interaction between the two factors appeared in results, so we can only conclude that the composition of the bacterial community change along the time in each group of the factor "condition"
margin












