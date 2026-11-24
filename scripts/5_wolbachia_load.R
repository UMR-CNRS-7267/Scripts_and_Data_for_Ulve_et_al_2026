# Wolbachia load analysis in ovaries, nerve cords and gut contents

library(data.table)
library(stringi)
library(lme4)
library(DHARMa)
library(emmeans)
library(mgcv)
library(car)
library(vegan)
library(patchwork)
library(ggplot2)
library(ggpmisc)

# 1. Data import and formatting----

setwd("/Users/.../")

wolbiall = read.table("5_1__Wolbachia_loads_NC_OV.txt", dec = ".", stringsAsFactors = T, h = T)
summary(wolbiall)
setDT(wolbiall)

wolbiall[, treatment := ifelse(grepl("ctrl", Name), "Control", "Exposed")] # rename treatment, Control = unexposed to the heatWave, Exposed = exposed to the heatwave

wolbiall[, organ := ifelse(grepl("CN", Name), "NC", "OV")] # create new column to rename organ: NC = nervous chain, OV = ovaries

wolbiall[, timepoint := stri_extract_first(Name, regex = "T[0-9]")] # extraction of timepoints
wolbiall[, timepoint_num := as.integer(stri_extract_first_regex(Name, "(?<=T)\\d+"))] # transform sampling timepoints as a quantitative value

wolbiall[, Experimental_day := fcase(
  timepoint_num == 1, 0,
  timepoint_num == 2, 8,
  timepoint_num == 3, 16,
  timepoint_num == 4, 32,
  timepoint_num == 5, 64
)] # and extraction of sampling timepoints as a real date

wolbiall[timepoint_num == 1, period := "pre-HW"] # splitting sampling timepoints into two categories (data from pre- and post- exposure to the heatwave - HW) to analyze Wolbachia load before and after HW. 
wolbiall[timepoint_num %in% c(2:5), period := "post-HW"]
wolbiall$period = factor(wolbiall$period, levels = c("pre-HW", "post-HW"))

wolbiall = wolbiall[!(Name == "T1_ctrl_CN_11") | (Plate_1 == 1)] # Remove individual used as internal control used in all qPCR plate
wolbiall[, (names(wolbiall)) := lapply(.SD, function(x) if (is.character(x)) as.factor(x) else x)] # Conversion of "character" variables to "factor" variables

  # 2. Calculate Wolbachia load and mean Cp (crossing point) by sample----

wolbiall[, load := 2^-(Cp_Wsp - Cp_Tleu)]#Cp
wolbiall = wolbiall[wolbiall$load<40,] # If the Cp is greater than 40, then we consider the individual to be uninfected.

MeanDose = wolbiall[, .(Wolbachia_Load = mean(load)), by = .(Name, treatment, timepoint,timepoint_num,Experimental_day, organ, Plate_1, period)]

  # 3. Statistical Analysis: Has the heatwave exposure affected Wolbachia load?----

    ## 3.1. In nerve cords (NC)----

NC = MeanDose[MeanDose$organ=="NC",]

gam_mod_NC1 = gam(Wolbachia_Load ~ treatment + treatment : timepoint_num + s(timepoint_num, k = 5),
              data = NC, method = "REML") # GAM with interaction treatment x experiment time

summary(gam_mod_NC1)
anova(gam_mod_NC1) # check whether the interaction is significant
gam.check(gam_mod_NC1)
qq.gam(gam_mod_NC1)

gam_mod_NC2 = gam(Wolbachia_Load ~ treatment + s(timepoint_num, k = 5), data = NC, method = "REML") # Test the effect of timepoint_num (non-linear smooth effect) on Wolbachia load, while controlling for treatment

summary(gam_mod_NC2)
anova(gam_mod_NC2)  # significant non-linear smooth effect -> different temporal patterns by treatment? (gam_mod_NC3)
gam.check(gam_mod_NC2)
qq.gam(gam_mod_NC2)

gam_mod_NC3 = gam(Wolbachia_Load ~ treatment + s(timepoint_num, k = 5, by = treatment), data = NC, method = "REML") # different temporal patterns by treatment?
summary(gam_mod_NC3)
anova(gam_mod_NC3)  # different temporal patterns by treatment
gam.check(gam_mod_NC3)
qq.gam(gam_mod_NC3)

anova(gam_mod_NC2, gam_mod_NC3, test = "F") # allowing treatment-specific smooth terms significantly improved the model

      ### 3.1.1. Plotting GAM predictions (Figure 4A)----

newdata = expand.grid(
  timepoint_num = seq(min(NC$timepoint_num), max(NC$timepoint_num), length.out = 200),
  treatment = unique(NC$treatment)
)
newdata$pred = predict(gam_mod_NC3, newdata, se.fit = TRUE)$fit
newdata$se = predict(gam_mod_NC3, newdata, se.fit = TRUE)$se.fit

levels_x = sort(unique(NC$timepoint_num))

Figure_4A = ggplot() +
  geom_boxplot(
    data = NC,
    aes(x = as.factor(timepoint_num), y = Wolbachia_Load, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = NC,
    aes(x = timepoint_num, y = Wolbachia_Load, color = treatment),
    position = position_dodge(width = 0.8)
  ) +
  geom_line(
    data = newdata,
    aes(x = timepoint_num, y = pred, color = treatment, group = treatment),
    size = 1.2
  ) +
  geom_ribbon(
    data = newdata,
    aes(
      x = timepoint_num,
      ymin = pred - 2 * se,
      ymax = pred + 2 * se,
      fill = treatment,
      group = treatment
    ),
    alpha = 0.15, color = NA
  ) +
  stat_summary(data = NC, 
    mapping = aes(x = timepoint_num, y = Wolbachia_Load, group = treatment, fill = treatment),
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
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
  scale_color_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  scale_fill_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a")) +
  labs(
    x = "Sampling point",
    y = expression(paste(italic("Wolbachia"), " load in nerve cords"))) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("A")

Figure_4A

    ## 3.2. In ovaries----

OV = MeanDose[MeanDose$organ=="OV",]# OV = ovaries

gam_mod_OV1 = gam(Wolbachia_Load ~ treatment + treatment : timepoint_num + s(timepoint_num, k = 5),
                  data = OV, method = "REML") # GAM with interaction treatment x experiment time

summary(gam_mod_OV1)
anova(gam_mod_OV1)  # check whether the interaction is significant
gam.check(gam_mod_OV1)
qq.gam(gam_mod_OV1)

gam_mod_OV2 = gam(Wolbachia_Load ~ treatment + s(timepoint_num, k = 5), data = OV, method = "REML") # Test the effect of timepoint_num (non-linear smooth effect) on Wolbachia load, while controlling for treatment

summary(gam_mod_OV2)
anova(gam_mod_OV2) # no significant effect of s(timepoint_num) -> used linear model (LM) instead of GAM
gam.check(gam_mod_OV2)
qq.gam(gam_mod_OV2)

lm1 = lm(Wolbachia_Load ~ treatment + treatment : timepoint_num, data = OV) 
summary(lm1)
Anova(lm1, type = "II") # no significant interaction -> LM without interaction (lm2)
plot(simulateResiduals(lm1))

lm2 = lm(Wolbachia_Load ~ treatment + timepoint_num, data = OV)
summary(lm2)
Anova(lm2)
plot(simulateResiduals(lm2))

    ### 3.2.1. Plotting data without GAM predictions (Figure 4B)----

Figure_4B = ggplot() +
  geom_boxplot(
    data = OV,
    aes(x = as.factor(timepoint_num), y = Wolbachia_Load, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = OV,
    aes(x = timepoint_num, y = Wolbachia_Load, color = treatment),
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(data = OV, 
               mapping = aes(x = timepoint_num, y = Wolbachia_Load, group = treatment, fill = treatment),
               fun = mean,
               geom = "point",
               shape = 23,
               size = 3,
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
    y = expression(paste(italic("Wolbachia"), " load in ovaries"))) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("B")

Figure_4B

lm3 = lm(Wolbachia_Load ~ treatment, data = OV[timepoint_num<2,])
summary(lm3)
Anova(lm3)
plot(simulateResiduals(lm3))

    ## 3.3. In gut contents----

      ### 3.3.1. Data import and formatting----

dt = fread("/Users/.../5_2_Wolbachia_reads_gut.txt") # those data come from the "alpha_beta_diversity.R"
setnames(dt, "total_reads", "counts")

dt = dt[, total_reads := sum(counts), by = samples]

freq_genus = dt[!Genus == "", .(total_counts = sum(counts), total_reads = sum(unique(total_reads)), n_samples = uniqueN(samples)), by = Genus] # to know most frequent Genus and their frequencies. Notice that unassigned reads to the Genus level were discarded
freq_genus[, freq_percent := (total_counts / total_reads) * 100]
setorder(freq_genus, -freq_percent)

dt = dt[Genus == "Wolbachia",]
dt = dt[order(total_reads, decreasing = T)]
dt[, samples := gsub(pattern = "*_F_filt.fastq.gz", replacement = "", x = samples)]

hist(dt$total_reads)

dt[grepl("AVUL", samples), timepoint := sub(".*(T[1-5]).*", "\\1", samples)] # add a column "timepoint" (T1, T2, T3, T4, T5)
dt[, timepoint_num := as.integer(sub(".*([1-5]).*", "\\1", timepoint))]
dt[grepl("AVUL", samples), condition := sub(".*(ctrl|test).*", "\\1", samples)] # add a column "condition" (ctrl/test)
dt[, condition := as.factor(condition)]
setnames(x = dt, "condition", "treatment")
dt[treatment == "ctrl", treatment := "Control"]
dt[treatment == "test", treatment := "Exposed"]

      ### 3.3.2. Means and standard errors of Wolbachia-assigned reads----

summary = dt[
  ,
  .(
    mean = round(mean(counts) * 100, 3),
    sd = round(sd(counts) * 100, 3),
    n = .N,
    se = round(sd(counts) / sqrt(.N) * 100, 3),
    lower = round(mean(counts) - qt(0.975, .N - 1) * sd(counts) / sqrt(.N), 5) * 100,
    upper = round(mean(counts) + qt(0.975, .N - 1) * sd(counts) / sqrt(.N), 5) * 100
  ),
  by = .(timepoint, treatment)
]

summary[order(mean, decreasing = T)]

levels_x = sort(unique(dt$timepoint_num))

      ### 3.3.3. Visualization of Wolbachia-assigned reads numbers----

ggplot() +
  geom_boxplot(
    data = dt,
    aes(x = as.factor(timepoint_num), y = counts, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = dt,
    aes(x = timepoint_num, y = counts, color = treatment),
    position = position_dodge(width = 0.8)
  ) +
  stat_summary(data = dt, 
               mapping = aes(x = timepoint_num, y = counts, group = treatment, fill = treatment),
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
  ylab("Reads assigned to Wolbachia in gut content") +
  xlab("Sampling point") +
  scale_color_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  scale_fill_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("C")

      ### 3.3.4. Has the heatwave exposure affected the Wolbachia density in gut contents?----

gam_mod_TD1 = gam(counts ~ treatment + treatment : timepoint_num + s(timepoint_num, k = 5), offset(log(total_reads)), family = nb(), method = "REML", data = dt) # GAM with interaction treatment x experiment time

summary(gam_mod_TD1)
anova(gam_mod_TD1) # check whether the interaction is significant
gam.check(gam_mod_TD1)
qq.gam(gam_mod_TD1)

gam_mod_TD2 = gam(counts ~ treatment + s(timepoint_num, k = 5), offset(log(total_reads)), family = nb(), method = "REML", data = dt) # Test the effect of timepoint_num (non-linear smooth effect) on Wolbachia load, while controlling for treatment

summary(gam_mod_TD2)
anova(gam_mod_TD2)# significant non-linear smooth effect -> different temporal patterns by treatment? (gam_mod_NC3)
gam.check(gam_mod_TD2)
qq.gam(gam_mod_TD2)

gam_mod_TD3 = gam(counts ~ treatment + s(timepoint_num, by = treatment, k = 5), offset(log(total_reads)), family = nb(), method = "REML", data = dt) # different temporal patterns between treatments?

summary(gam_mod_TD3)
anova(gam_mod_TD3) # different temporal patterns by treatment
gam.check(gam_mod_TD3)
qq.gam(gam_mod_TD3)

anova(gam_mod_TD2, gam_mod_TD3, test = "F") # allowing treatment-specific smooth terms significantly improved the model

      ### 3.3.5. Plotting GAM predictions (Figure 4C)----

newdata = expand.grid(
  timepoint_num = seq(min(dt$timepoint_num),
                      max(dt$timepoint_num),
                      length.out = 100),
  treatment = c("Control", "Exposed")
)

pred_link = predict(gam_mod_TD3, newdata = newdata, type = "link", se.fit = TRUE)

newdata$fit = gam_mod_TD3$family$linkinv(pred_link$fit)
newdata$lower = gam_mod_TD3$family$linkinv(pred_link$fit - 1.96 * pred_link$se.fit)
newdata$upper = gam_mod_TD3$family$linkinv(pred_link$fit + 1.96 * pred_link$se.fit)

levels_x = sort(unique(dt$timepoint_num))

Figure_4C = ggplot() +
  geom_boxplot(
    data = dt,
    aes(x = as.factor(timepoint_num), y = counts, fill = treatment),
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = dt,
    aes(x = timepoint_num, y = counts, color = treatment),
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
  stat_summary(data = dt, 
               mapping = aes(x = timepoint_num, y = counts, group = treatment, fill = treatment),
               fun = mean,
               geom = "point",
               shape = 23,
               size = 3,
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
  scale_color_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  scale_fill_manual("Treatment", values = c("Control" = "#7f7f7f","Exposed" = "#cd394a"), labels = c("Control", "Exposed")) +
  labs(
    x = "Sampling point",
    y = expression(paste("Reads counts assigned to ", italic("Wolbachia"), " in gut contents"))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 14, face = "bold"), plot.title = element_text(size=18, face="bold")) + ggtitle("C")

Figure_4A | Figure_4B | Figure_4C # Figure 4 in paper



