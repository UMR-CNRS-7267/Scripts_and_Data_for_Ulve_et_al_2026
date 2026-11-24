# Temperatures to construct the "summer temperatures" program and the "heatwave temperatures" program

library(data.table)
library(readxl)
library(patchwork)
library(ggplot2)

  # 1. Import data and formatting----

all = fread("/.../1_temperatures.txt")

all[, time := paste0(DATE, "_", hour)] # creation of the "time" column including the date and the hour
all$time = as.POSIXct(all$time, format = "%Y_%m_%d_%H") # conversion of the time column as POSIXct format
setorder(all, time) # sort by time for shift

transition_points = copy(all) # transition points: last line of each group before a change of period. transition points will allow to plot a continuous line with geom_line().
transition_points[, period_next := shift(period, type = "lead")] # extract lines where period change of value
transition_points = transition_points[period != period_next]
transition_points[, period := period_next]
transition_points[, period_next := NULL]

all_fixed = rbindlist(list(all, transition_points), use.names = TRUE, fill = TRUE)
setorder(all_fixed, time)

all_fixed[period == "pre_hw", head(unique(DATE), 5)]
all_fixed[period == "pre_hw", tail(unique(DATE), 5)]

all_fixed[period == "post_hw", head(unique(DATE), 5)]
all_fixed[period == "post_hw", tail(unique(DATE), 5)]

test = all_fixed[DATE %in% c("2025_03_03", "2025_03_04", "2025_03_05", "2025_03_06", "2025_03_07"),]

test = all_fixed[!DATE %in% c("2025_06_01", "2025_06_02", "2025_06_03", "2025_06_04", "2025_06_05"),]

all_fixed[, day := as.Date(time)]
last_pre = all_fixed[period == "pre_hw", max(day, na.rm = TRUE)]
last_pre

pre_keep = all_fixed[period == "pre_hw" & day >= (last_pre - 4)]

hw_keep = all_fixed[period == "hw"]

post_keep = all_fixed[period == "post_hw" & day <= min(day[period == "post_hw"]) + 10]

all_extract = rbind(pre_keep, hw_keep, post_keep)
origin = as.Date("2025-03-03")

  # 2. Plotting temperatures of the "heatwave temperatures" condition----

heatwave = ggplot(all_extract, aes(time, temp, col = period)) +
  geom_line(linewidth = 2) +
  xlab("Experiment time") +
  ylab("Temperatures (°C)") +
  scale_colour_manual(
    values = c(
      "pre_hw" = "#7f7f7f",
      "hw" = "#cd4f39",   
      "post_hw" = "#7f7f7f"
    )) +
  scale_x_datetime(
    breaks = as.POSIXct(c(
      "2025-03-21 00:00:00",  # J0
      "2025-03-26 00:00:00",  # J23
      "2025-04-03 00:00:00"  # J31
    )),
    labels = c("J0", "J23", "J31"), expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.line = element_line(linewidth = 1.1),
    axis.ticks = element_line(linewidth = 1.4),
    axis.ticks.length = unit(3, "mm"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
    )

annotate(
  "text",
  x = max(all_extract$time),
  y = -Inf,
  label = "J94",
  vjust = 2)

  # 3. Plotting temperatures of the "summer temperatures" condition----

ref = all_extract[24:47, temp]
dates = unique(all_extract[, .(day)])

dt_summer = data.table(temp = rep(ref, times = 24))
dt_summer[, time := as.numeric(rownames(dt_summer))]
dt_summer[, period := "summer_condition"]

summer = ggplot(dt_summer, aes(time, y = temp, group = period)) +
  geom_line(linewidth = 2, colour = "#7f7f7f") +
  xlab("Experiment time") +
  ylab("Temperatures (°C)") +
  scale_y_continuous(
    limits = c(13, 38),
    breaks = c(15, 20, 25, 30, 35),
    expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.line = element_line(linewidth = 1.1),
    axis.ticks = element_line(linewidth = 1.4),
    axis.ticks.length = unit(3, "mm"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

heatwave/summer # Figure 1 in paper

