# Sequencing depth in samples
library(data.table)
library(ggplot2)

depth = fread("/.../2_depth_AVUL.txt") # file generated from fastq files

depth[, sample_id := sub("_R[12]\\.fq\\.gz$", "", sample)]
depth = depth[, .(reads_R1_R2 = sum(reads)), by = sample_id]

p = ggplot(depth, aes(x = reorder(sample_id, reads_R1_R2), y = reads_R1_R2)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Sequencing depth by sample",
    x = "Sample",
    y = "Read counts"
  ) +
  theme_bw()

print(p)

sum(depth$reads_R1_R2) # number of sequenced reads
