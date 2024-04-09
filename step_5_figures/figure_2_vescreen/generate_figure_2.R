library(tidyverse)

data = read_tsv("expected_vecscreen_17000_output.tsv")

genomes_w_alignments = length(unique(data$query_file))
genomes_w_no_align = 17868 - genomes_w_alignments

final_data = tibble(labels = c("Genomes with alignments", "Genomes with no alignments"),
                    values = c(genomes_w_alignments, genomes_w_no_align))

ggplot(final_data, aes(x = labels, y = values)) +
  geom_col(fill = "dodgerblue2") +
  labs(y = "Counts", x = "") +
  geom_text(aes(label = values), vjust = -0.25) +
  theme_bw()
