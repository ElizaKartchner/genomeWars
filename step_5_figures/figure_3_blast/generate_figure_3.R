library(tidyverse) 
library(scales)

unique_sequences = read_csv("unique_sequences.csv")

print(unique_sequences)

unique_sequences %>%
  mutate(eVal = as.character(eVal)) -> unique_sequences

values <- c(0.1, 1.00E-10, 1.00E-20,
            1.00E-30,  1.00E-40,  1.00E-50, 
            1.00E-60, 1.00E-70, 1.00E-80)

big_graph = ggplot(unique_sequences, aes(x=eVal, y=num_in_bucket)) + 
  geom_col(fill="dodgerblue2") + 
  theme_bw(14) + 
  xlab("E-Value") + 
  ylab("Number of Genomes") + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1))

ggsave("big_graph.png", plot = big_graph, width = 6, height = 4, units = "in")


unique_sequences_tunc = read_csv("unique_sequences_tunc.csv")

small_graph = ggplot(unique_sequences_tunc, aes(x=eVal, y=num_in_bucket)) + 
  geom_col(fill="dodgerblue2") + 
  theme_bw(14) + 
  xlab("E-Value") + 
  ylab("Number of Genomes") + 
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1))

ggsave("small_graph.png", plot = small_graph, width = 6, height = 4, units = "in")


figure_5 = ggplot(unique_sequences, aes(x=eVal, y=Num_Unique_Sequences)) + 
  geom_point(color="dodgerblue2", size=2) + 
  theme_bw() + 
  ylab("Number of Unique Sequences") + 
  xlab("E-value") + 
  #scale_x_log10(labels = scales::scientific_format(), breaks = values) + 
  scale_x_continuous(trans  = compose_trans("log10", "reverse"), breaks = values)

figure_5

ggsave("figure_5.png", plot = figure_5, width = 6, height = 4, units = "in")