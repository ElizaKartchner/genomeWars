library(tidyverse)

gunc_file = read_tsv("expected_output_17000.tsv")

print(gunc_file)

tidy_gunc_file = gunc_file %>%
  select(genome, contamination_portion) %>%
  filter(!is.na(contamination_portion)) %>%
  mutate(contamination_portion = contamination_portion*100) %>%
  filter(contamination_portion > 25) %>%
  print()


figure_1 = ggplot(tidy_gunc_file, aes(x = genome, y = contamination_portion)) + 
  geom_jitter(position = position_jitter(width = 0.2), size = 2, color = "blue") +
  ylab("Percent (%) of Contamination") +
  xlab("Genome") + 
  ylim(0, 100) + 
  theme_bw(10) + 
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) 

figure_1 

ggsave("figure_1.png", plot = my_boxplot_small, width = 6, height = 4, units = "in")