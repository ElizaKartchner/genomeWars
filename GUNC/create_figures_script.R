library(tidyverse)

gunc_file = read_tsv("GUNC.progenomes_2.1.maxCSS_level.tsv")
gunc_file = read_tsv("fake_data.tsv")

print(gunc_file)

tidy_gunc_file = gunc_file %>%
  select(genome, contamination_portion) %>%
  filter(!is.na(contamination_portion)) %>%
  mutate(contamination_portion = contamination_portion*100) %>%
  print()


my_plot = ggplot(tidy_gunc_file, aes(x=genome, y=contamination_portion)) + 
  geom_col(fill="skyblue4") + 
  ylab("Average Percent (%) of Contamination") +
  xlab("Genus") + 
  ylim(0, 100) + 
  theme_bw(10) + 
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Percent of Genome Contamination per Genus According to GUNC")

my_plot

ggsave("my_plot.png", plot = my_plot, width = 6, height = 4, units = "in")


