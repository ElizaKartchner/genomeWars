library(tidyverse)

gunc_file = read_tsv("GUNC.progenomes_2.1.maxCSS_level.tsv")

print(gunc_file)

tidy_gunc_file = gunc_file %>%
  select(genome, contamination_portion) %>%
  filter(!is.na(contamination_portion)) %>%
  print()


my_plot = ggplot(tidy_gunc_file, aes(x=genome, y=contamination_portion)) + 
  geom_col(fill="skyblue4") + 
  ylab("Portion of Genome Contaminated") +
  xlab("Genome by RefSeq Id") + 
  ylim(0, 1) + 
  theme_bw(10) + 
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Portion of genome contamination according to GUNC")

ggsave("my_plot.png", plot = my_plot, width = 6, height = 4, units = "in")

