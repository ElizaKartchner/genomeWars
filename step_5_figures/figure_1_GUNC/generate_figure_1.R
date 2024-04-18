library(tidyverse)

gunc_file = read_tsv("expected_output_17000.tsv")

print(gunc_file)

tidy_gunc_file = gunc_file %>%
  select(genome, contamination_portion) %>%
  filter(!is.na(contamination_portion)) %>%
  mutate(contamination_portion = contamination_portion*100) %>%
  filter(contamination_portion > 25) %>%
  mutate(genome = factor(genome)) %>%
  mutate(genome = fct_recode(genome, 
                             "NC_005832.1" =  "NC_005832_1", 
                             "NC_010154.1" =  "NC_010154_1", 
                             "NC_008361.1" =  "NC_008361_1", 
                             "NC_003214.2" =  "NC_003214_2", 
                             "NC_002628.3" =  "NC_002628_3", 
                             "NC_010155.1" =  "NC_010155_1", 
                             "NC_010152.1" =  "NC_010152_1",
                             "NC_001731.1" =  "NC_001731_1", 
                             "NC_005946.1" =  "NC_005946_1", 
                             "NC_001902.1" =  "NC_001902_1", 
                             "NC_005336.1" =  "NC_005336_1", 
                             "NC_006268_1" =  "NC_006268_1")) %>%
  print()

genomes = c(NC_005832_1, 
            NC_010154_1, 
            NC_008361_1, 
            NC_003214_2, 
            NC_002628_3, 
            NC_010155_1, 
            NC_010152_1,
            NC_001731_1, 
            NC_005946_1, 
            NC_001902_1, 
            NC_005336_1, 
            NC_006268_1)


figure_1 = ggplot(tidy_gunc_file, aes(x = genome, y = contamination_portion)) + 
  geom_col(fill="dodgerblue2") + 
  ylab("Percent (%) of Contamination") +
  xlab("Genome") + 
  ylim(0, 100) + 
  theme_bw(10) + 
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) 

figure_1 

ggsave("figure_1.png", plot = figure_1, width = 6, height = 4, units = "in")