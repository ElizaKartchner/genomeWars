library(tidyverse)
# FCS-GX
summary = read_tsv("summary_17000_expected.tsv")
head(summary)

div_one_levels <- c(
  "prok",
  "anml",
  "arch",
  "fung",
  "plnt",
  "prst",
  "synt"
)

summary %>% 
  separate(div, sep=":", into=c("div_1", "div_2")) %>%
  group_by(div_1) %>%
  summarise(count = n()) %>%
  mutate(div_1 = factor(div_1, levels = div_one_levels)) %>%
  mutate(div_1 = fct_recode(div_1, "Prokaryotes (Bacteria)" = "prok",
                            "Animals (Metazoa)" = "anml",
                            "Archaea" = "arch", 
                            "Fungi" = "fung", 
                            "Plants (Viridiplantae)" = "plnt", 
                            "Protists (Other Eukaryota)" = "prst", 
                            "Synthetic" = "synt"
                            ))  -> div_1_file

print(div_1_file)


summary %>% 
  separate(div, sep=":", into=c("div_1", "div_2")) %>%
  group_by(div_2) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  pull(div_2) -> ordered_div_2 
print(ordered_div_2)

summary %>% 
  separate(div, sep=":", into=c("div_1", "div_2")) %>%
  group_by(div_2) %>%
  summarise(count = n()) %>%
  mutate(div_2 = factor(div_2, levels = ordered_div_2 ))-> div_2_file
  
print(div_2_file, n=27)

level_one = ggplot(div_1_file, aes(x=div_1, y=count)) + 
  geom_col(fill="#4169E1") + 
  theme_bw(14) + 
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) +
  ylab("Number of Genomes") + 
  xlab("")

level_one 

ggsave("figure_4a.png", plot = level_one, width = 8, height = 5, units = "in")


level_two = ggplot(div_2_file, aes(x=div_2, y=count)) + 
  geom_col(fill="#4169E1") + 
  theme_bw(14) + 
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  ylab("Number of Genomes") + 
  xlab("")

level_two

ggsave("figure_4b.png", plot = level_two, width = 8, height = 5, units = "in")