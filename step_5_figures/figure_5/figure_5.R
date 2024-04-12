library(tidyverse) 

unique_sequences = read_csv("unique_sequences.csv")

print(unique_sequences)

values <- c(0.1, 1.00E-10, 1.00E-20,
            1.00E-30,  1.00E-40,  1.00E-50, 
            1.00E-60, 1.00E-70, 1.00E-80)


figure_5 = ggplot(unique_sequences, aes(x=eVal, y=Num_Unique_Sequences)) + 
  geom_point(color="dodgerblue2", size=2) + 
  theme_bw() + 
  ylab("Number of Unique Sequences") + 
  xlab("E-value") + 
  scale_x_log10(labels = scales::scientific_format(), breaks = values) 

figure_5

ggsave("figure_5.png", plot = figure_5, width = 6, height = 4, units = "in")