library(tidyverse)

# Create a dataset a large dataset
set.seed(123)
large_fake_data <- data.frame(
  genome = rep(genera, each = 100),
  contamination_portion = c(rnorm(100, mean = 7, sd = 3),    
                            rnorm(100, mean = 3, sd = 5),    
                            rnorm(100, mean = 18, sd = 8),  
                            rnorm(100, mean = 26, sd = 10),  
                            rnorm(100, mean = 48, sd = 15),  
                            rnorm(100, mean = 2, sd = 2),    
                            rnorm(100, mean = 15, sd = 5),   
                            rnorm(100, mean = 65, sd = 20),  
                            rnorm(100, mean = 21, sd = 8)
                            )
)

print(large_fake_data)

# Create a boxplot on large fake dataset
my_boxplot_large = ggplot(large_fake_data, aes(x=genome, y=contamination_portion)) + 
  geom_boxplot(width = 0.6, fill = "lightblue", color = "blue") + 
  #geom_jitter(position = position_jitter(width = 0.2), size = 2, color = "red") +
  ylab("Percent (%) of Contamination") +
  xlab("Genus") + 
  ylim(0, 100) + 
  theme_bw(12) + 
  theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) 

my_boxplot_large

ggsave("my_boxplot_large.png", plot = my_boxplot_large, width = 6, height = 4, units = "in")

