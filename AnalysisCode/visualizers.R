# Developing boxplots for accuracy communication.
library(ggplot2)
library(tidyverse)

# Conumee Accuracies
synthStats <- data.frame(Chromosomes = A42W$Chromosomes, 
                         A42W = A42W$popAccAvg,
                         A48O = A48O$popAccAvg,
                         A48P = A48P$popAccAvg,
                         A6Z2 = A6Z2$popAccAvg,
                         A7EL = A7EL$popAccAvg,
                         A7EN = A7EN$popAccAvg,
                         A9HX = A9HX$popAccAvg) %>% 
  pivot_longer(
    cols = c(A42W, A48O, A48P, A6Z2, A7EL, A7EN, A9HX), 
    names_to = "Cases",
    values_to = "Accuracies",
  )
synthStats$Chromosomes <- as.factor(as.numeric(synthStats$Chromosomes))
  
chromConumee <- ggplot(synthStats, aes(x=Chromosomes, y=Accuracies)) + 
  geom_boxplot() + geom_jitter(position = position_jitter(0.01)) + labs(title = "Conumee Accuracies per Chromosome")
chromConumee

casesConumee <- ggplot(synthStats, aes(x=Cases, y=Accuracies)) + 
  geom_boxplot() + geom_jitter(position = position_jitter(0.01)) + labs(title = "Conumee Accuracies per Case")
casesConumee
