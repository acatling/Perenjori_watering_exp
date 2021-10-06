## Working through R for Data Science
# https://r4ds.had.co.nz/introduction.html
# 24/08/20

install.packages(c("nycflights13", "gapminder", "Lahman"))
library(tidyverse)
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy))
