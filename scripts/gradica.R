library(tidyverse)

grid1 <- read_csv("~/Finite-Populations-/r_n1.csv")
grid2 <- read_csv("~/Finite-Populations-/r_n2.csv")
media1 <- read_csv("~/Finite-Populations-/mean1.csv")
media2 <- read_csv("~/Finite-Populations-/mean2.csv")
tiempo1 <- read_csv("~/Finite-Populations-/runtime1.csv")[,2]
tiempo2 <- read_csv("~/Finite-Populations-/runtime2.csv")[,2]


full_grid <- bind_rows(grid1,grid2)
full_mean <- bind_rows(media1,media2)
full_time <- bind_rows(tiempo1,tiempo2)

max1 <- max(media1)
max2 <- max(media2)
min1 <- min(media1)
min2 <- min(media2)

modelo <- data_frame(r_n = full_grid$r_n,
                     mean = full_mean$mean,
                     Runtime = full_time$min) %>% 
                     mutate( label = as.factor(case_when( r_n <= 5000 ~ "Grid1",
                                                TRUE ~ "Grid2")))
                     
          


modelo %>% ggplot(aes(x = r_n, y = mean)) +
           geom_point(aes(  col = Runtime, shape = label), size = 3) +
           geom_hline(yintercept = max1, col = 'orange', size = 2, alpha = 0.7) +
           geom_hline(yintercept = min1, col = 'orange', size = 2, alpha = 0.7) +
           geom_hline(yintercept = max2, col = 'red', size = 2, alpha = 0.7) +
           geom_hline(yintercept = min2, col = 'red', size = 2, alpha = 0.7) +
           geom_line(col = "#79BC42") +
           scale_y_continuous(breaks = seq(from=min1,
                                           to = max1,
                                           length.out = 10))+
           theme_pro() +
           labs(title = "Posterior Total Mean",
                subtitle = "N = 100,000",
                caption = "Units in Thousands",
                x = "R - n",
                y = "Total Population Mean")
  
