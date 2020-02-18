library(tidyverse)

grid1 <- read_csv("~/Finite-Populations-/graph/r_n1.csv")
grid2 <- read_csv("~/Finite-Populations-/graph/r_n2.csv")
grid3 <- read_csv("~/Finite-Populations-/graph/r_n3.csv")
media1 <- read_csv("~/Finite-Populations-/graph/mean1.csv")
media2 <- read_csv("~/Finite-Populations-/graph/mean2.csv")
media3 <- read_csv("~/Finite-Populations-/graph/mean3.csv")
tiempo1 <- read_csv("~/Finite-Populations-/graph/runtime1.csv")
tiempo2 <- read_csv("~/Finite-Populations-/graph/runtime2.csv")
tiempo3 <- read_csv("~/Finite-Populations-/graph/runtime3.csv")

full_grid <- bind_rows(bind_rows(grid1,grid2),grid3)
full_mean <- bind_rows(bind_rows(media1,media2),media3)
full_time <- bind_rows(bind_rows(tiempo1,tiempo2),tiempo3)

max1 <- max(media1)
max2 <- max(media2)
max3 <- max(media3)
min1 <- min(media1)
min2 <- min(media2)
min3 <- min(media3)

modelo <- data_frame(r_n = full_grid$r_n,
                     mean = full_mean$mean,
                     Runtime = full_time$runtime) %>% 
                     mutate( label = as.factor(case_when( r_n <= 5000 ~ "Grid 1",
                                                          r_n > 5000 & r_n <= 10000 ~ "Grid 2",
                                                          TRUE ~ "Grid 3")))
                     
          


modelo %>% ggplot(aes(x = r_n, y = mean)) +
           geom_point(aes(  col = Runtime, shape = label), size = 3) +
           geom_hline(yintercept = max1, col = 'orange', size = 2, alpha = 0.7) +
           geom_hline(yintercept = min1, col = 'orange', size = 2, alpha = 0.7) +
           geom_hline(yintercept = max2, col = 'red', size = 2, alpha = 0.7) +
           geom_hline(yintercept = min2, col = 'red', size = 2, alpha = 0.7) +
           geom_hline(yintercept = max3, col = 'yellow', size = 2, alpha = 0.7) +
           geom_hline(yintercept = min3, col = 'yellow', size = 2, alpha = 0.7) +
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
  
