library(tidyverse)
library(dplyr)

grid1 <- read_csv("~/Finite-Populations-/graph/r_n1.csv")
grid2 <- read_csv("~/Finite-Populations-/graph/r_n2.csv")
grid3 <- read_csv("~/Finite-Populations-/graph/r_n3.csv")
grid4 <- read_csv("~/Finite-Populations-/graph/r_n4.csv")
grid5 <- read_csv("~/Finite-Populations-/graph/r_n5.csv")

media1 <- read_csv("~/Finite-Populations-/graph/mean1.csv")
media2 <- read_csv("~/Finite-Populations-/graph/mean2.csv")
media3 <- read_csv("~/Finite-Populations-/graph/mean3.csv")
media4 <- read_csv("~/Finite-Populations-/graph/mean4.csv")
media5 <- read_csv("~/Finite-Populations-/graph/mean5.csv")

tiempo1 <- read_csv("~/Finite-Populations-/graph/runtime1.csv")
tiempo2 <- read_csv("~/Finite-Populations-/graph/runtime2.csv")
tiempo3 <- read_csv("~/Finite-Populations-/graph/runtime3.csv")
tiempo4 <- read_csv("~/Finite-Populations-/graph/runtime4.csv")
tiempo5 <- read_csv("~/Finite-Populations-/graph/runtime5.csv")

full_grid <- bind_rows(
              bind_rows(
                bind_rows(
                  bind_rows(grid1,grid2),grid3),grid4),grid5)

full_mean <- bind_rows(
              bind_rows(
                bind_rows(
                  bind_rows(media1,media2),media3),media4),media5)

full_time <- bind_rows(
              bind_rows(
                bind_rows(
                  bind_rows(tiempo1,tiempo2),tiempo3),tiempo4),tiempo5)


max1 <- max(media1)
max2 <- max(media2)
max3 <- max(media3)
max4 <- max(media4)
max5 <- max(media5)

min1 <- min(media1)
min2 <- min(media2)
min3 <- min(media3)
min4 <- min(media4)
min5 <- min(media5)

modelo <- data_frame(r_n = full_grid$r_n,
                     mean = full_mean$mean,
                     Runtime = full_time$runtime) %>% 
                     mutate( Grid = as.factor(case_when( r_n <= 5000 ~ "Grid 1",
                                                          r_n > 5000 & r_n <= 10000 ~ "Grid 2",
                                                          r_n > 10000 & r_n <= 30000 ~ "Grid 3",
                                                          r_n > 30000 & r_n <= 60000 ~ "Grid 4",
                                                          TRUE ~ "Grid 5")),
                             Max = case_when(r_n <= 5000 ~ max1,
                                             r_n > 5000 & r_n <= 10000 ~ max2,
                                             r_n > 10000 & r_n <= 30000 ~ max3,
                                             r_n > 30000 & r_n <= 60000 ~ max4,
                                             TRUE ~ max5),
                             Min = case_when(r_n <= 5000 ~ min1,
                                             r_n > 5000 & r_n <= 10000 ~ min2,
                                             r_n > 10000 & r_n <= 30000 ~ min3,
                                             r_n > 30000 & r_n <= 60000 ~ min4,
                                             TRUE ~ min5))
                             
                             
                     
modelo1 <- modelo %>% 
           gather("variable", "value", -c(r_n,Runtime,Grid,Max,Min))

modelo1 %>% ggplot(aes(x = r_n, y = value, group = variable)) +
            geom_point(aes(shape = Grid, fill = Runtime), size = 3) +
            scale_shape_manual(values=c(21, 22, 23,24,25))+
            geom_vline(xintercept = 1600, show.legend = 'R-N TLC', col = "purple", size = 1)+
            geom_text(aes(x=7000, 
                label= paste("R-N CLT ~",1600), 
                y=1.5628), 
            text=element_text(size=5),
            family = 'Arial') +
            geom_line(col = "#79BC42") +
            geom_line(aes(y = Max,group = Grid, linetype = Grid))+
            geom_line(aes(y = Min, group = Grid,  linetype = Grid)) +
            scale_linetype_manual(values = c("solid", "dashed", "dotted", "twodash","dotdash"))+
            scale_y_continuous(breaks = signif(seq(from=min1,
                                                   to = max1,
                                                   length.out = 10),5))+
            scale_x_continuous(breaks = seq(from = 0,
                                  to = 90000,
                                  by = 2000)) +
            theme_pro() +
            theme(
              axis.text.x = element_text(angle = 60)
            ) +
            labs(title = "Posterior Total Mean",
                 subtitle = "N = 100,000",
                 caption = "Units in Thousands",
                 x = "R - n",
                 y = "Total Population Mean")

ecm <- read_csv('graph/ecm.csv', col_names = FALSE)

ecm <- ecm^2

graph <- data_frame(r_n = seq(from = 50, to = 2500, by = 50),
                    ecm = sqrt(colMeans(ecm)))

source('scripts/helpers.R')
graph %>% ggplot(aes(x = r_n, y = ecm)) + 
          geom_line( size = 0.65, col = '#79BC42') +
          labs(title = 'Error Cuadrático Medio',
               x = 'R - n',
               y = 'ECM')+
          scale_x_continuous(breaks = seq(from = 50, to = 2600, by = 250)) +
          theme_pro()
