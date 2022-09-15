library(easysorter)
library(tidyverse)
library(mclust)
library(cowplot)

## load sorter files
dirs <- "raw/20220222"
raw <- easysorter::read_data(dirs)
score <- raw[[1]]

plate <- score %>%
  dplyr::filter(call50 == "object" & !condition =="Wash") %>%
  dplyr::mutate(hour = as.numeric(condition))

#readr::write_csv(plate, path = "processed/rawdata.csv")

#################################################
# performing clustering using Mclust package ####
#################################################

cleanData <- function(data, s) {
  X <- data
  m <- X %>%
    dplyr::mutate(logTOF = log(TOF), logEXT = log(EXT)) %>%
    dplyr::select(logTOF, logEXT) %>%
    Mclust(., G = 3)
  medians <- X %>%
    bind_cols(.,as_tibble(m$classification)) %>% 
    dplyr::group_by(value) %>%
    dplyr::summarize(cluster = median(TOF), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(cluster)) %>%
    dplyr::pull()
  raw <- X %>%
    dplyr::bind_cols(.,as_tibble(m$classification)) %>% 
    dplyr::group_by(value) %>%
    dplyr::mutate(median = median(TOF), cluster_median = median) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(median)) %>%
    dplyr::mutate_at(vars(cluster_median),
                     ~dplyr::case_when(. == medians[1] ~ 1,
                                       . == medians[2] ~ 2,
                                       . == medians[3] ~ 3,
                                       TRUE ~ 4)) %>%
    dplyr::select(-value)
  
  return(raw)
}

clustr <- plate %>%
  dplyr::mutate(tp=hour, strn = strain) %>%
  dplyr::group_by(tp, strn) %>% 
  tidyr::nest()

t_final <- purrr::map_dfr(clustr$data, ~cleanData(.x)) %>%
  dplyr::group_by(hour, strain, cluster_median, median) %>%
  tidyr::nest()

t_clean <- t_final %>%
  tidyr::unnest(., cols = c(data)) %>%
  ungroup()

# plot to view
tof <- t_clean %>%  
  dplyr::mutate(cluster_median = as_factor(cluster_median)) %>%
  dplyr::filter(TOF < 750, !cluster_median == 3) %>%
  ggplot2::ggplot() + aes(x = hour, y = TOF, color = cluster_median) +
  geom_jitter(size=0.2, width=0.25, alpha=0.4) +
  theme_bw(12) + labs(x="Time (hours)") + theme(legend.position="bottom") +
  facet_wrap(~strain, ncol = 1) + panel_border() 

norm.ext <- t_clean %>%  
  dplyr::mutate(cluster_median = as_factor(cluster_median)) %>%
  dplyr::filter(TOF < 750, !cluster_median == 3) %>%
  ggplot2::ggplot() + aes(x = hour, y = norm.EXT, color = cluster_median) +
  geom_jitter(size=0.2, width=0.25, alpha=0.4) +
  theme_bw(12) + labs(x="Time (hours)") + theme(legend.position="bottom") +
  facet_wrap(~strain, ncol = 1) + panel_border() 

cowplot::plot_grid(tof, norm.ext, nrow = 1, ncol = 2, align = "hv")


######################################
# filtering out clusters/outliers ####
######################################
# function to remove outliers (will use for norm.red)
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

clean_raw <- t_clean %>%
  dplyr::filter(!cluster_median == 3,
                !(strain == "dpy-1" & cluster_median == 1 & hour %in% c(4,6,7,8,35:42)),
                !(strain == "dpy-1" & cluster_median == 2 & hour %in% c(2,3,12,13,14,16,17,18,20,22,24,25,27:34)),
                !(strain == "dpy-5" & cluster_median == 1 & hour %in% c(6,18)),
                !(strain == "dpy-5" & cluster_median == 2 & hour %in% c(2,5,10,12,13,15,16,17,19,21,23,25,28:42)),
                !(strain == "lon-3" & cluster_median == 1 & hour %in% c(1,2,4,9,34)),
                !(strain == "lon-3" & cluster_median == 2 & hour %in% c(5,8,11,12,15,16:21,23:28,30:33,35:42)),
                !(strain == "N2" & cluster_median == 1 & hour %in% c(4,7,31,33)),
                !(strain == "N2" & cluster_median == 2 & hour %in% c(8:12,14,16:29,34:42))) %>%
  group_by(hour, strain) %>%
  dplyr::mutate(red_outlier = (remove_outliers(norm.red)),
                red_outlier = ifelse(is.na(red_outlier), TRUE, FALSE),
                TOF_outlier = (remove_outliers(TOF)),
                TOF_outlier = ifelse(is.na(TOF), TRUE, FALSE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!red_outlier == T, !TOF_outlier)

# plot to view
tof <- clean_raw %>%  
  dplyr::mutate(cluster_median = as_factor(cluster_median)) %>%
  ggplot2::ggplot() + aes(x = hour, y = TOF, color = cluster_median) +
  geom_jitter(size=0.2, width=0.25, alpha=0.4) +
  theme_bw(12) + labs(x="Time (hours)") + theme(legend.position="bottom") +
  facet_wrap(~strain, ncol = 1) + panel_border() 

norm.ext <- clean_raw %>%  
  dplyr::mutate(cluster_median = as_factor(cluster_median)) %>%
  ggplot2::ggplot() + aes(x = hour, y = norm.EXT, color = cluster_median) +
  geom_jitter(size=0.2, width=0.25, alpha=0.4) +
  theme_bw(12) + labs(x="Time (hours)") + theme(legend.position="bottom") +
  facet_wrap(~strain, ncol = 1) + panel_border() 

cowplot::plot_grid(tof, norm.ext, nrow = 1, ncol = 2, align = "hv")

#readr::write_csv(clean_raw, path = "processed/pruneddata.csv")



















