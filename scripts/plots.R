library(tidyverse)
library(cowplot)

#load data
#rawdata <- readr::read_csv("processed/rawdata.csv") 
pruneddata <- readr::read_csv("processed/pruneddata.csv") %>%
  dplyr::mutate(TOF = 1.67110915*TOF + 100.31575759,
                norm.EXT = 15.68506057*norm.EXT -1.04595184,
                volume = (pi/4)*TOF*(norm.EXT)^2,
                timepoint = case_when(hour == 1~"01",
                                      hour == 2~"02",
                                      hour == 3~"03",
                                      hour == 4~"04",
                                      hour == 5~"05",
                                      hour == 6~"06",
                                      hour == 7~"07",
                                      hour == 8~"08",
                                      hour == 9~"09",
                                      hour %in% c(10:42)~as.character(hour)),
                stage = case_when(strain %in% c("N2", "dpy-1", "dpy-5", "lon-3") & hour %in% c(1:10) ~ "L1",
                                          strain == "N2" & hour %in% c(16:22) ~ "L2",
                                          strain == "N2" & hour %in% c(27:31) ~ "L3",
                                          strain == "lon-3" & hour %in% c(16:22) ~ "L2",
                                          strain == "lon-3" & hour %in% c(27:32) ~ "L3",
                                          strain == "dpy-1" & hour %in% c(19:26) ~ "L2",
                                          strain == "dpy-1" & hour %in% c(31:35) ~ "L3",
                                          strain == "dpy-5" & hour %in% c(19:24) ~ "L2",
                                          strain == "dpy-5" & hour %in% c(30:34) ~ "L3")) %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2", "dpy-1","dpy-5", "lon-3")))



## length vs width plot
pruneddata %>%
  dplyr::filter(TOF < 1000) %>%
  ggplot(.) +
  aes(x = TOF, y = norm.EXT, color = stage) +
  geom_point(alpha = 0.4, size = 0.009) +
  theme_cowplot(14) +
  facet_wrap(~strain, scales = "free") + panel_border() +
  scale_color_manual(values = c('#4477AA', '#117733', '#DDCC77',"grey")) +
  labs(x="Median Length", y = "Median Width", title = "Population Data") + guides(color = "none")

# medians
pruneddata %>%
  dplyr::group_by(strain, hour, timepoint, row, col) %>%
  dplyr::summarize(median.Length = median(TOF),
                   median.Width = median(norm.EXT), .groups = "drop") %>%
  dplyr::mutate(stage = case_when(strain %in% c("N2", "dpy-1", "dpy-5", "lon-3") & hour %in% c(1:10) ~ "L1",
                                  strain == "N2" & hour %in% c(16:22) ~ "L2",
                                  strain == "N2" & hour %in% c(27:31) ~ "L3",
                                  strain == "lon-3" & hour %in% c(16:22) ~ "L2",
                                  strain == "lon-3" & hour %in% c(27:32) ~ "L3",
                                  strain == "dpy-1" & hour %in% c(19:26) ~ "L2",
                                  strain == "dpy-1" & hour %in% c(31:35) ~ "L3",
                                  strain == "dpy-5" & hour %in% c(19:24) ~ "L2",
                                  strain == "dpy-5" & hour %in% c(30:34) ~ "L3")) %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2", "dpy-1","dpy-5", "lon-3"))) %>%
  ggplot(.) +
  aes(x = median.Length, y = median.Width) +
  geom_point(size = 0.3) +
  theme_cowplot(18) +
  scale_color_manual(values = c('#4477AA', '#117733', '#DDCC77',"grey")) +
  facet_wrap(~strain, scales = "free") + panel_border() +
  labs(x="Median Length", y = "Median Width") + guides(color= "none")


## plots split out to compare to N2
plot.red <- function(strains, colorfill) {
  pruneddata %>%
    dplyr::filter(strain %in% strains) %>%
    dplyr::group_by(strain, hour, timepoint, row, col) %>%
    dplyr::summarize(median.red = median(red/volume), .groups = "drop") %>%
    ggplot(.) +
    aes(x = hour, y = median.red, fill = strain, color = strain) + 
    geom_boxplot(aes(x = timepoint),outlier.alpha = 0, alpha = 0.7) +
    geom_jitter(size=0.1, width=0.2, alpha=0.2) +
    geom_smooth(span = 0.18, se = F, size = 0.8, method = "loess", aes(color = strain)) +
    theme_cowplot(font_size = 22, rel_small = 16/22) +
    labs(x="Time (hours)", y = "Median Fluorescence") +
    scale_fill_manual(values = colorfill) +
    scale_color_manual(values = colorfill) +
    theme(legend.position="top", legend.background = element_rect(fill="#FFFFFF")) +
    scale_x_discrete(breaks = c("05", 10, 15, 20, 25, 30, 35, 40, 45, 50))
  
}

plot.length <- function(strains, colorfill) {
  pruneddata %>%
    dplyr::filter(strain %in% strains) %>%
    dplyr::group_by(strain, hour, timepoint, row, col) %>%
    dplyr::summarize(median.Length = median(TOF), .groups = "drop") %>%
    ggplot(.) +
    aes(x = hour, y = median.Length, fill = strain, color = strain) + 
    geom_boxplot(aes(x = timepoint),outlier.alpha = 0, alpha = 0.7) +
    geom_jitter(size=0.1, width=0.2, alpha=0.2) +
    theme_cowplot(font_size = 22, rel_small = 16/22) +
    labs(x="Time (hours)", y = "Median Length") +
    scale_fill_manual(values = colorfill) +
    scale_color_manual(values = colorfill) +
    theme(legend.position="top", legend.background = element_rect(fill="#FFFFFF")) +
    scale_x_discrete(breaks = c("05", 10, 15, 20, 25, 30, 35, 40, 45, 50))
}

plot.width <- function(strains, colorfill) {
  pruneddata %>%
    dplyr::filter(strain %in% strains) %>%
    dplyr::group_by(strain, hour, timepoint, row, col) %>%
    dplyr::summarize(median.Width = median(norm.EXT), .groups = "drop") %>%
    ggplot(.) +
    aes(x = hour, y = median.Width, fill = strain, color = strain) + 
    geom_boxplot(aes(x = timepoint),outlier.alpha = 0, alpha = 0.7) +
    geom_jitter(size=0.1, width=0.2, alpha=0.2) +
    theme_cowplot(font_size = 22, rel_small = 16/22) +
    labs(x="Time (hours)", y = "Median Width") +
    scale_fill_manual(values = colorfill) +
    scale_color_manual(values = colorfill) +
    theme(legend.position="top", legend.background = element_rect(fill="#FFFFFF")) +
    scale_x_discrete(breaks = c("05", 10, 15, 20, 25, 30, 35, 40, 45, 50))
}

plot.volume <- function(strains, colorfill) {
  pruneddata %>%
    dplyr::filter(strain %in% strains) %>%
    dplyr::group_by(strain, hour, timepoint, row, col) %>%
    dplyr::summarize(median.Width = median(volume), .groups = "drop") %>%
    ggplot(.) +
    aes(x = hour, y = median.Width, fill = strain, color = strain) + 
    geom_boxplot(aes(x = timepoint),outlier.alpha = 0, alpha = 0.7) +
    geom_jitter(size=0.1, width=0.2, alpha=0.2) +
    theme_cowplot(font_size = 22, rel_small = 16/22) +
    labs(x="Time (hours)", y = "Median Volume") +
    scale_fill_manual(values = colorfill) +
    scale_color_manual(values = colorfill) +
    theme(legend.position="top", legend.background = element_rect(fill="#FFFFFF"))  +
    scale_x_discrete(breaks = c("05", 10, 15, 20, 25, 30, 35, 40, 45, 50)) +
    scale_y_log10()
}

#1200x700
a <- plot.volume(c("N2", "lon-3", "dpy-1", "dpy-5"), c("#e69f00","#009e73","#0072b2","black"))
b <- plot.volume(c("N2", "lon-3"), c("#0072b2","black"))
c <- plot.volume(c("N2", "dpy-1"), c("#e69f00","black"))
d <- plot.volume(c("N2", "dpy-5"), c("#009e73","black"))

cowplot::plot_grid(a, b, c, d, ncol = 2, nrow = 2, align = "hv")




a <- plot.red(c("N2", "lon-3", "dpy-1", "dpy-5"), c("#e69f00","#009e73","#0072b2","black"))
b <- plot.length(c("N2", "lon-3", "dpy-1", "dpy-5"), c("#e69f00","#009e73","#0072b2","black"))
c <- plot.width(c("N2", "lon-3", "dpy-1", "dpy-5"), c("#e69f00","#009e73","#0072b2","black"))

cowplot::plot_grid(a, b, c, ncol = 2, nrow = 2, align = "hv")



clean %>% 
  dplyr::mutate(TOF = 1.67110915*TOF + 100.31575759,
                norm.EXT = 15.68506057*norm.EXT -1.04595184) %>%
  #dplyr::filter(strain == "dpy-5") %>%
  dplyr::group_by(hour, strain, plate, row, col) %>%
  #dplyr::group_by(hour, strain) %>%
  dplyr::summarize(median.TOF = median(TOF), median.norm.EXT = median(norm.EXT)) %>%
  ungroup() %>%
  ggplot() + 
  #aes(x = TOF, y = norm.EXT) +
  aes(x = median.TOF, y = median.norm.EXT) + 
  geom_point(aes(color = hour), size = 0.5) +
  #geom_point(size = 0.01, color="grey", alpha=0.5) + 
  theme_bw() + labs(x="Median Length", y="Median Width") + 
  facet_wrap(~strain, scales = "free") 
  geom_smooth(se = F, span = 0.05)
  geom_path(aes(color = hour))



