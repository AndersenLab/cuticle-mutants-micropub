library(tidyverse)
library(cowplot)
library(gganimate)
library(gifski)

#load data
#rawdata <- readr::read_csv("processed/rawdata.csv") 
pruneddata <- readr::read_csv("data/pruneddata.csv") %>%
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
                                      hour %in% c(10:42)~as.character(hour))) %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2", "dpy-1","dpy-5", "lon-3"))) %>%
  dplyr::mutate(Genotype = factor(strain, levels = c("N2", "dpy-1","dpy-5", "lon-3"), labels = c("WT", "dpy-1","dpy-5", "lon-3")))



#mean tp -- no facet
anim <- pruneddata %>%
  dplyr::group_by(strain, Genotype, hour, timepoint) %>%
  dplyr::summarize(median.TOF = median(TOF),
                   median.norm.EXT= median(norm.EXT), 
                   median.norm.red = median(norm.red), .groups = "drop") %>%
  dplyr::mutate(stage = case_when(strain %in% c("N2", "dpy-1", "dpy-5", "lon-3") & hour %in% c(1:10) ~ "L1",
                                  strain == "N2" & hour %in% c(16:22) ~ "L2",
                                  strain == "N2" & hour %in% c(26:32) ~ "L3",
                                  strain == "lon-3" & hour %in% c(16:22) ~ "L2",
                                  strain == "lon-3" & hour %in% c(26:30) ~ "L3",
                                  strain == "dpy-1" & hour %in% c(19:26) ~ "L2",
                                  strain == "dpy-1" & hour %in% c(31:35) ~ "L3",
                                  strain == "dpy-5" & hour %in% c(19:24) ~ "L2",
                                  strain == "dpy-5" & hour %in% c(29:34) ~ "L3"))  %>%
  dplyr::mutate(hour = as.numeric(timepoint),
                `Molting?` = dplyr::case_when(median.norm.red <= 0.06 & timepoint > 12 ~ "yes", TRUE ~ "no")) %>%
  dplyr::filter(hour <= 55) %>%
  ggplot(.) +
  aes(x = median.TOF, y = median.norm.EXT, color = Genotype) +
  geom_point(size = 2, alpha = 1) + 
  scale_color_manual(values =  c("#e69f00","#009e73","#0072b2","black"), labels = c("WT", expression(italic("dpy-1")),
                                                                                     expression(italic("dpy-5")),
                                                                                     expression(italic("lon-3")))) +
  theme(legend.position="top", legend.background = element_rect(fill="#FFFFFF")) +
  #guides(color = "none") +
  theme_cowplot(font_size = 20, rel_small = 16/20) + 
  #facet_wrap(~strain, scales = "free") + panel_border() +
  transition_time(hour) + shadow_mark(future = F, alpha = alpha/2, size = size/2) +
  #transition_manual(hour, cumulative = TRUE) +
  labs(x = expression(paste("Median Animal Length")), y = expression(paste("Median Animal Width")))

animate(anim, height = 1920, width = 2400, units = "px", res = 300)

anim_save("figures/sumGIF.gif")

#facet anim
anim <- pruneddata %>%
  dplyr::group_by(strain, Genotype, hour, timepoint, row, col) %>%
  dplyr::summarize(median.TOF = median(TOF),
                   median.norm.EXT= median(norm.EXT), 
                   median.norm.red = median(norm.red), .groups = "drop") %>%
  dplyr::mutate(stage = case_when(strain %in% c("N2", "dpy-1", "dpy-5", "lon-3") & hour %in% c(1:10) ~ "L1",
                                  strain == "N2" & hour %in% c(16:22) ~ "L2",
                                  strain == "N2" & hour %in% c(26:32) ~ "L3",
                                  strain == "lon-3" & hour %in% c(16:22) ~ "L2",
                                  strain == "lon-3" & hour %in% c(26:30) ~ "L3",
                                  strain == "dpy-1" & hour %in% c(19:26) ~ "L2",
                                  strain == "dpy-1" & hour %in% c(31:35) ~ "L3",
                                  strain == "dpy-5" & hour %in% c(19:24) ~ "L2",
                                  strain == "dpy-5" & hour %in% c(29:34) ~ "L3"))  %>%
  dplyr::mutate(hour = as.numeric(timepoint),
                `Molting?` = dplyr::case_when(median.norm.red <= 0.06 & timepoint > 12 ~ "yes", TRUE ~ "no")) %>%
  dplyr::filter(hour <= 55) %>%
  dplyr::mutate(Genotype = factor(strain, levels = c("N2", "dpy-1","dpy-5", "lon-3"), labels = c("WT", expression(italic("dpy-1")),
                                                                                                 expression(italic("dpy-5")),
                                                                                                 expression(italic("lon-3"))))) %>%
  ggplot(.) +
  aes(x = median.TOF, y = median.norm.EXT, color = Genotype) +
  geom_point(size = 1, alpha = 1) +
  scale_color_manual(values = c("#e69f00","#009e73","#0072b2","black")) +
  theme(legend.position="top", legend.background = element_rect(fill="#FFFFFF")) +
  guides(color = "none") +
  theme_cowplot(font_size = 20, rel_small = 16/20) + 
  facet_wrap(~Genotype, scales = "free", labeller = label_parsed) + panel_border() +
  transition_time(hour) + shadow_mark(future = F, alpha = alpha/4, size = size/4) +
  labs(x = expression(paste("Median Animal Length")), y = expression(paste("Median Animal Width")))

animate(anim, height = 1920, width = 2400, units = "px", res = 300)

anim_save("figures/facetGIF.gif")
