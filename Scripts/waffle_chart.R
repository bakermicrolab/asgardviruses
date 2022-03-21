#Create waffle plots for VPF-Class VPF family membership ratios
#Ian Rambo 2021

library(tidyverse)
library(waffle)
library(rcartocolor)

setwd("../Data")

vpf_family <- read.table("vpf_family_virname.tsv",
                         header = TRUE, sep = "\t")

family <- read.table("family.tsv", header = TRUE)

virus_name <- c("Meg22_1012_scaffold_2kb_scaffold_226", "Meg22_1214_scaffold_2kb_scaffold_313",
               "Meg22_1012_scaffold_2kb_scaffold_91", "Meg22_1214_scaffold_2kb_scaffold_152",
               "Meg22_1012_scaffold_2kb_scaffold_548", "Meg22_1214_scaffold_2kb_scaffold_2849",
               "Guay1_scaffold_2kb_Guay1_scaffold00443")
proposed_name <- c("Fenrir Meg22_1012", "Fenrir Meg22_1214", "Nidhogg Meg22_1012",
                   "Nidhogg Meg22_1214", "Ratatoskr", "Skoll", "Vedfolnir")

vnames <- data.frame(virus_name, proposed_name)

vfam <- family %>% dplyr::left_join(vnames)



vpf_family <- vpf_family %>% mutate(membership_pct = membership_ratio * 100,
                                    membership_round = round(membership_ratio,
                                                             digits = 4)) %>%
  filter(virus_proposed_name != "Cluster")


palette <- c("#003f5c","#ffdb58","#665191","#a05195",
             "#d45087","#f95d6a","#ff7c43","#ffa600",
             "#004c6d","#0e7194","#2099ba","#2099ba",
             "#38c3de","#56efff")


waffle_filename <- "vpf_waffle_grid.pdf"

waffle_plot <- vfam %>%
  ggplot(aes(fill = class_name, values = membership_ratio)) +
  theme_enhance_waffle() +
  geom_waffle(make_proportional = TRUE,
              nrows = 10,
              colour = "white",
              flip = TRUE,
              size = 0.33,
              alpha = 0.7) +
  #scale_fill_manual(values = palette, name = "VPF Classification") +
  rcartocolor::scale_fill_carto_d(name = "VPF Classification", palette = "Safe") +
  theme(legend.text=element_text(size=10),
        legend.title=element_text(size=11),
        strip.text.x = element_text(size=12)) +
  facet_wrap(~proposed_name, nrow = 2)

ggsave(waffle_filename, plot = waffle_plot,
       dpi = 400, device = "pdf")

#=============================================================================
png(filename="treemap.png")
    #width=850, height=850, res=300)
treemap(vpf_family %>% filter(virus_proposed_name == "Nidhogg-1214"),
        index = "class_name",
        vSize = "membership_round",
        type = "index",
        # Main
        title="",
        palette="Dark2",
        # Borders:
        border.col=c("black"),
        border.lwds=1,
        # Labels
        inflate.labels=TRUE,
        fontcolor.labels="white",
        align.labels=c("center", "center")
        )
dev.off()
