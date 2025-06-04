library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggthemes)
library(ggprism)
library(envalysis)
library(readxl)
library(tidyr)

setwd("~/Desktop/PhD/2025 First Term/QPCR/20250604/")

#read the data
data <- read_excel("stimulation_ra_final_list.xlsx")

#have the logarithmic fold changes
data$Log2FC_SOD2 <- log2(data$FC_SOD2)
data$Log2FC_HIF1a <- log2(data$FC_HIF1a)
data$Log2FC_GLUT1 <- log2(data$FC_GLUT1)

#filter for only a stimulation group and patient group to show it in a graph
selected_stimulation <- "IFNg"
selected_group <- "RA"
selected_data <- data %>%
  filter(Stimulation == selected_stimulation & Patient_Group == selected_group) %>%
  select(Patient, Patient_Group, Log2FC_HIF1a, Log2FC_SOD2, Log2FC_GLUT1)
#or select(Patient, Patient_Group, Log2FC_CD74, Log2FC_HLA_DRA)

#reshape the dataframe
selected_data <- selected_data %>%
  pivot_longer(
    cols = starts_with("Log2FC"),
    names_to = "Gene",
    values_to = "Log2FC"
  ) %>%
  drop_na(Log2FC) %>% #one of the GLUT1 measurement is removed (it is NA)
  mutate(Gene = recode(Gene,
                       "Log2FC_HIF1a" = "HIF1a",
                       "Log2FC_SOD2" = "SOD2",
                       "Log2FC_GLUT1" = "GLUT1"))


# calculate p values with t-test
pvals <- selected_data %>%
  group_by(Gene) %>%
  summarise(p.value = t.test(Log2FC, mu = 0)$p.value) %>%
  mutate(label = paste0("p = ", signif(p.value, 3)))

#draw the plot
p <- ggplot(selected_data, aes(x = Gene, y = Log2FC)) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    aes(fill = Gene),
    width = 0.5,
    linewidth = 0.5,
    color = "black")+
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.3,
    linewidth = 0.5
  ) + 
  geom_violin(alpha = 0.3, position = "dodge",aes(fill=Gene)) +
  geom_jitter(width = 0.15, shape = 21, fill = "white", stroke = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  labs(title = paste0("Relative Expressions\n", "post-", selected_stimulation, " Stimulation\n","In ",selected_group), 
        x = "Genes", y = "Log Fold Change") +
  theme_publish() +
  theme(legend.position = "none",plot.background = element_rect(fill = "transparent",colour = NA),
         axis.text.x = element_text(angle = 0)) +
  scale_fill_prism(palette = "candy_soft") +
  geom_text(data = pvals, aes(x = Gene, y = max(selected_data$Log2FC) + 0.6, label = label)) 

#see the plot
print(p)

#saving the plot
ggsave(
  plot = p,
  filename = paste0(selected_stimulation, "_stimulated_RA.png"),
  bg = "transparent",
  dpi = 1000,
  width = 5, height = 4
)

