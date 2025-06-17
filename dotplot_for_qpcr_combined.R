library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggthemes)
library(ggprism)
library(envalysis)
library(readxl)
library(tidyr)
library(ggrepel)

setwd("~/Desktop/PhD/2025 First Term/QPCR/20250520")

#read the data
data <- read_excel("stimulation_oa_healthy_final_list.xlsx")

#have the logarithmic fold changes
data$Log2FC_SOD2 <- log2(data$FC_SOD2)
data$Log2FC_HIF1a <- log2(data$FC_HIF1a)
data$Log2FC_GLUT1 <- log2(data$FC_GLUT1)
data$Log2FC_CD74 <- log2(data$FC_CD74)
data$Log2FC_IL6 <- log2(data$FC_IL6)
data$Log2FC_HLA_DRA <- log2(data$FC_HLA_DRA)
data$Log2FC_CHI3L1 <- log2(data$FC_CHI3L1)
data$Log2FC_EPAS1 <- log2(data$FC_EPAS1)

#filter for only a stimulation group and patient group to show it in a graph
selected_stimulation <- "IFNg"
selected_data <- data %>%
  filter(Stimulation == selected_stimulation ) %>%
  select(Patient, Patient_Group, Log2FC_HIF1a, Log2FC_SOD2, 
         Log2FC_GLUT1,Log2FC_IL6,Log2FC_CHI3L1, Log2FC_EPAS1)
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
                       "Log2FC_GLUT1" = "GLUT1",
                       "Log2FC_CD74" = "CD74",
                       "Log2FC_IL6" = "IL6", 
                       "Log2FC_HLA_DRA" = "HLA-DRA",
                       "Log2FC_EPAS1" = "EPAS1",
                       "Log2FC_CHI3L1" = "CHI3L1"))

# calculate p values with t-test
pvals <- selected_data %>%
  group_by(Gene,Patient_Group) %>%
  summarise(p.value = t.test(Log2FC, mu = 0)$p.value) %>% 
  mutate(label = paste0("p = ", signif(p.value, 2)))

#draw the plot
p <- ggplot(selected_data, aes(x = Patient_Group, y = Log2FC)) +
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
  facet_wrap(~Gene)+
  geom_violin(alpha = 0.6, position = "dodge",aes(fill=Patient_Group)) +
  geom_jitter(width = 0.15, shape = 21, fill = "white", stroke = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  labs(title = paste0("Relative Expressions\n", "post-", selected_stimulation), 
       x = NULL, y = "Log Fold Change") +
  theme_publish() +
  theme(legend.position = "none",plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(angle = 0)) +
  scale_fill_manual(values = c("Healthy" = "#9FC2CC", "OA" ="#694D75")) +
  geom_text(data = pvals, aes(x = Patient_Group, y = max(selected_data$Log2FC) + 0.6, label = label)) 
# geom_label_repel(
 #    aes(label = Patient),
 #    label.padding = 0.5,
 #    box.padding = 0.5,
 #    point.padding = 0.5,
 #    size = 3, alpha =0.8) 
  

#see the plot
print(p)

#saving the plot
ggsave(
  plot = p,
  filename = paste0(selected_stimulation,"_stimulated_w_violin_healthy_and_OA_1.png"),
  bg = "transparent",
  dpi = 1000,
  width = 6, height = 6
)

