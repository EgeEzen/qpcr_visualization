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

setwd("~/Desktop/PhD/2025 First Term/QPCR/20250619/")

#read the data
data <- read_excel("thp1_stimulation_final_list.xlsx")

#have the logarithmic fold changes
data$Log2FC_GLUT1 <- log2(data$FC_GLUT1)
data$Log2FC_EPAS1 <- log2(data$FC_EPAS1)
data$Log2FC_SPP1 <- log2(data$FC_SPP1)
data$Log2FC_IRF7 <- log2(data$FC_IRF7)
data$Log2FC_SOCS3 <- log2(data$FC_SOCS3)

selected_data <- data %>% filter(Stimulation %in% c("LPS","SPP1")) %>%
  select(Hours_After, Stimulation,Log2FC_GLUT1, Log2FC_EPAS1,Log2FC_SPP1,Log2FC_IRF7,Log2FC_SOCS3) %>%
  mutate(Hours_After = factor(Hours_After, c("3","6","12","24")))

#reshape the dataframe
selected_data <- selected_data %>%
  pivot_longer(
    cols = starts_with("Log2FC"),
    names_to = "Gene",
    values_to = "Log2FC"
  ) %>%
  drop_na(Log2FC) %>% #one of the GLUT1 measurement is removed (it is NA)
  mutate(Gene = recode(Gene,
                       "Log2FC_EPAS1" = "EPAS1",
                       "Log2FC_SPP1" = "SPP1",
                       "Log2FC_GLUT1" = "GLUT1",
                       "Log2FC_IRF7" = "IRF7",
                       "Log2FC_SOCS3" = "SOCS3"))

# calculate p values with t-test (only one biological replicate, we cannot do this?)
# pvals <- selected_data %>%
#   group_by(Gene,Hours_After) %>%
#   summarise(p.value = t.test(Log2FC, mu = 0)$p.value) %>% 
#   mutate(label = paste0("p = ", signif(p.value, 2)))

#draw the plot
p <- ggplot(selected_data, aes(x = Hours_After, y = Log2FC, color = Stimulation, group=Stimulation)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Gene,scales ="free_x")+
  geom_point(size=2, shape = 21, fill = "white", stroke = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_color_manual(values = c("LPS" = "#D95F02", "SPP1" = "#1B9E77")) +
  labs(title= "Log2FC Compared to CTRL over Time",
    x = "Hours After Stimulation",
    y = "Log2 Fold Change"
  ) +
  theme_publish() +
  theme(legend.position = "right", plot.title = element_text(size=14,face="bold"))
#see the plot
print(p)

#saving the plot
ggsave(
  plot = p,
  filename = "thp1_stimulation_over_time.png",
  bg = "transparent",
  dpi = 1000,
  width = 6, height = 5
)

