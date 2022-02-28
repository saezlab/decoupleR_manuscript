library(ggplot2)
library(tidyverse)
library(patchwork)

# Create dir
path_figs <- file.path('figures')
dir.create(path_figs, showWarnings = F, recursive = T)

vals <- c(10,50,250,1250)

df <- as_tibble(read.csv('data/prc/scale_results.txt', sep = ' ')) %>%
  filter(!is.na(time)) %>%
  mutate(memr=memr/1e6, time=time/60)
df$mat <- df$mat %>% map(function(m){vals[[m]]}) %>% unlist()
df$avg <- df$time * 60 / df$mat / 250
m_time <- df %>% group_by(lang) %>% summarize(m=median(avg))
df$mat <- factor(df$mat, levels = c('10', '50', '250', '1250'))

# X seconds per sample and regulator


colnames(df) <- c('Language', '# samples', 'Method', 'Time(m)', 'Memory(GB)')
time_plot <- ggplot(data=df, aes(x=`# samples`, y=`Time(m)`, color=Language, group = Language)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Method, ncol = 11, strip.position = 'top') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#fed140", '#2269bb')) +
  theme(axis.title.x=element_blank(), strip.background=element_rect(fill="white"))

memo_plot <- ggplot(data=df, aes(x=`# samples`, y=`Memory(GB)`, color=Language, group = Language)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Method, ncol = 11, strip.position = 'top') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#fed140", '#2269bb')) +
  theme(strip.background=element_rect(fill="white"))

# Merge together and save
pdf(file = file.path(path_figs, 'supp_fig_7.pdf'),
    width = 9, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches
time_plot / memo_plot +
  plot_layout(guides = 'collect')
dev.off()
