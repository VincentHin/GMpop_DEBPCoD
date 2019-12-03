## Prepare data
FirstRepro_sub <- subset(popbif.DIST_RESAMP.R0, par1 %in% c(0,10,30), c('par1', 'dist_resamp', 'IDagefirstrepro_zero', 'IDagefirstweaning_zero'))
FirstReproWean.melted <- melt(FirstRepro_sub, id.vars = c('dist_resamp', 'par1'))
FirstReproWean.melted$par1 <- factor(FirstReproWean.melted$par1, levels = c('0','10','30'))
FirstRepro_means <- aggregate(list(value = FirstReproWean.melted$value), by = list('par1' = FirstReproWean.melted$par1, 'dist_resamp' = FirstReproWean.melted$dist_resamp, 'variable' = FirstReproWean.melted$variable), mean, na.rm = TRUE)
FirstRepro_means$labels <- round(FirstRepro_means$value,1)

dodge <- position_dodge2(width=0.75)

## Create Figure 5
Figure5 <-
  ggplot(data = FirstReproWean.melted,
    mapping = aes(x = par1, y = value, fill = variable)) +
  geom_boxplot(varwidth = FALSE, size = 0.2, outlier.size = 0.2) +
  scale_fill_manual(name = element_blank(),
    values = c(GMcolours$class3, GMcolours$class1),
    labels = c('When first calf is born', 'When first calf is weaned')) +
  guides(fill = guide_legend(nrow = 1, title = element_blank())) +
  facet_grid(. ~ dist_resamp, labeller = labeller(dist_resamp = labels_dist_resamp)) +
  # coord_flip() +
  ylab('Age (years)') + xlab('Disturbance period') +
  ggopts +
  theme(legend.box.margin = margin(-7,0,-3,0))

## Add mean values as points (shape = 4)
Figure5_plusMeans <-
  Figure5 + geom_point(stat = "summary", fun.y = mean, size = 0.2, position = dodge, shape = 4)

## Add mean values as labels next to boxplots
Figure5_plusText <- 
  Figure5_plusMeans + 
  geom_text(data = subset(FirstRepro_means, variable == "IDagefirstrepro_zero"), 
    mapping = aes(x = as.factor(par1), label = labels), nudge_x = 0.15, size = 1.7) +
  geom_text(data = subset(FirstRepro_means, variable == "IDagefirstweaning_zero"), 
    mapping = aes(x = as.factor(par1), label = labels), nudge_x = 0.55, size = 1.7) +
  expand_limits(x = 3.8)
##

ggsave(filename = "Figure5.pdf",
  # ggsave(filename = "Figure5.png",
  plot = Figure5_plusText, width = 5, height = 3)
##
