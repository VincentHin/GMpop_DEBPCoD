## Figure 6
plot.df.melted <- melt(popbif.DIST_RESAMP.R0.stats, id.vars = c('dist_resamp','par1','no.fems'))
plot.df.melted$stat <- as.factor(ldply(strsplit(x = as.character(plot.df.melted$variable), ".", fixed = TRUE), '[', 1)[,1])
plot.df.melted$variable <- as.factor(ldply(strsplit(x = as.character(plot.df.melted$variable), ".", fixed = TRUE), '[', 2)[,1])
plot.df.melted <- spread(plot.df.melted, stat, value)
head(plot.df.melted); unique(plot.df.melted$variable); unique(plot.df.melted$dist_resamp)
##

pop_labels <- c('ageatdeath' = 'Female age at death',
                'agefirstrepro' = 'Female AfR',
                'agefirstweaning' = 'Female AfW',
                'totalfemcalves_oldfems' = 'LRO females age > 10 yrs')
#
let_labels <- data.frame(x = 1.5,
                         y = c(10.2,7.2,13,3.35),
                         letter = letters[1:4],
                         variable = c('ageatdeath', 'agefirstrepro', 'agefirstweaning', 'totalfemcalves_oldfems'))

Figure6 <-
  ggplot(data = subset(plot.df.melted, dist_resamp != 'Winter_0' & variable != 'totalfemcalves' & variable != "weanedcalves" & variable != "totalcalves_eol"),
    mapping = aes(x = par1, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = dist_resamp), alpha = 0.2) +
  geom_line(aes(col = dist_resamp), alpha = 0.8, size = 0.5) + 
  geom_point(aes(col = dist_resamp), alpha = 0.8, size = 0.5) +
  scale_colour_manual(name = element_blank(),
    values = c('Summer_0' = GMcolours$class3,
      'Summer_0.25' = GMcolours$class2,
      'Winter_0.25' = GMcolours$class1),
    labels = c("No seasonality",
      "Summer disturbance &\n seasonality = 0.25",
      "Winter disturbance &\n seasonality = 0.25")) +
  scale_fill_manual(name = element_blank(),
    values = c('Summer_0' = GMcolours$class3,
      'Summer_0.25' = GMcolours$class2,
      'Winter_0.25' = GMcolours$class1),
    labels = c("No seasonality",
      "Summer disturbance &\n seasonality = 0.25",
      "Winter disturbance &\n seasonality = 0.25")) +
  geom_text(data = let_labels,
    mapping = aes(x = x, y = y, label = letter)) +
  facet_wrap(. ~ variable, scales = 'free', ncol = 1, labeller = labeller(variable = pop_labels)) +
  guides(col = guide_legend(nrow = 3)) +
  xlab('Disturbance duration (days)') + ylab('Mean (+/- 95% CI)') +
  ggopts + theme(legend.justification = 'left', legend.key.height = unit(0.5,"cm"), legend.key.width = unit(0.5,"cm"))
##
ggsave(filename = 'Figure6.png',
  # ggsave(filename = 'Figure6.pdf',
  plot = Figure6,
  width = 2.5, height = 7, dpi = 600)
##

## PLOT LRO OF ALL FEMALES
pop_labels <- c('totalfemcalves' = 'LRO all females')
##
Figure_S3 <-
  ggplot(data = subset(plot.df.melted, dist_resamp != 'Winter_0' & variable == 'totalfemcalves'),
    mapping = aes(x = par1, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = dist_resamp), alpha = 0.2) +
  geom_line(aes(col = dist_resamp), alpha = 0.8, size = 0.5) + 
  geom_point(aes(col = dist_resamp), alpha = 0.8, size = 0.5) +
  scale_colour_manual(name = element_blank(),
    values = c('Summer_0' = GMcolours$class3,
      'Summer_0.25' = GMcolours$class2,
      'Winter_0.25' = GMcolours$class1),
    labels = c("No seasonality",
      "Summer disturbance &\n seasonality = 0.25",
      "Winter disturbance &\n seasonality = 0.25")) +
  scale_fill_manual(name = element_blank(),
    values = c('Summer_0' = GMcolours$class3,
      'Summer_0.25' = GMcolours$class2,
      'Winter_0.25' = GMcolours$class1),
    labels = c("No seasonality",
      "Summer disturbance &\n seasonality = 0.25",
      "Winter disturbance &\n seasonality = 0.25")) +
  facet_wrap(. ~ variable, scales = 'free', ncol = 1, labeller = labeller(variable = pop_labels)) +
  guides(col = guide_legend(nrow = 3)) +
  xlab('Disturbance duration (days)') + ylab('Mean (+/- 95% CI)') +
  ggopts + theme(legend.justification = 'left', legend.key.height = unit(0.5,"cm"), legend.key.width = unit(0.5,"cm"))
##
# ggsave(filename = 'Figure_S3.png',
ggsave(filename = 'Figure_S3.pdf',
  plot = Figure_S3,
  width = 2.5, height = 3, dpi = 600)
