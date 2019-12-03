ageclass_labels <- c('all' = 'Proportion starved\n(all ages)', 'calf' = 'Calves', 'agemin10' = 'Weaning until\nage 10 yrs', 'ageplus10' = 'Age > 10 yrs')

##
head(popbif.DIST_RESAMP.R0.deathcause)
df.spread <- spread(
  subset(popbif.DIST_RESAMP.R0.deathcause, dist_resamp != "Winter_0"), 
  deathcause, value)
names(df.spread)[4:6] <- paste0("deathcause", c(3:1))
df.spread$deathcauseall = with(df.spread, deathcause1 + deathcause2 + deathcause3)
df.spread$deathcausestarved = with(df.spread, deathcause2 + deathcause3)
df.spread$combined <- df.spread$deathcauseall
df.spread$combined[df.spread$ageclass == "all"] <- df.spread$deathcausestarved[df.spread$ageclass == "all"]
df.spread$par1 <- factor(df.spread$par1, levels = seq(40,0,by=-2))
df.spread$ageclass <- factor(df.spread$ageclass, levels = c('calf', 'agemin10', 'ageplus10', 'all'))

##
Figure3 <-
  ggplot(data = df.spread, 
    mapping = aes(y = as.factor(par1), x = combined, colour = ageclass, label = round(combined,2))) +
  geom_segment(data = subset(df.spread, ageclass != "all"),
    aes(x = 0, y = as.factor(par1), xend = deathcause1, yend = as.factor(par1)), colour = GMcolours$class1, size = 1.5) +
  geom_segment(data = subset(df.spread, ageclass != "all"),
    aes(x = deathcause1, y = as.factor(par1), xend = deathcause1 + deathcause2, yend = as.factor(par1)), colour = GMcolours$class2, size = 1.5) +
  geom_segment(data = subset(df.spread, ageclass != "all"),
    aes(x = deathcause1 + deathcause2, y = as.factor(par1), xend = deathcauseall, yend = as.factor(par1)), colour = GMcolours$class3, size = 1.5) +
  geom_segment(data = subset(df.spread, ageclass == "all"),
    aes(x = 0, y = as.factor(par1), xend = combined, yend = as.factor(par1)), colour = "black", size = 1.5) +
  # geom_text(nudge_x = 0.2, size = 3, colour = 'black') +
  facet_grid(dist_resamp ~ ageclass, labeller = labeller(dist_resamp = labels_dist_resamp, ageclass = ageclass_labels)) +
  # geom_point(size = 1.5, colour = 'black') +
  xlab('') +
  scale_y_discrete("Disturbance period", seq(40,0,by=-4)) +
  coord_cartesian(xlim = c(0,0.6)) +
  ## to create legend only:
  geom_line(data = subset(df.spread, ageclass != "all"), size = 0) +
  scale_colour_manual(name = element_blank(),
    values = c('calf' = "#00AFBB", 'agemin10' = "#E7B800", 'ageplus10' = "#FC4E07"),
    labels = c('Not reduced', 'Reduced & not starving', 'Reduced & starving')) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  ggopts + ggopts_lollipop +
  theme(legend.box.margin = margin(-15,0,5,0),
    axis.text = element_text(size = 7))
##
ggsave(filename = 'Figure3.pdf',
  # ggsave(filename = 'Figure3.png',
  Figure3, heigh = 4.5, width = 6)
##
