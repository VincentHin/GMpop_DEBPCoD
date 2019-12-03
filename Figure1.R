whales.df <- subset(GMpop_default$out, select = c('time_yrs', 'total', 'male', 'female', 'resting', 'waiting', 'pregnant', 'lactating', 'waitlact', 'preglact'))
whales.df$species <- factor(x = 'Whales', levels = c('Prey', 'Whales'))
whales.melt <- melt(data = whales.df, id.vars = c('time_yrs', 'species'))
head(whales.melt); tail(whales.melt)

resource.df <- subset(GMpop_default$out, select = c('time_yrs', 'resource'))
resource.df$species <- factor(x = 'Prey', levels = c('Prey', 'Whales'))
resource.melt <- melt(data = resource.df, id.vars = c('time_yrs', 'species'))
head(resource.melt); head(resource.melt)

out.melt <- rbind(whales.melt, resource.melt)
rm(whales.df, whales.melt, resource.df, resource.melt)

vars <- c('total', 'male', 'female', 'resting', 'waiting', 'pregnant', 'lactating', 'waitlact', 'preglact')
name <- c('Total', 'Male', 'Female', 'Resting', 'Waiting', 'Pregnant', 'Lactating', 'Waiting & Lactating', 'Pregnant & Lactating')
cols <- c('total' =  GMcolours$total,
  'resource' = GMcolours$resource,
  'male' = GMcolours$male,
  'female' = GMcolours$female,
  'resting' = GMcolours$resting_trans,
  'waiting' = GMcolours$waiting_trans,
  'pregnant' = GMcolours$pregnant_trans,
  'lactating' = GMcolours$lactating_trans,
  'waitlact' = GMcolours$waitlact_trans,
  'preglact' = GMcolours$preglact_trans)
##

Figure1_timeseries <-
  ggplot(data = out.melt,
    mapping = aes(x = time_yrs, y = value, col = variable)) +
  geom_line(size = 0.2) +
  scale_colour_manual(values = cols, breaks = vars, labels = name, name = element_blank()) +
  facet_grid(species ~ ., scales = 'free_y') + 
  xlab('Time (years)') + ylab('Density') +
  xlim(c(0,200)) +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3)) +
  ggopts
##
ggsave(filename = 'Figure1_timeseries.png', plot = Figure1_timeseries, width = 4, height = 6)
##