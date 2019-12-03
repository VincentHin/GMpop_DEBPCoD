## PREPARE DATA FOR GGPLOTTING
avg <- list(popbif.DISTs_RESAMP0$avg.out, popbif.DISTs_RESAMP0.25$avg.out, popbif.DISTw_RESAMP0.25$avg.out)
avg[[1]]$dist_resamp <- 'Summer_0'
avg[[2]]$dist_resamp <- 'Summer_0.25'
avg[[3]]$dist_resamp <- 'Winter_0.25'

## create df for avg.out
whales <- ldply(avg, function(x) { subset(x = x, select = c('bifpar','total','female','male','dist_resamp')) })
whales.gath <- gather(whales, class, avg, 2:4)
whales.gath$species <- 'whales'
resource <- ldply(avg, function(x) { subset(x = x, select = c('bifpar','resource','dist_resamp')) })
resource.gath <- gather(resource, class, avg, 2)
resource.gath$species <- 'resource'

minmax <- list(popbif.DISTs_RESAMP0$minmax.out, popbif.DISTs_RESAMP0.25$minmax.out, popbif.DISTw_RESAMP0.25$minmax.out)
minmax[[1]]$dist_resamp <- 'Summer_0'
minmax[[2]]$dist_resamp <- 'Summer_0.25'
minmax[[3]]$dist_resamp <- 'Winter_0.25'

## create dfs for minmax.out
min.whales <- ldply(minmax, function(x) { subset(x = x, stat == 'min', select = c('bifpar','total','female','male','dist_resamp')) })
min.whales.gath <- gather(min.whales, class, min, 2:4)
max.whales <- ldply(minmax, function(x) { subset(x = x, stat == 'max', select = c('bifpar','total','female','male','dist_resamp')) })
max.whales.gath <- gather(max.whales, class, max, 2:4)

min.resource <-  ldply(minmax, function(x) { subset(x = x, stat == 'min', select = c('bifpar','resource','dist_resamp')) })
min.resource.gath <- gather(min.resource, class, min, 2)
max.resource <-  ldply(minmax, function(x) { subset(x = x, stat == 'max', select = c('bifpar','resource','dist_resamp')) })
max.resource.gath <- gather(max.resource, class, max, 2)

# cbind minmax columns to avg df
whales.gath$min <- min.whales.gath$min
whales.gath$max <- max.whales.gath$max
resource.gath$min <- min.resource.gath$min
resource.gath$max <- max.resource.gath$max

# rbind whales and resource data
popbif.DIST_RESAMP <- rbind(whales.gath, resource.gath)
popbif.DIST_RESAMP$species <- factor(popbif.DIST_RESAMP$species, levels = c('whales','resource'))
head(popbif.DIST_RESAMP)
## Clean intermediate objects
rm(avg,minmax,whales.gath,resource.gath,min.whales.gath,max.whales.gath,min.whales,max.whales,
  min.resource,max.resource,min.resource.gath,max.resource.gath)
##

##
vars <- c('total', 'male', 'female','resource')
name <- c('Total', 'Male', 'Female','Prey')
cols <- c('total'=GMcolours$total,
  'male'=GMcolours$male,
  'female'=GMcolours$female,
  'resource'=GMcolours$resource)
##
labels_dist_resamp <- c(Summer_0 = 'No seasonality',
                        Summer_0.25 = 'Seasonality = 0.25\nSummer disturbance',
                        Winter_0.25 = 'Seasonality = 0.25\nWinter disturbance')

Figure2 <-
  ggplot.Density.DIST_RESAMP <-
  ggplot(data = popbif.DIST_RESAMP,
    mapping = aes(bifpar, avg)) + 
  geom_line(aes(colour = class), size = 0.5) + geom_point(aes(colour = class), size = 0.5) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = class), alpha = 0.3) +
  scale_fill_manual(values=cols, breaks=vars, labels=name, name=element_blank()) +
  scale_colour_manual(values=cols, breaks=vars, labels=name, name=element_blank()) +
  guides(colour=guide_legend(override.aes = list(size = 1))) +
  facet_grid(species ~ dist_resamp, scales = 'free_y', labeller = labeller(dist_resamp = labels_dist_resamp)) +
  xlab('Disturbance duration (days)') + ylab('Mean (min/max) density') +
  ggopts +
  theme(strip.text.y = element_blank(), strip.background.y = element_blank())
##
ggsave(filename = 'Figure2.pdf',
  # ggsave(filename = 'Figure2.png',
  plot = Figure2, width = 6, height = 4.5)
##