.libPaths('~/R/x86_64-pc-linux-gnu-library/3.6')
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)

mutant23 <- list.dirs("/hosts/linuxhome/mutant23/tmp/bramve/", full.names = TRUE, recursive = FALSE)
evolving_mutation <- grep('evolving', mutant23, value = TRUE)

mutant7 <- list.dirs("/hosts/linuxhome/mutant7/tmp/bramve/", full.names = TRUE, recursive = FALSE)
evolving_production <- grep('new_seeding', mutant7, value = TRUE)

popfile_em_vector <- c()
popfile_nem_vector <- c()
costfile_em_vector <- c()
costfile_nem_vector <- c()

timepointlist <- c()
cycle_points <- c("2500","5000","7500")
for(i in cycle_points){
  timepoints <- seq(paste0(8600,as.numeric(i)),paste0(9500,as.numeric(i)),by = 1000000)
  timepointlist <- c(timepointlist, timepoints)
}

for(i in evolving_production){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_type",as.character(timepoints),".txt")))
      popfile_nem_vector <- c(popfile_nem_vector,file_name)
    }
  }
}

for(i in evolving_production){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_costs",as.character(timepoints),".txt")))
      costfile_nem_vector <- c(costfile_nem_vector, file_name)
    }
  }
}

for(i in evolving_mutation){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_type",as.character(timepoints),".txt")))
      popfile_em_vector <- c(popfile_em_vector, file_name)
    }
  }
}


for(i in evolving_mutation){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_costs",as.character(timepoints),".txt")))
      costfile_em_vector <- c(costfile_em_vector, file_name)
    }
  }
}

pop_name_vector <- c('time', 'seed','type', 'popsize_wt','popsize_mu', 'run')
cost_name_vector <- c('time', 'seed', 'type', 'cost','run')

popframe_em <- data.frame(time=numeric(0),seed=character(0),type=character(0),popsize_wt=numeric(0),popsize_mu=numeric(0),run=character(0))

for(i in popfile_em_vector){
  temp <- fread(i)
  colnames(temp) <- c('type', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*type|\\.txt$', '', head(i)))
  popsize_wt <- sum(temp$type == "1")
  popsize_mu <- sum(temp$type == "2")
  type <- as.character('evolving mutation')
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],type,popsize_wt,popsize_mu,run[2])
  colnames(run_values) <- pop_name_vector
  popframe_em <- bind_rows(popframe_em,run_values)
}

popframe_nem <- data.frame(time=numeric(0),seed=character(0),type=character(0),popsize_wt=numeric(0),popsize_mu=numeric(0),run=character(0))

for(i in popfile_nem_vector){
  temp <- fread(i)
  colnames(temp) <- c('type', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*type|\\.txt$', '', head(i)))
  popsize_wt <- sum(temp$type == "1")
  popsize_mu <- sum(temp$type == "2")
  type <- as.character("no evolving mutation")
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],type,popsize_wt,popsize_mu,run[2])
  colnames(run_values) <- pop_name_vector
  popframe_nem <- bind_rows(popframe_nem,run_values)
}

popframe <- bind_rows(popframe_em, popframe_nem)

costframe_em <- data.frame(time=numeric(0), seed = character(0), type=character(0),cost=numeric(0),run=character(0))

for(i in costfile_em_vector){
  temp <- fread(i)
  colnames(temp) <- c('cost', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*costs|\\.txt$', '', head(i)))
  cost <- mean(temp$cost[temp$cost!=0])
  type <- as.character('evolving mutation')
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],type,cost,run[2])
  colnames(run_values) <- cost_name_vector
  costframe_em <- bind_rows(costframe_em,run_values)
}

costframe_nem <- data.frame(time=numeric(0),seed=character(0), type=character(0),cost=numeric(0),run=character(0))

for(i in costfile_nem_vector){
  temp <- fread(i)
  colnames(temp) <- c('cost', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*costs|\\.txt$', '', head(i)))
  cost <- mean(temp$cost[temp$cost!=0])
  type <- as.character('no evolving mutation')
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint, seed[2], type,cost,run[2])
  colnames(run_values) <- cost_name_vector
  costframe_nem <- bind_rows(costframe_nem,run_values)
}

costframe <- bind_rows(costframe_em,costframe_nem)

pop_cost_frame <- merge(popframe,costframe)

long_pop_cost_frame <- gather(pop_cost_frame,key = 'wt_or_mu', value = 'popsize',c(popsize_wt,popsize_mu))

long_pop_cost_frame %>%
  filter(time %in% seq(86002500, 95002500, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = cost, color = interaction(wt_or_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('cost vs popsize quarter cycle') +
  theme_bw() +
    
ggsave('~/Documents/strepto/scatter_popsize_cost_em_nem/scatter_2500_revisited.png',
      height = 210,
      width  = 297,
      units= 'mm',
      dpi= 150)

long_pop_cost_frame %>%
  filter(time %in% seq(86005000, 95005000, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = cost, color = interaction(wt_or_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label=seed, hjust = -1))+
  geom_smooth(method = "lm") +
  ggtitle('cost vs popsize half cycle') +
  theme_bw() +

ggsave('~/Documents/strepto/scatter_popsize_cost_em_nem/scatter_5000_revisited.png',
       height = 210,
       width  = 297,
       units= 'mm',
       dpi= 150)

long_pop_cost_frame %>%
  filter(time %in% seq(86007500, 95007500, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = cost, color = interaction(wt_or_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('cost vs popsize three-quarter cycle') +
  theme_bw() +

ggsave('~/Documents/strepto/scatter_popsize_cost_em_nem/scatter_7500_revisited.png',
       height = 210,
       width  = 297,
       units= 'mm',
       dpi= 150)

  
  