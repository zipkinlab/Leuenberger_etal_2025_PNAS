#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combining abundance, richness, diversity for figures and tables
# Using spAbundance (Jeff's package)
# Wendy Leuenberger
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load packages ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
  library(tidyverse)
  library(magrittr)
  library(spAbundance)
  library(pracma)
  library(spOccupancy)
  library(coda)
  library(abind)
  library(WendyFunctions)
  library(MetBrewer)
}

# File path
RemoteHPCC <- TRUE
if(RemoteHPCC == TRUE){
  path = 'Y:/ButterflyChapterOne/'
} else {
  path = ''
}

# Load object with dimensions 
load(file = 'Output/DataAndDimensions.RData')
# Load object with county information
load(file = paste0(path, 'Data/CountyKey.RData'))
# Source code for species traits
source('Code/R/SpeciesFunctionalTraits.R')
# Load SpeciesLists for info on which species are in which group
# Only changes if abundances change - if which species are rare/common would be 
# affected
load(file = paste0(path, 'Output/12_AbundanceTrends/SpeciesLists.RData'))

NSpp <- SpeciesLists %>% 
  lengths 
NSppGroup <- tibble(Name = names(NSpp),
                    NSpp = NSpp)
NSppGroup %<>% 
  mutate(Name = case_when(Name == 'AllSpp' ~ 'Overall', 
                          .default = Name))

# Group metrics ####

# Abundance data
load(paste0(path, 'Output/12_AbundanceTrends/Communities.RData'))
post.hoc.df.abund <- post.hoc.df
# Richness data
load(paste0(path, 'Output/10_CombinedRichness/CommunityRichness.RData'))
# Diversity data
load(paste0(path, 'Output/11_CombinedHprime/CommunityDiversity.RData'))
# Overall from below
load('Output/trend-figures/OverallMetrics.RData')


post.hoc.df.abund %<>% 
  mutate(Metric = 'Abundance',
         trend.est = trend.est * 100,
         trend.025 = trend.025 * 100,
         trend.975 = trend.975 * 100)

## Percent change for richness ####
if( exists('RichnessList') == FALSE ) {
  load(file = paste0(path, 'Output/10_CombinedRichness/RichnessList.RData'))
}
RichnessStartEnd <- tibble(Name = names(RichnessList),
                           Start = NA,
                           End = NA)
for ( nn.list in 1:length(RichnessList) ) {
  RichnessStartEnd[nn.list, 'Start'] <- 
    RichnessList[[nn.list]][ , , 1:5] %>% median
  RichnessStartEnd[nn.list, 'End'] <- 
    RichnessList[[nn.list]][ , , 28:32] %>% median
}
post.hoc.df.rich %<>% 
  filter(MeanCompare == 'Mean') %>% 
  mutate(TotalDecline = trend.est * 32) %>% 
  left_join(RichnessStartEnd) %>% 
  rename(SppPerYear.est = trend.est,
         SppPerYear.025 = trend.025,
         SppPerYear.975 = trend.975) %>% 
  # select(Name, trend.est, TotalDecline, Start) %>% 
  mutate(trend.est = SppPerYear.est/Start * 100,
         trend.025 = SppPerYear.025/Start * 100,
         trend.975 = SppPerYear.975/Start * 100) 


## Percent change for diversity ####
if (exists('DiversityList') == FALSE) {
  load(paste0(path, 'Output/11_CombinedHprime/DiversityList.RData'))
}
DiversityStartEnd <- tibble(Name = names(DiversityList),
                           Start = NA,
                           End = NA)
for ( nn.list in 1:length(DiversityList) ) {
  DiversityStartEnd[nn.list, 'Start'] <- 
    DiversityList[[nn.list]][ , , 1:5] %>% median(na.rm = TRUE)
  DiversityStartEnd[nn.list, 'End'] <- 
    DiversityList[[nn.list]][ , , 28:32] %>% median(na.rm = TRUE)
}
post.hoc.df.div %<>% 
  filter(MeanCompare == 'Mean') %>% 
  mutate(TotalDecline = trend.est * 32) %>% 
  left_join(DiversityStartEnd) %>% 
  rename(SppPerYear.est = trend.est,
         SppPerYear.025 = trend.025,
         SppPerYear.975 = trend.975) %>% 
  # select(Name, trend.est, TotalDecline, Start) %>% 
  mutate(trend.est = SppPerYear.est/Start * 100,
         trend.025 = SppPerYear.025/Start * 100,
         trend.975 = SppPerYear.975/Start * 100) %>% 
  mutate(Metric = 'Evenness')



## Overall ####
post.hoc.df.all.combinable <- post.hoc.df.all %>% 
  rename(raw.trend.est = trend.est,
         raw.trend.025 = trend.025,
         raw.trend.975 = trend.975) %>% 
  rename(trend.est = trend.est.real,
         trend.025 = trend.025.real,
         trend.975 = trend.975.real) %>% 
  mutate(MeanCompare = 'Mean',
         Name = 'Overall',
         trend.prob.neg = 1-trend.prob.pos,
         Metric = gsub(' trend', '', Metric))

# Table S1 ####
TableRichness <- post.hoc.df.rich %>% 
  bind_rows(post.hoc.df.all.combinable %>% 
              filter(Metric == 'Richness')) %>% 
  left_join(NSppGroup) %>% 
  mutate(CredInt = paste0("(", trend.025 %>% round(2), " to ",
                          trend.975 %>% round(2), ")")) %>% 
  select(Name, NSpp, trend.est, CredInt, trend.prob.neg)
# TableRichness %>%
#   slice(14, 1:8, 10:13, 9) %>%
#   write.csv(file = 'Output/10_CombinedRichness/GroupMeanRich.csv',
#                         row.names = FALSE)

# Table S2 ####
TableDiversity <- post.hoc.df.div %>% 
  bind_rows(post.hoc.df.all.combinable %>% 
              filter(Metric == 'Evenness')) %>% 
  left_join(NSppGroup) %>% 
  mutate(CredInt = paste0("(", trend.025 %>% round(2), " to ",
                          trend.975 %>% round(2), ")")) %>% 
  select(Name, NSpp, trend.est, CredInt, trend.prob.neg)
# TableDiversity %>%
#   slice(14, 1:8, 10:13, 9) %>%
#   write.csv(file = 'Output/11_CombinedHprime/GroupMeanDiv.csv',
#                         row.names = FALSE)

# Table S3 ####
TableAbundance <- post.hoc.df.abund %>% 
  filter(MeanCompare == 'Mean') %>% 
  bind_rows(post.hoc.df.all.combinable %>% 
              filter(Metric == 'Abundance')) %>% 
  left_join(NSppGroup) %>% 
  mutate(CredInt = paste0("(", trend.025 %>% round(2), " to ",
                          trend.975 %>% round(2), ")")) %>% 
  select(Name, NSpp, trend.est, CredInt, trend.prob.neg)
# TableAbundance %>%
#   slice(14, 1:8, 10:13, 9) %>%
#   write.csv(file = 'Output/12_AbundanceTrends/GroupMean.csv',
#                         row.names = FALSE)




post.hoc.df <- post.hoc.df.abund %>% 
  bind_rows(post.hoc.df.rich) %>% 
  bind_rows(post.hoc.df.div) %>% 
  bind_rows(post.hoc.df.all.combinable) %>% 
  select(trend.est, trend.prob.pos, trend.prob.neg, trend.025, trend.975,
         Name, MeanCompare, Metric)

Comparisons <- tibble(Name = c("Migratory", "Residents",
                                 "Generalists", "Specialists", 
                                 "Multivoltine", "Univoltine",
                                 "Common", "Rare",
                                 "OverwinterAdult", "OverwinterEgg", 
                                 "OverwinterLarva", "OverwinterPL",
                                 "OverwinterPupa",
                               "Overall"),
                      Comparison = c(rep('Residency', 2),
                                     rep('Specialty', 2),
                                     rep('Voltinism', 2),
                                     rep('Prevalence', 2),
                                     rep('Overwinter', 5),
                                     "Overall"))

  
post.hoc.df %<>% 
  left_join(Comparisons)

post.hoc.df %<>% 
  mutate(Name = gsub('Overwinter', '', Name))

# post.hoc.df %<>% 
#   mutate(Subheading = ifelse(Metric == 'Abundance', 
#                              'Percent change',
#                              'Species/year change'),
#          Metric = ifelse(Metric == 'Richness',
#                          'Richness trend', 
#                          'Abundance'))

post.hoc.df %>% 
  filter(MeanCompare == 'Mean',
         Metric == 'Richness') %>% 
  select(Name, trend.est, trend.prob.neg) %>% 
  arrange(trend.est) 

post.hoc.df %>% 
  filter(MeanCompare == 'Mean') %>% 
  # select(-trend.25, -trend.75, -MeanCompare) %>% 
  print(digits = 3)

post.hoc.df.mean <- post.hoc.df %>% filter(MeanCompare == 'Mean')

post.hoc.df.mean %<>%
  mutate(Name = case_when(Name == 'PL' ~ 'Larva/Pupa',
                          .default = Name),
         Name = factor(Name, levels = c("Migratory", "Residents", 
                                        "Generalists", "Specialists",
                                        "Multivoltine", "Univoltine",
                                        "Common", "Rare",
                                        "Egg", "Larva", 
                                        "Larva/Pupa", "Pupa",
                                        "Adult",
                                        "Overall")))
post.hoc.df.mean$Metric %<>% factor(levels = c("Richness", "Evenness", "Abundance"))
post.hoc.df.mean$Comparison %<>% 
  factor(levels = c("Overall", 
                    "Prevalence",
                    "Overwinter",
                    "Residency",
                    "Voltinism",
                    "Specialty"))

colors <- met.brewer('Benedictus', n = 8)


ggplot(post.hoc.df.mean,
# SaveForMike <- ggplot(post.hoc.df.mean,
# ggplot(post.hoc.df.mean %>% filter(Comparison %in% c('Overall', 'Abundance', 'Overwinter')),
# ggplot(post.hoc.df.mean %>% filter(Comparison %in% c('Residency', 'Voltinism', 'Specialty')),
       aes(x = Name, y = trend.est, color = trend.prob.neg)) +
  scale_color_gradient2(midpoint = 0.5,
                        low = colors[8], mid = 'grey', high = colors[1],
                        limits = c(0,1)) +
  geom_point() +
  geom_linerange(aes(ymin = trend.025, ymax = trend.975)) + 
  facet_grid(Metric ~ Comparison, scales = 'free', space = 'free_x') +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, # rotate axis
                                   vjust = 0.25, # line up with tick mark
                                   hjust = 1), # right align) 
        legend.position = 'bottom') +
  labs(y = 'Annual percent change (95% CI)',
       x = NULL,
       color = 'Probability of decline') +
  geom_hline(yintercept = 0)
ggsave(filename = 'Output/trend-figures/AbundanceRichnessDiversity.jpg',
       width = 180, height = 120, units = 'mm')
ggsave(filename = 'Output/trend-figures/AbundanceRichnessDiversity2.jpg',
       width = 180, height = 120, units = 'mm')
ggsave(filename = 'Output/trend-figures/AbundanceRichnessDiversity3.jpg',
       width = 180, height = 120, units = 'mm')

# Save data and plot for Mike
save(list = c('post.hoc.df.mean', 'colors', 'SaveForMike'), file = 'Output/trend-figures/CommunityFigData.RData')
  
# Vertical
# ggplot(post.hoc.df.mean,
# SaveForMike <- ggplot(post.hoc.df.mean,
ggplot(post.hoc.df.mean %>% filter(Comparison %in% c('Overall', 'Prevalence', 'Overwinter')),
# ggplot(post.hoc.df.mean %>% filter(Comparison %in% c('Residency', 'Voltinism', 'Specialty')),
  aes(y = Name, x = trend.est, color = trend.prob.neg)) +
  scale_color_gradient2(midpoint = 0.5,
                        low = colors[8], mid = 'grey', high = colors[1],
                        limits = c(0,1)) +
  geom_point() +
  geom_linerange(aes(xmin = trend.025, xmax = trend.975)) + 
  facet_grid(Comparison ~ Metric, scales = 'free', space = 'free_y') +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, # rotate axis
                                   vjust = 0.25, # line up with tick mark
                                   hjust = 1), # right align) 
        legend.position = 'bottom') +
  labs(y = 'Annual percent change (95% CI)',
       x = NULL,
       color = 'Probability of decline') +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  theme_minimal()


ggplot(post.hoc.df.mean %>% filter(Comparison %in% c('Abundance')), 
       aes(y = Name, x = trend.est)) +
  geom_point(size = 2) +
  geom_linerange(aes(xmin = trend.025, xmax = trend.975), size = 1) + 
  facet_wrap( ~ Metric + Subheading, scales = 'free_x') +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank(),
        # axis.text.x = element_text(angle = 90, # rotate axis
        #                            vjust = 0.25, # line up with tick mark
        #                            hjust = 1), # right align
        legend.position = 'none') +
  labs(x = 'Median effect (95% CI)',
       y = NULL) +
  geom_vline(xintercept = 0) 
  
ggsave(filename = 'Output/trend-figures/AbundanceRareCommon.jpg',
       width = 13.3, height = 7.5, units = 'in')
# ggsave(filename = 'Output/trend-figures/AbundanceRichness.jpg',
#        width = 180, height = 100, units = 'mm')


# All species combined ####
# Abundance
load(file = paste0(path, 'Output/12_AbundanceTrends/AllAbundance.RData'))
# Richness
load(file = paste(path, 'Output/trend-figures/richness-posthoc-for-figure.rda',
                  sep = ''))
out.posthoc.rich <- out.posthoc
# Diversity
load(file = paste(path, 'Output/trend-figures/Hprime-posthoc-for-figure.rda',
                  sep = ''))
out.posthoc.div <- out.posthoc
rm(new.data.list)

RichnessBeta <- out.posthoc.rich$beta.samples[,'years']
DiversityBeta <- out.posthoc.div$beta.samples[,'years']

post.hoc.df.all <- post.hoc.df.all.abund %>% 
  add_row(trend.est = RichnessBeta %>% quantile(0.5, names = FALSE),
          trend.prob.pos = mean(RichnessBeta > 0),
          trend.025 = RichnessBeta %>% quantile(0.025, names = FALSE),
          trend.25 = RichnessBeta %>% quantile(0.25, names = FALSE),
          trend.75 = RichnessBeta %>% quantile(0.75, names = FALSE),
          trend.975 = RichnessBeta %>% quantile(0.975, names = FALSE),
          Name = 'All',
          Metric = 'Richness trend') %>% 
  add_row(trend.est = DiversityBeta %>% quantile(0.5, names = FALSE),
          trend.prob.pos = mean(DiversityBeta > 0),
          trend.025 = DiversityBeta %>% quantile(0.025, names = FALSE),
          trend.25 = DiversityBeta %>% quantile(0.25, names = FALSE),
          trend.75 = DiversityBeta %>% quantile(0.75, names = FALSE),
          trend.975 = DiversityBeta %>% quantile(0.975, names = FALSE),
          Name = 'All',
          Metric = 'Evenness trend')

# Currently coefficients are on log scale
# Put on real scale
# Multiply by 100 for percents (for graphing)
post.hoc.df.all %<>% 
  mutate(
    trend.025.real = (exp(trend.025) - 1) * 100,
    trend.est.real = (exp(trend.est) - 1) * 100,
    trend.975.real = (exp(trend.975) - 1) * 100,
    MetricSimple = str_remove(Metric, ' trend') %>% 
      factor(levels = c('Abundance', 'Richness', 'Evenness'))
  )

post.hoc.df.all %>% 
  mutate(TotalChange = (1+(exp(trend.est) - 1))**32)

# Save trend for all for figure above
# save(post.hoc.df.all, file = 'Output/trend-figures/OverallMetrics.RData')

ggplot(post.hoc.df.all,
       # aes(x = MetricSimple, y = trend.est.real)) +
       aes(x = trend.est.real, y = MetricSimple)) +
  geom_point(size = 2) +
  scale_y_discrete(limits = rev) +
  geom_linerange(aes(xmin = trend.025.real,
                     xmax = trend.975.real),
                 size = 1) +
  # geom_linerange(aes(ymin = trend.025.real,
  #                    ymax = trend.975.real)) +
  # theme_bw(base_size = 12) +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank()) +
  # geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(x = NULL,
       y = 'Annual percent change (95% CI)')
ggsave(filename = 'Output/trend-figures/AllSppMetrics.pdf',
       width = 90, height = 75, units = 'mm')
ggsave(filename = 'Output/trend-figures/AllSppMetricsVertical.pdf',
       width = 90, height = 75, units = 'mm')
ggsave(filename = 'Output/trend-figures/AllSppMetricsVertical.jpg',
       width = 13.3, height = 7.5, units = 'in')
# ggsave(filename = 'Output/trend-figures/AllSppMetrics.jpg',
#        width = 90, height = 75, units = 'mm')

richBounds <- c(-2, 0.5)
evenBounds <- c(-2, 6)
abundBounds <- c(-4, 0.5)

# Mike's improved figure! ####
p <- post.hoc.df.mean %>% 
  ggplot() +
  theme_void() +
  geom_vline(xintercept = 0, color = "grey30", linetype = "dotted") +
  geom_hline(yintercept = 0.3, color = "grey90" ) +
  aes(y=Name, color = trend.prob.neg) +
  geom_point(aes(x=trend.est), size = 1) +
  geom_segment(aes(x = trend.025, xend = trend.975)) +
  scale_color_gradient2(midpoint = 0.5, mid = "grey50", low = colors[8], high = colors[1]) +
  labs(color = 'Probability\nof decline',
       x = 'Annual percent change (95% CI)') +
  facet_grid(
    Comparison ~ Metric, 
    scales = "free",
    switch = "y",
    space = "free_y") +
  ggh4x::facetted_pos_scales(
    x = list(
      Metric == "Richness\ntrend"  ~ scale_x_continuous(limits = c(richBounds[1],richBounds[2]),
                                                        breaks = c(-2,-1,0), expand = c(0.1,0.2)),
      Metric == "Evenness\ntrend"  ~ scale_x_continuous(limits = c(evenBounds[1],evenBounds[2]),
                                                        expand = c(0.1,0.2)),
      Metric == "Abundance\ntrend" ~ scale_x_continuous(limits = c(abundBounds[1], abundBounds[2]),
                                                        expand = c(0.1,0.2))
    )
  ) +
  scale_y_discrete() +
  coord_cartesian(clip = 'off') +
  theme(
    strip.placement = "outside",
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(angle = 90, size = 10, vjust = 1),
    # panel.border = element_rect(color = "blue", fill = NA),
    panel.spacing.x=unit(1, "lines"),
    panel.spacing.y=unit(0, "lines"),
    axis.title.x = element_text(vjust = -2.5),
    axis.text.x = element_text(vjust = -2.5),
    legend.key.size = unit(0.5, 'cm'),
    #legend.position = c(0.1, 0.6),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    plot.margin = margin(t = 0,
                         r = 0,
                         b = 9,
                         l = 2,
                         "pt")
  ) 

p

ggsave(filename = 'Output/trend-figures/PNAS_Figure2.eps',
       device = 'eps',
       plot = p,
       width = 178, height = 208, units = 'mm')

ggsave(filename = 'Output/trend-figures/PNAS_Figure2.jpg',
       device = 'jpg',
       plot = p,
       width = 178, height = 208, units = 'mm')

ggsave(filename = "Output/trend-figures/CommunityMetrics.jpg", plot = p, width = 6, height = 7)

#' i would like to have the annual percent change (95% CI) in the saveable plot
#' I would also like to have breaks on x-axis to signify where the richness trend ends 
#' and the evenness trend begins
