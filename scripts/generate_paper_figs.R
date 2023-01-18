###################################
#### Generate figures in paper ####
###################################
file.cxn <- file("logs/generate_paper_figs_LOGFILE.txt", open = "wt")
sink(file = file.cxn, type = "message")

message('START SCRIPT')
message(date(), ': Loading / installing required libraries for figure generation')
if (!require('DescTools', quietly = TRUE))
  install.packages('DescTools')
library(DescTools)

if (!require('stringr', quietly = TRUE))
  install.packages('stringr')
library(stringr)

if (!require('dplyr', quietly = TRUE))
  install.packages('dplyr')
library(dplyr)

if (!require('ggplot2', quietly = TRUE))
  install.packages('ggplot2')
library(ggplot2)

if (!require('ggExtra', quietly = TRUE))
  install.packages('ggExtra')
library(ggExtra)

message(date(), ': Generating figures in manuscript that rely on data')

## Figure 1: Study Overview
message(date(), ': Skipping Figure 1: Study Overview...')

## Figure 2: Venn Diagram
message(date(), ': Skipping Figure 2: Venn Diagram...')
## Optional - uncomment to check that the numbers reported in Venn diagram are correct
# moderator_data_dir <- '../data/results_with_moderators/'
# gw.moderator.df <- do.call(rbind, lapply(22:1, readFile <- function(x) {
#   read.csv(paste0(mediator_data_dir, x, '/all_pp_max.csv'))
# }))
# ps <- 1
# match.vec <- gw.moderator.df[[paste0('TOP_',ps)]] == gw.moderator.df[[paste0('STAB_',ps)]]
# table(match.vec)

## Figure 3: Matching vs Non-matching
message(date(), ': Generating Figure 3: Matching vs Non-matching...')
moderator_data_dir <- '../data/results_with_moderators/'
gw.moderator.df <- do.call(rbind, lapply(22:1, readFile <- function(x) {
  read.csv(paste0(mediator_data_dir, x, '/all_pp_max.csv'))
}))
gw.moderator.df$GENE <- sapply(gw.moderator.df$GENE,
                              function(x){
                                str_split(x, '_')[[1]][2]
                              })
annots_data_dir <- '../data/all_func_annots/'
ps.combined.top <- do.call(rbind,
                           lapply(paste0(annots_data_dir,
                                         1:22,'/top_annots.csv'), 
                                  read.csv))
ps.combined.stab <- do.call(rbind,
                            lapply(paste0(annots_data_dir,
                                          1:22,'/stab_annots.csv'), 
                                   read.csv))
ps <- 1 # using Potential Set 1

# Right subfigure
track <- 'CADD.PHRED'
ps.track.top <- ps.combined.top %>% 
  select(c('GENE', paste0(track,'_',ps))) %>%
  `colnames<-`(c('GENE','TOP'))
ps.track.stab <- ps.combined.stab %>% 
  select(c('GENE', paste0(track,'_',ps))) %>%
  `colnames<-`(c('GENE','STAB'))
merged.df <- merge(merge(gw.moderator.df,
                         ps.track.top,by='GENE'),
                   ps.track.stab,
                   by='GENE')
merged.df$MATCH <- (merged.df[[paste0('TOP_',ps)]] == merged.df[[paste0('STAB_',ps)]])
relevant.ids <- complete.cases(merged.df %>% select(c(paste0('TOP_',ps),paste0('STAB_',ps),'TOP','STAB')))
collect.vars.df <- merged.df[relevant.ids,] %>% select(c('TOP','STAB','MATCH'))

getProbsWithBernCI <- function(cutoff) {
  require(DescTools)
  # Compute CI for genes with matching top and stable variant
  matching.variants <- collect.vars.df %>% subset(MATCH) %>% select(c('TOP','STAB'))
  match.tf.df <- as.data.frame(matching.variants > cutoff)
  C <- nrow(match.tf.df)
  counts.vec <- table(match.tf.df$TOP) %>% as.numeric()
  #match.multinom.ci <- MultinomCI(x=counts.vec)
  match.multinom.ci <- matrix(NA,nrow=2,ncol=3)
  match.multinom.ci[,1] <- counts.vec/sum(counts.vec)
  se <- sqrt(prod(counts.vec/sum(counts.vec))/sum(counts.vec))
  match.multinom.ci[,2] <- match.multinom.ci[,1] - 1.96*se
  match.multinom.ci[,3] <- match.multinom.ci[,1] + 1.96*se
  match.prob.vec <- counts.vec / C
  
  # Compute CI for genes with non-matching top and stable variant
  nonmatching.variants <- collect.vars.df %>% subset(!MATCH) %>% select(c('TOP','STAB'))
  nonmatch.tf.df <- as.data.frame(nonmatching.variants > cutoff)
  NC <- nrow(nonmatch.tf.df)
  top.counts.vec <- table(nonmatch.tf.df$TOP) %>% as.numeric()
  stab.counts.vec <- table(nonmatch.tf.df$STAB) %>% as.numeric()
  top.multinom.ci <- matrix(NA,nrow=2,ncol=3)
  top.multinom.ci[,1] <- top.counts.vec/sum(top.counts.vec)
  se <- sqrt(prod(top.counts.vec/sum(top.counts.vec))/sum(top.counts.vec))
  top.multinom.ci[,2] <- top.multinom.ci[,1] - 1.96*se
  top.multinom.ci[,3] <- top.multinom.ci[,1] + 1.96*se
  
  stab.multinom.ci <- matrix(NA,nrow=2,ncol=3)
  stab.multinom.ci[,1] <- stab.counts.vec/sum(stab.counts.vec)
  se <- sqrt(prod(stab.counts.vec/sum(stab.counts.vec))/sum(stab.counts.vec))
  stab.multinom.ci[,2] <- stab.multinom.ci[,1] - 1.96*se
  stab.multinom.ci[,3] <- stab.multinom.ci[,1] + 1.96*se
  
  top.prob.vec <- top.counts.vec / NC
  stab.prob.vec <- stab.counts.vec / NC
  
  return(data.frame(MATCH = match.prob.vec[2],
                    MATCH.upp = match.multinom.ci[2,3],
                    MATCH.low = match.multinom.ci[2,2],
                    TOP = top.prob.vec[2],
                    TOP.upp = top.multinom.ci[2,3],
                    TOP.low = top.multinom.ci[2,2],
                    STAB = stab.prob.vec[2],
                    STAB.upp = stab.multinom.ci[2,3],
                    STAB.low = stab.multinom.ci[2,2]))
}

prob.df <- do.call(rbind, 
                   lapply(seq(from=10,to=20,length.out=200),getProbsWithBernCI))
prob.df$cutoff <- seq(from=10,to=20,length.out=200)
point.est <- prob.df[,c(1,4,7,10)] %>% reshape2::melt(id.vars = 'cutoff')
upper <- prob.df[,c(2,5,8,10)] %>%
  `colnames<-`(c('MATCH','TOP','STAB','cutoff')) %>%
  reshape2::melt(id.vars = 'cutoff')
colnames(upper)[3] <- 'upper'
lower <- prob.df[,c(3,6,9,10)] %>%
  `colnames<-`(c('MATCH','TOP','STAB','cutoff')) %>%
  reshape2::melt(id.vars = 'cutoff')
colnames(lower)[3] <- 'lower'
point.est$upper <- upper$upper; point.est$lower <- lower$lower

fig3.right.plot <- ggplot(point.est, aes(x=cutoff,y=100*value)) +
  geom_line(aes(color = variable)) +
  theme_bw() +
  ylab('Percent of All X Variants') +
  xlab(paste0('Deleteriousness cutoff')) +
  geom_vline(xintercept = 10, lty ='dashed') +
  geom_vline(xintercept = 20, lty ='dashed') +
  xlim(c(10,20))+
  ylim(c(0,15)) +
  geom_ribbon(aes(ymin = 100*lower, ymax = 100*upper, fill = variable), alpha=0.4) +
  scale_colour_manual('',labels=c('X=Matching',
                                  'X=Non-match Top',
                                  'X=Non-match Stable'),
                      values=c('#e41a1c',
                               '#4daf4a',
                               '#377eb8')) +
  scale_fill_manual('',labels=c('x=Matching',
                                'x=Non-match Top',
                                'x=Non-match Stable'),
                    values=c('#e41a1c',
                             '#4daf4a',
                             '#377eb8')) +
  guides(fill='none') +
  ggtitle('PHRED-scaled CADD') +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = 'bold'),
        legend.position=c(0.55, 0.75))

# Left subfigure
track <- 'CADD.RawScore' 
ps.track.top <- ps.combined.top %>% 
  select(c('GENE', paste0(track,'_',ps))) %>%
  `colnames<-`(c('GENE','TOP'))
ps.track.stab <- ps.combined.stab %>% 
  select(c('GENE', paste0(track,'_',ps))) %>%
  `colnames<-`(c('GENE','STAB'))
merged.df <- merge(merge(gw.moderator.df,
                         ps.track.top,by='GENE'),
                   ps.track.stab,
                   by='GENE')
merged.df$MATCH <- (merged.df[[paste0('TOP_',ps)]] == merged.df[[paste0('STAB_',ps)]])
relevant.ids <- complete.cases(merged.df %>% select(c(paste0('TOP_',ps),paste0('STAB_',ps),'TOP','STAB')))
collect.vars.df <- merged.df[relevant.ids,] %>% select(c('TOP','STAB','MATCH'))
match.scores <- (collect.vars.df %>% subset(MATCH))$STAB
nonmatch.top.scores <- (collect.vars.df %>% subset(!MATCH))$TOP
nonmatch.stab.scores <- (collect.vars.df %>% subset(!MATCH))$STAB
df.for.plot <- data.frame(SCORE = c(match.scores,
                                    nonmatch.top.scores,
                                    nonmatch.stab.scores),
                          TYPE = c(rep('Matching',length(match.scores)),
                                   rep('Non-match Top',length(nonmatch.top.scores)),
                                   rep('Non-match Stable',length(nonmatch.stab.scores))))

fig3.left.plot <- ggplot(df.for.plot, aes(x=SCORE)) +
  geom_density(aes(fill = TYPE, color = TYPE),alpha=0.45) +
  xlim(c(1,10)) +
  theme_bw() +
  scale_color_manual('',labels=c('Matching',
                                 'Non-match Stable',
                                 'Non-match Top'),
                     values=c('#e41a1c',
                              '#377eb8',
                              '#4daf4a')) +
  scale_fill_manual('',labels=c('Matching',
                                'Non-match Stable',
                                'Non-match Top'),
                    values=c('#e41a1c',
                             '#377eb8',
                             '#4daf4a')) +
  xlab('Score') +
  ylab('Density') +
  geom_vline(xintercept = 1, lty ='dashed') +
  ggtitle('Raw CADD') +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = 'bold'),
        legend.position=c(0.55, 0.75))
combined.plot <- gridExtra::grid.arrange(fig3.left.plot, fig3.right.plot, ncol = 2)  

## Uncomment to save figure in output directory
# ggsave(combined.plot,
#        file = paste0('../outputs/manuscript_figs/fig3.jpg'),
#        width = 8, height = 4,
#        dpi = 400)

## Figure 4: Top vs Stable (amongst all non-matching)
message(date(), ': Generating Figure 4: Top vs Stable...')
# Left subfigure
track <- 'CADD.RawScore' 
nonmatching.df <- merged.df %>% subset(!MATCH)
rawcadd.plot <- ggplot(nonmatching.df, aes(x=TOP,y=STAB)) +
  geom_point(alpha=0.15) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, lty= 'dashed') +
  xlab('Top Variant Score') +
  ylab('Stable Variant Score') + 
  xlim(c(-2,6))+
  ylim(c(-2,6))+
  coord_fixed() +
  ggtitle('Raw CADD') +
  theme(plot.title=element_text(face = 'bold',
                                hjust = 0.5),
        aspect.ratio=1)

fig4.left.plot <- ggExtra::ggMarginal(rawcadd.plot, type='density')

# Right subfigure
track <- 'CADD.PHRED'
ps.track.top <- ps.combined.top %>% 
  select(c('GENE', paste0(track,'_',ps))) %>%
  `colnames<-`(c('GENE','TOP'))
ps.track.stab <- ps.combined.stab %>% 
  select(c('GENE', paste0(track,'_',ps))) %>%
  `colnames<-`(c('GENE','STAB'))
merged.df <- merge(merge(gw.moderator.df,
                         ps.track.top,by='GENE'),
                   ps.track.stab,
                   by='GENE')
merged.df$MATCH <- (merged.df[[paste0('TOP_',ps)]] == merged.df[[paste0('STAB_',ps)]])
relevant.ids <- complete.cases(merged.df %>% select(c(paste0('TOP_',ps),paste0('STAB_',ps),'TOP','STAB')))
collect.vars.df <- merged.df[relevant.ids,] %>% select(c('TOP','STAB','MATCH'))

getProbsWithMultinomCI <- function(cutoff) {
  require(DescTools)
  tf.df <- as.data.frame(collect.vars.df > cutoff)
  C <- nrow(tf.df)
  count.TOP.notSTAB <- (tf.df %>% subset(TOP & !STAB) %>% nrow())
  count.STAB.notTOP <- (tf.df %>% subset(!TOP & STAB) %>% nrow())
  count.NEITHER <- (tf.df %>% subset(!TOP & !STAB) %>% nrow())
  count.BOTH <- (tf.df %>% subset(TOP & STAB) %>% nrow())
  
  multinom.ci <- MultinomCI(x=c(count.TOP.notSTAB,
                                count.STAB.notTOP,
                                count.NEITHER,
                                count.BOTH))
  prob.TOP.notSTAB <- (tf.df %>% subset(TOP & !STAB) %>% nrow()) / C
  prob.STAB.notTOP <- (tf.df %>% subset(!TOP & STAB) %>% nrow()) / C
  prob.NEITHER <- (tf.df %>% subset(!TOP & !STAB) %>% nrow()) / C
  prob.BOTH <- (tf.df %>% subset(TOP & STAB) %>% nrow()) / C
  return(data.frame(TOP.not.STAB = prob.TOP.notSTAB,
                    TOP.not.STAB.upp = multinom.ci[1,3],
                    TOP.not.STAB.low = multinom.ci[1,2],
                    STAB.not.TOP = prob.STAB.notTOP,
                    STAB.not.TOP.upp = multinom.ci[2,3],
                    STAB.not.TOP.low = multinom.ci[2,2],
                    NEITHER = prob.NEITHER,
                    NEITHER.upp = multinom.ci[3,3],
                    NEITHER.low = multinom.ci[3,2],
                    BOTH = prob.BOTH,
                    BOTH.upp = multinom.ci[4,3],
                    BOTH.low = multinom.ci[4,2]))
}
prob.df <- do.call(rbind, 
                   lapply(seq(from=10,to=20,length.out=200),getProbsWithMultinomCI))
prob.df$cutoff <- seq(from=10,to=20,length.out=200)
point.est.2 <- prob.df[,c(1,4,7,10,13)] %>% reshape2::melt(id.vars = 'cutoff')
upper <- prob.df[,c(2,5,8,11,13)] %>% 
  `colnames<-`(c('TOP.not.STAB','STAB.not.TOP','NEITHER','BOTH','cutoff')) %>%
  reshape2::melt(id.vars = 'cutoff')
colnames(upper)[3] <- 'upper'

lower <- prob.df[,c(3,6,9,12,13)] %>% 
  `colnames<-`(c('TOP.not.STAB','STAB.not.TOP','NEITHER','BOTH','cutoff')) %>%
  reshape2::melt(id.vars = 'cutoff')
colnames(lower)[3] <- 'lower'
point.est.2$upper <- upper$upper; point.est.2$lower <- lower$lower
point.est.2 <- point.est.2 %>% subset(variable != 'NEITHER')

fig4.right.plot <- ggplot(point.est.2, aes(x=cutoff,y=100*value)) +
  geom_line(aes(color = variable)) +
  theme_bw() +
  ylab('Percent of Genes') +
  xlab('Deleteriousness cutoff') +
  geom_vline(xintercept = 10, lty ='dashed') +
  geom_vline(xintercept = 20, lty ='dashed') +
  xlim(c(10,20))+
  ylim(c(0,10)) +
  geom_ribbon(aes(ymin = 100*lower, ymax = 100*upper, fill = variable), alpha=0.4) +
  scale_colour_manual('',labels=c('Only Top SNP deleterious',
                                  'Only Stable SNP deleterious',
                                  'Both deleterious'),
                      values=c('#31a354',
                               '#3182bd',
                               '#c51b8a')) +
  scale_fill_manual('',labels=c('Only Top SNP deleterious',
                                'Only Stable SNP deleterious',
                                'Both deleterious'),
                    values=c('#31a354',
                             '#3182bd',
                             '#c51b8a')) +
  
  ggtitle('PHRED-scaled CADD') +
  theme(plot.title = element_text(face = 'bold',
                                  hjust = 0.5)) +
  ggtitle('PHRED-scaled CADD') +
  theme(plot.title=element_text(face = 'bold',
                                hjust = 0.5),
        legend.key.size = unit(0.3, 'cm'),
        legend.position=c('0.6','0.75'))+
  guides(fill='none')

combined.plot <- gridExtra::grid.arrange(fig4.left.plot,
                                         fig4.right.plot,
                                         widths=c(0.5, 0.5),
                                         ncol=2)
## Uncomment to save figure in output directory
# ggsave(combined.plot,
#        file = paste0('../outputs/manuscript_figs/fig4.jpg'),
#        width = 8, height = 4,
#        dpi = 400)

## Figure 5: PICS2 Overview
message(date(), ': Skipping Figure 5: PICS2 Overview...')

message('END OF SCRIPT')
sink()