library(tidyverse)
library(ggridges)
library(ggpmisc)
#library(bbplot)

# load data
df = read.csv('ObsPredv3.csv')
head(df)

# Figures ----

# Figure 1
# plotting density ridges between observed and predicted values
obs_pred_dist = ggplot(df)+
    labs(subtitle = '',x = 'Predicted and Observed values', y = 'Model')+
    geom_density_ridges(aes(x=Obs,y=as.factor(Dataset), fill='Obs'), alpha=0.5, scale=1)+
    geom_density_ridges(aes(x=Pred,y=as.factor(Dataset), fill='Pred'), alpha=0.5, scale=1)+
    theme_bw()+ #bbc_style()+
    theme(legend.position = c(0.9,0.85))+
    scale_fill_discrete(name = ' ')

obs_pred_dist

ggsave('Fig1-Obs_Pred_distribuitions.png',obs_pred_dist,dpi=300, height = 4, width = 3, scale = 2)

# Figure 2
# plotting density ridges between observed and predicted values
difference_dist = ggplot(df)+
    labs(subtitle = 'Density distribuition of the difference for each Model',x = 'Residuals (Observed - Predicted)', y = 'Model')+
    geom_density_ridges(aes(x=(Obs-Pred),y=as.factor(Dataset)), alpha=0.5, scale=1)+
    theme_bw()+ #bbc_style()+
    theme(legend.position = 'null')

difference_dist

ggsave('Fig2-Differences_distribuitions.png',difference_dist,dpi=300, height = 4, width = 3, scale = 2)

# Figure 3
#plotting obs vs pred
obs_pred_scatterplot = ggplot(df,aes(x=Pred,y=Obs))+
    labs(subtitle = 'Scatterplots for each Model',x = 'Predicted', y = 'Observed')+
    geom_point(aes(col = as.factor(Dataset)), size = 1)+
    geom_abline(slope = 1, intercept = 0, alpha = 0.25)+
    xlim(c(0,21))+ylim(c(0,21))+
    facet_wrap(~Dataset)+
    theme_bw()+ #bbc_style()+
    theme(legend.position = 'null')+
    stat_poly_line() +
    stat_poly_eq(aes(label = after_stat(eq.label))) +
    stat_poly_eq(label.y = 0.9)

obs_pred_scatterplot

ggsave('Fig3-Obs_Pred_scatterplot.png',obs_pred_scatterplot,dpi=300, height = 2, width = 3, scale = 2)


# Figure 4
# Residuals
residual_plot = ggplot(df)+
    labs(subtitle = '',x = 'Predicted', y = 'Residuals (Observed - Predicted)')+
    xlim(c(0,20))+
    geom_hline(aes(yintercept = 0), col = 'red')+
    geom_hline(aes(yintercept = -1.5), col = 'red', linetype=2)+
    geom_hline(aes(yintercept = 1.5), col = 'red', linetype=2)+
    geom_point(aes(x=Pred,y=(Obs-Pred), col = as.factor(Dataset)), alpha=0.8)+
    theme_bw()+ #bbc_style()+
    theme(legend.position = 'null')+
    facet_wrap(~Dataset, ncol=1)
    
residual_plot
ggsave('Fig4-Residuls.png',residual_plot,dpi=300, height = 4, width = 2, scale = 2)

# Tables ----
# Table 2
sum_table = df %>% group_by(Dataset) %>%
    summarise(obs_n = length(Obs),
            obs_min = min(Obs),
            obs_max = max(Obs),
            obs_mean = mean(Obs),
            obs_sd = sd(Obs),
            model_n = length(Pred),
            model_min = min(Pred),
            model_max = max(Pred),
            model_mean = mean(Pred),
            model_sd = sd(Pred)
            )
sum_table

# Table 3
source('Metric_calculator.R')

recomended_metrics_calc(df, Obs, Pred, Dataset)
