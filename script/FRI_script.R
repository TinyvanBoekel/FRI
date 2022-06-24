# This script contains the R code for the manuscript of van Boekel & Roux, Food Research International

# Libraries used:
library(papaja)
library(rmarkdown)
library(bookdown)
library(ggplot2)
library(dplyr)
library(brms)
library(tidyverse)
library(rstan)
library(broom)
library(GGally)
library(tidybayes)
library(here)
library(patchwork)
library(kableExtra)
library(ggridges)
library(cowplot)
library(cmdstanr)
theme_set(theme_bw())

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

knitr::opts_chunk$set(echo = FALSE, 
                      fig.align = "center",
                      fig.width = 6, 
                      fig.height = 4,
                      dev = "png",
                      warning = FALSE,
                      message = FALSE,
                      results= "asis",
                      tidy=TRUE,
                      comment = NA)
# Loading data

asc70 <- read.csv(file=here("data","asc_70_me.csv"), header=TRUE, sep=";")
asc70 <- rename(asc70, run=serial)
asc70$run <- as.factor(asc70$run)
asc_avg_norm <- read.csv(file=here("data", "asc_avg_norm.csv"), header=TRUE, sep=";")

# Code for Figure 1:

(asc70plot <- ggplot(data=asc70, aes(x=time, y=conc, colour=run))+
    geom_point()+
    geom_line()+
    labs(x="time (h)", y="AA (mM)"))+
  theme(legend.key.size = unit(0.5, 'cm'))

# brms code for regression of averaged normalized concentration:
nlform<-bf(avg ~ exp(-kr*time),
           kr~1, nl=TRUE)

nlprior<-c(prior(normal(0.5,0.3), nlpar="kr", lb=0),
           prior(cauchy(0,10), class="sigma")
)           

asc_avg_norm_bayes <- brm(formula=nlform, data=asc_avg_norm, family = gaussian(), 
                          prior = nlprior, iter = 4000, warmup=2000, 
                          control = list(adapt_delta = 0.99),
                          file=here("fits", "asc_avg_norm_bayes"))
asc_avg_norm_bayes_post <- as_draws_df(asc_avg_norm_bayes) # store posterior distribution in a dataframe for later use

# Code for Table 1:

asc_avg_norm_summary1 <- summary(asc_avg_norm_bayes)
asc_avg_norm_summary2 <- rbind(data.frame(asc_avg_norm_summary1$fixed), data.frame(asc_avg_norm_summary1$spec_pars) )
rownames(asc_avg_norm_summary2) <- c("$k_r (\\text {h}^{-1})$","$\\sigma_e \\text {(mM)}$")
colnames(asc_avg_norm_summary2) <- c("mean","SE", "lower bound", "upper bound")

asc_avg_norm_summary2[1:2, 1:4] %>% rownames_to_column(var = "parameter") %>% kbl(digits = c(3,3,3), escape=F, booktabs=TRUE, caption = "Numerical summary of the posterior distribution resulting from Bayesian nonlinear regression for the normalized and averaged ascorbic acid data using the first-order model. SE= standard error, the lower and upper bound represent the 95% credible interval") %>% kable_styling(position="center", full_width = F)

# Code for Figure 2A:

asc_avg_norm_cor <- dplyr::select(asc_avg_norm_bayes_post,b_kr_Intercept:sigma)
asc_avg_norm_cor <- setNames(asc_avg_norm_cor, c(c(expression(paste("k"[r])), expression(sigma[e]))))

asc_avg_norm_corplot1 <-asc_avg_norm_cor  %>% ggpairs(diag=list(continuous="densityDiag"),
                                                      mapping=aes(fill="red"),
                                                      upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
                                                      labeller=label_parsed)+ 
  theme(strip.text.x = element_text(size = 9, color = "red"),
        strip.text.y = element_text(size = 9, color = "red"))+
  theme(axis.text.x = element_text(angle=70, hjust=1)
  )

# Code for Figure 2B:

newavg <- data.frame(time = seq(from = 0, to = 7, by = 0.1))

fitavg1_norm <- cbind(newavg, fitted(asc_avg_norm_bayes, newdata = newavg, re_formula = NA)[,-2])
#[,-2] means that the second column is omitted from fitted, containing the est.error and leaving estimate and Q2.5 and Q97.5, re_formula=NA means that group-level effects are omitted, so fitting only the grand mean population parameters
names(fitavg1_norm) <- c("time", "avg", "lower", "upper")
#same for prediction intervals
predavg1_norm <- cbind(newavg, predict(asc_avg_norm_bayes, newdata=newavg, re_formula=NA)[,-2])
names(predavg1_norm) <- c("time", "avg", "lower", "upper")
names(asc_avg_norm) <- c("time", "avg", "std")

plot_fe_asc_norm <- ggplot(data=asc_avg_norm, aes(x = time, y = avg)) +
  geom_point(size=2, shape=21,stroke=1, fill="red") +
  geom_errorbar(aes(ymin=avg-std, ymax=avg+std), width=.2)+
  geom_line(data=fitavg1_norm, aes(y=avg)) +
  geom_ribbon(data = fitavg1_norm, aes(ymin = lower, ymax = upper), fill = "blue", alpha = .3)+
  geom_ribbon(data = predavg1_norm, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = .6)+
  coord_cartesian(ylim = c(0,1))+
  labs(x="time (h)", y= expression(c/c[0]))

# Code for Figure 2C, D, E:

# calculating the nonlinear regression line: 
time.seq <- data.frame(time = seq(from = 0, to = 7, by = 0.1))

#residuals, QQ and lag plots
resid_asc_avg = resid(asc_avg_norm_bayes)[, "Estimate"]
residplot_asc_avg <-  ggplot() + geom_point(data = NULL, size = 2, aes(x = asc_avg_norm$time, y = resid_asc_avg))+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  xlim(0,7)+
  labs(x = "time (min)", y = "residual")

# calculating the QQ plot for residuals:
res_asc_avg <- as.data.frame(residuals(asc_avg_norm_bayes, summary=T))
qqplot_asc_avg <-  ggplot(res_asc_avg, aes(sample = Estimate)) +
  stat_qq() +
  stat_qq_line()

# calculating the lag plot

lagx_asc_avg <- resid_asc_avg[-1]
lagy_asc_avg <- resid_asc_avg[-length((resid_asc_avg))]
lag_asc_avg <- data.frame(lagx_asc_avg,lagy_asc_avg)
lagplot_asc_avg <- ggplot(lag_asc_avg, aes(x=lagx_asc_avg, y=lagy_asc_avg))+
  geom_point()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(x="residual(i-1)",y="residual(i)")

top_row <- plot_grid(ggmatrix_gtable(asc_avg_norm_corplot1),plot_fe_asc_norm,labels = c('A', 'B'))

bottom_row <- plot_grid(residplot_asc_avg,qqplot_asc_avg,lagplot_asc_avg, ncol = 3,labels = c('C', 'D','E'))

combi_plot <- plot_grid(top_row, bottom_row, nrow=2, align="v")
combi_plot

# Code for brms regression of pooled data:

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), c0~1, k~1, nt~1,nl=TRUE)
nlprior<-c(prior(normal(6,1), nlpar = "c0"),
           prior(normal(0.3,0.3), nlpar="k"),
           prior(normal(1,0.5), nlpar="nt"),
           prior(cauchy(0,10), class="sigma")
)           

asc_bayes_all <- brm(formula=nlform, data=asc70, family = gaussian(), prior = nlprior, warmup=3000, iter=6000, chains=4, control = list(adapt_delta = 0.999), file=here("fits","asc_bayes_all"))

asc_bayes_all_post <- as_draws_df(asc_bayes_all)

# Code for Figure 3:

asc_cor <- dplyr::select(asc_bayes_all_post,b_c0_Intercept:sigma)
asc_cor <- setNames(asc_cor, c(expression(paste("c"[0])), expression(paste("k"[r])),expression(paste("n"[t])), expression(sigma[e])))

asc_corplot <-asc_cor  %>% ggpairs(diag=list(continuous="densityDiag"),
                                   mapping=aes(fill="red"),
                                   upper = list(continuous = wrap("cor", size = 3, stars=FALSE, title="corr. coef")), 
                                   labeller=label_parsed)+ 
  theme(strip.text.x = element_text(size = 8, color = "red"),
        strip.text.y = element_text(size = 8, color = "red"))+
  theme(axis.text.x = element_text(size=8, angle=70, hjust = 1))+
  theme(axis.text.y = element_text(size=6))

# plot with pooled bayesian regression line

newavg <- data.frame(time = seq(from = 0, to = 7, by = 0.1))

fitavg1 <- cbind(newavg, fitted(asc_bayes_all, newdata = newavg, re_formula = NA)[,-2])
#[,-2] means that the second column is omitted from fitted, containing the est.error and leaving estimate and Q2.5 and Q97.5, re_formula=NA means that group-level effects are omitted
names(fitavg1) <- c("time", "conc", "lower", "upper")
#same for prediction intervals
predavg1 <- cbind(newavg, predict(asc_bayes_all, newdata=newavg, re_formula=NA)[,-2])
names(predavg1) <- c("time", "conc", "lower", "upper")

plot_fe_asc <- ggplot(asc70, aes(x = time, y = conc)) +
  geom_point(size=2, shape=21, stroke=1, fill="red") +
  geom_line(data=fitavg1, aes(y=conc)) +
  geom_ribbon(data = fitavg1, aes(ymin = lower, ymax = upper), fill = "blue", alpha = .3)+
  geom_ribbon(data = predavg1, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = .6)+
  coord_cartesian(ylim = c(0,7.5))+
  labs(x="time (h)", y= "[AA] (mM)")

resid_asc1 = resid(asc_bayes_all)[, "Estimate"]

#residuals, QQ and lag plots

residplot_asc1 <-  ggplot() + geom_point(data = NULL, size = 2, aes(x = asc70$time, y = resid_asc1))+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  labs(x = "time (h)", y = "residual")

# calculating the QQ plot for residuals:
res_asc1 <- as.data.frame(residuals(asc_bayes_all, summary=T))

qqplot_asc1 <-  ggplot(res_asc1, aes(sample = Estimate)) +
  stat_qq() +
  stat_qq_line() 

#lagplot
lagx <- resid_asc1[-1]
lagy <- resid_asc1[-length((resid_asc1))]
lag <- data.frame(lagx,lagy)
lagplot_asc1 <- ggplot(lag, aes(x=lagx, y=lagy))+
  geom_point(size=1)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(x="residual(i-1)" , y="residual(i)")

top_row2 <- plot_grid(ggmatrix_gtable(asc_corplot),plot_fe_asc,labels = c('A', 'B'))

bottom_row2 <- plot_grid(residplot_asc1,qqplot_asc1,lagplot_asc1, ncol = 3,labels = c('C', 'D','E'))

combi_plot2 <- plot_grid(top_row2, bottom_row2, nrow=2, align="v")
combi_plot2

# Code for Table 2:

asc_summary1 <- summary(asc_bayes_all)
asc_summary2 <- rbind(data.frame(asc_summary1$fixed), data.frame(asc_summary1$spec_pars) )
rownames(asc_summary2) <- c("$c_0 \\text {(mM)}$", "$k_r (\\text {h}^{-1})$", "$n_t (-)$", "$\\sigma_e \\text {(mM)}$")
colnames(asc_summary2) <- c("mean","SE", "lower bound", "upper bound")

asc_summary2[1:4, 1:4] %>% rownames_to_column(var = "parameter") %>% kbl(digits = c(2,2,2,2), escape=F, booktabs=TRUE, caption = "Numerical summary of the posterior distribution resulting from nonlinear Bayesian regression of the $n^{th}$-order model for the pooled ascorbic acid data. SE= standard error, the lower and upper bound represent the 95% credible interval") %>% kable_styling(position="center", full_width = F)


#calculating individual fits for all runs 
time.seq <- data.frame(time = seq(from = 0, to = 7, by = 0.1))
#run 1

asc_run1 <- subset(asc70, run==1, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model1 <- brm(formula=nlform, data=asc_run1, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits", "asc_model1"))
post_asc1 <- as_draws_df(asc_model1)

# composing a dataframe containing fitted values for the individual regressions:

fit_asc1 <- cbind(time.seq, fitted(asc_model1, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc1$run <- as.factor(rep(1, length(nrow(time.seq))))

#predline1 <-
#  predict(asc_model1, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline1, file=here("fits","predline1"))

predline1 <- readRDS(file=here("fits", "predline1"))
predline1$run <- as.factor(rep(1, length(nrow(time.seq))))

regrline1 <- ggplot(data = asc_run1, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline1, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc1, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 1")+
  ylim(0,6.2)


#run 2

asc_run2 <- subset(asc70, run==2, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model2 <- brm(formula=nlform, data=asc_run2, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model2"))

post_asc2 <- as_draws_df(asc_model2)

fit_asc2 <- cbind(time.seq, fitted(asc_model2, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc2$run <- as.factor(rep(2, length(nrow(time.seq))))

#predline2 <-
#  predict(asc_model2, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline2, file=here("fits","predline2"))

predline2 <- readRDS(file=here("fits", "predline2"))
predline2$run <- as.factor(rep(2, length(nrow(time.seq))))

regrline2 <- ggplot(data = asc_run2, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline2, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc2, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 2")+
  ylim(0,6.2)

#run 3

asc_run3 <- subset(asc70, run==3, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model3 <- brm(formula=nlform, data=asc_run3, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model3"))
print(asc_model3)
post_asc3 <- as_draws_df(asc_model3)
fit_asc3 <- cbind(time.seq, fitted(asc_model3, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc3$run <- as.factor(rep(3, length(nrow(time.seq))))

#predline3 <-
#  predict(asc_model3, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline3, file=here("fits","predline3"))

predline3 <- readRDS(file=here("fits", "predline3"))
predline3$run <- as.factor(rep(3, length(nrow(time.seq))))

regrline3 <- ggplot(data = asc_run3, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline3, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc3, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 3")+
  ylim(0,6.2)

#run 4

asc_run4 <- subset(asc70, run==4, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model4 <- brm(formula=nlform, data=asc_run4, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model4"))

post_asc4 <- as_draws_df(asc_model4)

fit_asc4 <- cbind(time.seq, fitted(asc_model4, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc4$run <- as.factor(rep(4, length(nrow(time.seq))))

#predline4 <-
#  predict(asc_model4, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline4, file=here("fits","predline4"))

predline4 <- readRDS(file=here("fits", "predline4"))
predline4$run <- as.factor(rep(4, length(nrow(time.seq))))

regrline4 <- ggplot(data = asc_run4, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline4, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc4, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 4")+
  ylim(0,6.2)

#run 5
asc_run5 <- subset(asc70, run==5, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model5 <- brm(formula=nlform, data=asc_run5, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model5"))
post_asc5 <- as_draws_df(asc_model5)

fit_asc5 <- cbind(time.seq, fitted(asc_model5, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc5$run <- as.factor(rep(5, length(nrow(time.seq))))

#predline5 <-
#  predict(asc_model5, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline5, file=here("fits","predline5"))

predline5 <- readRDS(file=here("fits", "predline5"))
predline5$run <- as.factor(rep(5, length(nrow(time.seq))))

regrline5 <- ggplot(data = asc_run5, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline5, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc5, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 5")+
  ylim(0,6.2)

#run 6

asc_run6 <- subset(asc70, run==6, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model6 <- brm(formula=nlform, data=asc_run6, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model6"))

post_asc6 <- as_draws_df(asc_model6)

fit_asc6 <- cbind(time.seq, fitted(asc_model6, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc6$run <- as.factor(rep(6, length(nrow(time.seq))))

#predline6 <-
#  predict(asc_model6, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline6, file=here("fits","predline6"))

predline6 <- readRDS(file=here("fits", "predline6"))
predline6$run <- as.factor(rep(6, length(nrow(time.seq))))

regrline6 <- ggplot(data = asc_run6, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline6, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc6, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 6")+
  ylim(0,6.2)

#run 7

asc_run7 <- subset(asc70, run==7, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model7 <- brm(formula=nlform, data=asc_run7, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model7"))
post_asc7 <- as_draws_df(asc_model7)

fit_asc7 <- cbind(time.seq, fitted(asc_model7, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc7$run <- as.factor(rep(7, length(nrow(time.seq))))

#predline7 <-
#  predict(asc_model7, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline7, file=here("fits","predline7"))

predline7 <- readRDS(file=here("fits", "predline7"))
predline7$run <- as.factor(rep(7, length(nrow(time.seq))))

regrline7 <- ggplot(data = asc_run7, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline7, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc7, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 7")+
  ylim(0,6.2)

#run 8

asc_run8 <- subset(asc70, run==8, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model8 <- brm(formula=nlform, data=asc_run8, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model8"))
post_asc8 <- as_draws_df(asc_model8)

fit_asc8 <- cbind(time.seq, fitted(asc_model8, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc8$run <- as.factor(rep(8, length(nrow(time.seq))))

#predline8 <-
#  predict(asc_model8, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline8, file=here("fits","predline8"))

predline8 <- readRDS(file=here("fits", "predline8"))
predline8$run <- as.factor(rep(8, length(nrow(time.seq))))

regrline8 <- ggplot(data = asc_run8, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline8, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc8, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 8")+
  ylim(0,6.2)

#run 9

asc_run9 <- subset(asc70, run==9, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model9 <- brm(formula=nlform, data=asc_run9, family=gaussian(),
                  prior=nlprior, control = list(adapt_delta=0.99),
                  file=here("fits","asc_model9"))
post_asc9 <- as_draws_df(asc_model9)

fit_asc9 <- cbind(time.seq, fitted(asc_model9, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc9$run <- as.factor(rep(9, length(nrow(time.seq))))

#predline9 <-
#  predict(asc_model9, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline9, file=here("fits","predline9"))

predline9 <- readRDS(file=here("fits", "predline9"))
predline9$run <- as.factor(rep(9, length(nrow(time.seq))))

regrline9 <- ggplot(data = asc_run9, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline9, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc9, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 9")+
  ylim(0,6.2)

#run 10

asc_run10 <- subset(asc70, run==10, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model10 <- brm(formula=nlform, data=asc_run10, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","asc_model10"))

post_asc10 <- as_draws_df(asc_model10)

fit_asc10 <- cbind(time.seq, fitted(asc_model10, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc10$run <- as.factor(rep(10, length(nrow(time.seq))))

#predline10 <-
#  predict(asc_model10, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline10, file=here("fits","predline10"))

predline10 <- readRDS(file=here("fits", "predline10"))
predline10$run <- as.factor(rep(10, length(nrow(time.seq))))

regrline10 <- ggplot(data = asc_run10, 
                     aes(x = time, y = conc)) +
  geom_ribbon(data = predline10, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc10, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 10")+
  ylim(0,6.2)

#run 11

asc_run11 <- subset(asc70, run==11, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model11 <- brm(formula=nlform, data=asc_run11, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","asc_model11"))
post_asc11 <- as_draws_df(asc_model11)

fit_asc11 <- cbind(time.seq, fitted(asc_model11, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc11$run <- as.factor(rep(11, length(nrow(time.seq))))

#predline11 <-
#  predict(asc_model11, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline11, file=here("fits","predline11"))

predline11 <- readRDS(file=here("fits", "predline11"))
predline11$run <- as.factor(rep(11, length(nrow(time.seq))))

regrline11 <- ggplot(data = asc_run11, 
                     aes(x = time, y = conc)) +
  geom_ribbon(data = predline11, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc11, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 11")+
  ylim(0,6.2)

#run 12

asc_run12 <- subset(asc70, run==12, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model12 <- brm(formula=nlform, data=asc_run12, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","asc_model12"))

post_asc12 <- as_draws_df(asc_model12)

fit_asc12 <- cbind(time.seq, fitted(asc_model12, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc12$run <- as.factor(rep(12, length(nrow(time.seq))))

#predline12 <-
#  predict(asc_model12, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline12, file=here("fits","predline12"))

predline12 <- readRDS(file=here("fits", "predline12"))
predline12$run <- as.factor(rep(12, length(nrow(time.seq))))

regrline12 <- ggplot(data = asc_run12, 
                     aes(x = time, y = conc)) +
  geom_ribbon(data = predline12, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc12, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 12")+
  ylim(0,6.2)

#run 13

asc_run13 <- subset(asc70, run==13, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model13 <- brm(formula=nlform, data=asc_run13, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","asc_model13"))
post_asc13 <- as_draws_df(asc_model13)

fit_asc13 <- cbind(time.seq, fitted(asc_model13, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc13$run <- as.factor(rep(13, length(nrow(time.seq))))

#predline13 <-
#  predict(asc_model13, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline13, file=here("fits","predline13"))

predline13 <- readRDS(file=here("fits", "predline13"))
predline13$run <- as.factor(rep(13, length(nrow(time.seq))))

regrline13 <- ggplot(data = asc_run13, 
                     aes(x = time, y = conc)) +
  geom_ribbon(data = predline13, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc13, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 13")+
  ylim(0,6.2)

#run 14

asc_run14 <- subset(asc70, run==14, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model14 <- brm(formula=nlform, data=asc_run14, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","asc_model14"))
post_asc14 <- as_draws_df(asc_model14)

fit_asc14 <- cbind(time.seq, fitted(asc_model14, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc14$run <- as.factor(rep(14, length(nrow(time.seq))))

#predline14 <-
#  predict(asc_model14, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline14, file=here("fits","predline14"))

predline14 <- readRDS(file=here("fits", "predline14"))
predline14$run <- as.factor(rep(14, length(nrow(time.seq))))

regrline14 <- ggplot(data = asc_run14, 
                     aes(x = time, y = conc)) +
  geom_ribbon(data = predline14, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc14, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 14")+
  ylim(0,6.2)

#run 15

asc_run15 <- subset(asc70, run==15, time:conc)

nlform <- bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
             c0~1,
             nt~1,
             kr~1, 
             nl=TRUE)

nlprior <- c(prior(normal(6,1),nlpar="c0" ),
             prior(normal(1,1), nlpar="nt"),
             prior(normal(0.3,0.3), nlpar="kr"),
             prior(cauchy(0,1), class="sigma")
)
asc_model15 <- brm(formula=nlform, data=asc_run15, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","asc_model15"))

post_asc15 <- as_draws_df(asc_model15)

fit_asc15 <- cbind(time.seq, fitted(asc_model15, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_asc15$run <- as.factor(rep(15, length(nrow(time.seq))))

#predline15 <-
#  predict(asc_model15, 
#          newdata = time.seq) %>%
# as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline15, file=here("fits","predline15"))

predline15 <- readRDS(file=here("fits", "predline15"))
predline15$run <- as.factor(rep(15, length(nrow(time.seq))))

regrline15 <- ggplot(data = asc_run15, 
                     aes(x = time, y = conc)) +
  geom_ribbon(data = predline15, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_asc15, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="AA (mM)", subtitle = "run 15")+
  ylim(0,6.2)


# collection of unpooled individual fits:
fit_unpooled_all <- rbind(predline1, predline2, predline3, predline4, predline5, predline6, predline7, predline8, predline9, predline10, predline11, predline12, predline13, predline14, predline15)

# Code for Figure 4:

fitplot <- ggplot(asc70, aes(x=time, y=conc))+
  geom_point()+
  geom_ribbon(data=fit_unpooled_all,aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),fill = "lightblue", alpha=0.8) +
  geom_line(data=fit_unpooled_all, aes(x=time, y=Estimate))+
  facet_wrap(~run, ncol=8)+
  labs(x="time (h)", y="AA, mmol/L")

resid <- resid(asc_model1)[, "Estimate"]
resid_run1 <- data.frame(resid) %>% bind_cols(time=asc_run1$time)
resid_run1$run <- as.factor(rep(1, length(nrow(time.seq))))

resid = resid(asc_model2)[, "Estimate"]
resid_run2 <- data.frame(resid) %>% bind_cols(time=asc_run2$time)
resid_run2$run <- as.factor(rep(2, length(nrow(time.seq))))

resid = resid(asc_model3)[, "Estimate"]
resid_run3 <- data.frame(resid) %>% bind_cols(time=asc_run3$time)
resid_run3$run <- as.factor(rep(3, length(nrow(time.seq))))

resid = resid(asc_model4)[, "Estimate"]
resid_run4 <- data.frame(resid) %>% bind_cols(time=asc_run4$time)
resid_run4$run <- as.factor(rep(4, length(nrow(time.seq))))

resid = resid(asc_model5)[, "Estimate"]
resid_run5 <- data.frame(resid) %>% bind_cols(time=asc_run5$time)
resid_run5$run <- as.factor(rep(5, length(nrow(time.seq))))

resid = resid(asc_model6)[, "Estimate"]
resid_run6 <- data.frame(resid) %>% bind_cols(time=asc_run6$time)
resid_run6$run <- as.factor(rep(6, length(nrow(time.seq))))

resid = resid(asc_model7)[, "Estimate"]
resid_run7 <- data.frame(resid) %>% bind_cols(time=asc_run7$time)
resid_run7$run <- as.factor(rep(7, length(nrow(time.seq))))

resid = resid(asc_model8)[, "Estimate"]
resid_run8 <- data.frame(resid) %>% bind_cols(time=asc_run8$time)
resid_run8$run <- as.factor(rep(8, length(nrow(time.seq))))

resid = resid(asc_model9)[, "Estimate"]
resid_run9 <- data.frame(resid) %>% bind_cols(time=asc_run9$time)
resid_run9$run <- as.factor(rep(9, length(nrow(time.seq))))

resid = resid(asc_model10)[, "Estimate"]
resid_run10 <- data.frame(resid) %>% bind_cols(time=asc_run10$time)
resid_run10$run <- as.factor(rep(10, length(nrow(time.seq))))

resid = resid(asc_model11)[, "Estimate"]
resid_run11 <- data.frame(resid) %>% bind_cols(time=asc_run11$time)
resid_run11$run <- as.factor(rep(11, length(nrow(time.seq))))

resid = resid(asc_model12)[, "Estimate"]
resid_run12 <- data.frame(resid) %>% bind_cols(time=asc_run12$time)
resid_run12$run <- as.factor(rep(12, length(nrow(time.seq))))

resid = resid(asc_model13)[, "Estimate"]
resid_run13 <- data.frame(resid) %>% bind_cols(time=asc_run13$time)
resid_run13$run <- as.factor(rep(13, length(nrow(time.seq))))

resid = resid(asc_model14)[, "Estimate"]
resid_run14 <- data.frame(resid) %>% bind_cols(time=asc_run14$time)
resid_run14$run <- as.factor(rep(14, length(nrow(time.seq))))

resid = resid(asc_model15)[, "Estimate"]
resid_run15 <- data.frame(resid) %>% bind_cols(time=asc_run15$time)
resid_run15$run <- as.factor(rep(15, length(nrow(time.seq))))

resid_all <- rbind(resid_run1,resid_run2,resid_run3,resid_run4,resid_run5,resid_run6,resid_run7,resid_run8,resid_run9,resid_run10,resid_run11,resid_run12,resid_run13,resid_run14,resid_run15)

residplot <-  ggplot(data = resid_all, aes(x=time, y=resid)) + 
  geom_point()+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  xlim(0,7)+
  facet_wrap(~run, ncol = 8)+
  labs(x="time (h)", y="residual")

p_nt <-
  bind_rows(
    post_asc1,
    post_asc2,
    post_asc3,
    post_asc4,
    post_asc5,
    post_asc6,
    post_asc7,
    post_asc8,
    post_asc9,
    post_asc10,
    post_asc11,
    post_asc12,
    post_asc13,
    post_asc14,
    post_asc15
  )
iter <- 4000

p_nt <- 
  p_nt %>% 
  mutate(run = rep(c("1","2","3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"),each = iter)) %>%  mutate(run=fct_relevel(run, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))

p_c0_plot <- p_nt %>% 
  ggplot(aes(x = b_c0_Intercept, y = run)) +
  geom_halfeyeh(fill = "green4", alpha=0.5,
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("c"[0])), y="run")

p_nt_plot <- p_nt %>% 
  ggplot(aes(x = b_nt_Intercept, y = run)) +
  geom_halfeyeh(fill = "green4", alpha=0.5,
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("n"[t])), y="run")

p_kr_plot <- p_nt %>% 
  ggplot(aes(x = b_kr_Intercept, y = run)) +
  geom_halfeyeh(fill = "green4", alpha=0.5,
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("k"[r])), y="run")

parplot <- p_c0_plot+p_nt_plot+p_kr_plot 

total_plot <- fitplot/residplot/parplot+plot_annotation(tag_levels = 'A')
total_plot

# Code for brms regression of partial pooled data (multilevel modeling)

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt)), 
           c0 ~1+(1|ID|run),
           nt~1+(1|ID|run),
           kr ~ 1+(1|ID|run), 
           nl=TRUE)

nlprior<-c(prior(normal(6.0,0.5), nlpar = "c0", lb=0),
           prior(normal(0.5,0.1), nlpar="kr", lb=0),
           prior(normal(1,0.1), nlpar="nt", lb=0),
           prior(cauchy(0,1), class="sigma"),
           prior(cauchy(0,1), class="sd", nlpar="c0"),
           prior(cauchy(0,1), class="sd", nlpar="nt"),
           prior(cauchy(0,1), class="sd", nlpar="kr"),
           prior(lkj(1), class="cor")
)           

inits <- list(
  c0 = 6.0,
  kr= 0.6,
  nt=0.92
)

list_of_inits <- list(inits, inits, inits, inits)

asc_multi_bayes <- brm(formula=nlform, data=asc70, family=gaussian(),
                       prior=nlprior, iter=20000, warmup=10000, seed=1234,
                       chains=4,
                       cores=4,
                       inits = list_of_inits,
                       control = list(adapt_delta=0.9999, max_treedepth=15), backend="cmdstanr",
                       file=here("fits","asc_multi_bayes"))

asc_multi_bayes_post <- as_draws_df(asc_multi_bayes) # store posterior distribution in a dataframe for later use

# Code for Figure 5:

asc_cor2 <- dplyr::select(asc_multi_bayes_post,b_c0_Intercept:sigma)

asc_cor2 <- setNames(asc_cor2, c(expression(paste("c"[0])),expression(paste("n"[t])), expression(paste("k"[r])), expression(sigma[c[0]]),expression(sigma[n[t]]), expression(sigma[k[r]]),expression(rho[c[0]-n[t]]),expression(rho[c[0]-k[r]]),expression(rho[n[t]-k[r]]), expression(sigma[e])))

asc_corplot2 <-asc_cor2  %>% ggpairs(diag=list(continuous="densityDiag"),upper = list(continuous = wrap("cor", size = 2, stars=FALSE, title="cor.coef")), mapping=aes(fill="red"), labeller=label_parsed)+ 
  theme(text=element_text(size=7),strip.text.x = element_text(size = 7, color = "red"),
        strip.text.y = element_text(size = 7, color = "red")+
          theme(axis.text = element_text(size=1))
  )

intercept_run <- asc_multi_bayes %>% spread_draws(b_c0_Intercept, r_run__c0[run,term]) %>% 
  mutate(b_c0_Intercept = r_run__c0 + b_c0_Intercept) %>% filter(term=="Intercept")

intercept_run$run <- as.character(intercept_run$run)

# Average effect: the grand mean
intercept_average <- spread_draws(asc_multi_bayes, b_c0_Intercept) %>% 
  mutate(run = "Grand mean")

# Combine average and run-specific effects in a data frame

intercept_all <- bind_rows(intercept_run, intercept_average) %>% 
  ungroup() %>% mutate(run=fct_relevel(run,"Grand mean","1","2","3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" ))

# Data frame of summary numbers
intercept_all_sum <- group_by(intercept_all, run) %>% 
  mean_qi(b_c0_Intercept)
# forest/ridgeline plot for c0
forest_c0 <- intercept_all %>% ggplot(aes(b_c0_Intercept, run))+
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1,fill="Lightblue"
  ) +
  geom_pointintervalh(
    data = intercept_all_sum, size = 1
  ) +
  labs(x=expression(paste("c"[0])), y="run", subtitle="C")

#procedure for nt:
intercept_run <- asc_multi_bayes %>% spread_draws(b_nt_Intercept, r_run__nt[run,term]) %>% 
  mutate(b_nt_Intercept = r_run__nt + b_nt_Intercept) %>% filter(term=="Intercept")

intercept_run$run <- as.character(intercept_run$run)

# Average effect: the grand mean for nt
intercept_average <- spread_draws(asc_multi_bayes, b_nt_Intercept) %>% 
  mutate(run = "Grand mean")

# Combine average and study-specific effects data frames for nt

intercept_all <- bind_rows(intercept_run, intercept_average) %>% 
  ungroup() %>% mutate(run=fct_relevel(run,"Grand mean","1","2","3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" ))

# Data frame of summary numbers for nt
intercept_all_sum <- group_by(intercept_all, run) %>% 
  mean_qi(b_nt_Intercept)

# forest/ridgeline plot for nt:
forest_nt <- intercept_all %>% ggplot(aes(b_nt_Intercept, run))+
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1,fill="Lightblue"
  ) +
  geom_pointintervalh(
    data = intercept_all_sum, size = 1
  ) +
  labs(x=expression(paste("n"[t])), y="run", subtitle="D")


# procedure for kr:
intercept_run <- asc_multi_bayes %>% spread_draws(b_kr_Intercept, r_run__kr[run,term]) %>% 
  mutate(b_kr_Intercept = r_run__kr + b_kr_Intercept) %>% filter(term=="Intercept")

intercept_run$run <- as.character(intercept_run$run)

# Average effect: the grand mean
intercept_average <- spread_draws(asc_multi_bayes, b_kr_Intercept) %>% 
  mutate(run = "Grand mean")

# Combine average and study-specific effects data frames

intercept_all <- bind_rows(intercept_run, intercept_average) %>% 
  ungroup() %>% mutate(run=fct_relevel(run,"Grand mean","1","2","3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" ))

# Data frame of summary numbers
intercept_all_sum <- group_by(intercept_all, run) %>% 
  mean_qi(b_kr_Intercept)

# Draw forest/ridgeline plot for kr:

forest_kr <- intercept_all %>% ggplot(aes(b_kr_Intercept, run))+
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1,fill="Lightblue"
  ) +
  geom_pointintervalh(
    data = intercept_all_sum, size = 1
  ) +
  labs(x=expression(paste("k"[r])), y="run", subtitle="E")

multi_forest <- forest_c0+forest_nt+forest_kr

newavg <- data.frame(time = seq(from = 0, to = 7, by = 0.1))

fitavg2 <- cbind(newavg, fitted(asc_multi_bayes, newdata = newavg, re_formula = NA)[,-2])
#[,-2] means that the second column is omitted from fitted, containing the est.error and leaving estimate and Q2.5 and Q97.5, re_formula=NA means that group-level effects are omitted
names(fitavg2) <- c("time", "conc", "lower", "upper")

predavg2 <- cbind(newavg, predict(asc_multi_bayes, newdata=newavg, re_formula=NA)[,-2])
names(predavg2) <- c("time", "conc", "lower", "upper")

plot_me_fitavg <- ggplot(asc70, aes(x = time, y = conc)) +
  geom_point(size=2, shape=21, stroke=1, fill="red")+
  geom_line(data=fitavg2, aes(y=conc), size=1.5, colour="black") +
  geom_ribbon(data = fitavg2, aes(ymin = lower, ymax = upper), fill = "blue", alpha =0.5)+
  geom_ribbon(data = fitavg1, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.6)+
  geom_line(data=fitavg1, aes(y=conc), size=1.5, linetype=2, alpha=0.5) +
  coord_cartesian(ylim = c(0,7))+
  labs(x="time (h)", y="AA (mM)")

combi_plot3 <- plot_grid(ggmatrix_gtable(asc_corplot2),plot_me_fitavg,labels = c('A', 'B'),rel_widths = c(1.6, 1))
combi_plot4 <- plot_grid(combi_plot3, multi_forest, ncol=1)
combi_plot4

# Code for Table 3:

asc_var_cor <- summary(asc_multi_bayes)
asc_var_cor2 <- rbind(data.frame(asc_var_cor$fixed), data.frame(asc_var_cor$random$run), data.frame(asc_var_cor$spec_pars) )
rownames(asc_var_cor2) <- c("$c_0 \\text { (mM)}$", "$n_t (-)$", "$k_r (\\text {h}^{-1})$", "$\\sigma_u \\text {(mM)}$","$\\sigma_v (-)$", "$\\sigma_w (\\text {h}^{-1})$","$\\rho_{u,v}$", "$\\rho_{u,w}$","$\\rho_{v,w}$","$\\sigma_{e} \\text { (mM)}$")

colnames(asc_var_cor2) <- c("mean","SE", "lower bound", "upper bound")

asc_var_cor2[1:10,1:4] %>% 
  rownames_to_column(var = "parameter") %>% kbl(booktabs=T, escape=F, digits = c(3,3,3,3), caption = "Numerical summary of the posterior distribution resulting from bayesian multilevel modeling with varying initial concentration, reaction order, rate constant and correlation between them of the ascorbic acid data. SE = standard error, the lower and upper bound reflect the 95% credible interval") %>% kable_styling(position="center", full_width = F)

# First-order single level for comparison:
nlform<-bf(conc ~ c0*exp(-k*time), 
           c0~1, k~1,nl=TRUE)

nlprior<-c(prior(normal(6,1), nlpar = "c0"),
           prior(normal(0.3,0.3), nlpar="k"),
           prior(cauchy(0,10), class="sigma")
)           

asc_bayes_all_fo <- brm(formula=nlform, data=asc70, family = gaussian(), prior = nlprior, warmup=2000, iter=4000, chains=4, control = list(adapt_delta = 0.95), file=here("fits","asc_bayes_all_fo"))

# First-order multilevel for comparison:
nlform<-bf(conc ~ c0*exp(-kr*time), 
           c0 ~1+(1|ID|run),
           kr ~ 1+(1|ID|run), 
           nl=TRUE)

nlprior<-c(prior(normal(6.0,1.0), nlpar = "c0", lb=0),
           prior(normal(0.3,0.3), nlpar="kr", lb=0),
           prior(cauchy(0,10), class="sigma"),
           prior(cauchy(0,10), class="sd", nlpar="c0"),
           prior(cauchy(0,10), class="sd", nlpar="kr"),
           prior(lkj(2), class="cor")
)           

asc_multi_fo_bayes <- brm(formula=nlform, data=asc70, family=gaussian(),
                          prior=nlprior, iter=8000, warmup=4000, seed=1234, chains=4, control = list(adapt_delta=0.999, max_treedepth=15),
                          file=here("fits","asc_multi_fo_bayes"))

# Code for Figure 6:

pp1 <- pp_check(asc_bayes_all, ndraws = 100, type="ecdf_overlay")+ labs(x="AA (mM)",subtitle = "single level nth-order model")
pp2 <- pp_check(asc_multi_bayes, ndraws = 100,type="ecdf_overlay")+ labs(x="AA (mM)",subtitle="multilevel nth-order model")
pp3 <- pp_check(asc_bayes_all_fo, ndraws = 100,type="ecdf_overlay")+ labs(x="AA (mM)", subtitle="single level first-order model")
pp4 <- pp_check(asc_multi_fo_bayes, ndraws = 100,type="ecdf_overlay")+ labs(x="AA (mM)", subtitle="multilevel first-order model")

pp1 + pp2 + pp3 + pp4

# Code for model comparison:

asc_bayes_all <- add_criterion(asc_bayes_all, c("loo","waic"), file=here("fits","asc_bayes_all"))
asc_multi_bayes <- add_criterion(asc_multi_bayes, c("loo","waic"), file=here("fits","asc_multi_bayes"))
asc_multi_fo_bayes <- add_criterion(asc_multi_fo_bayes, c("loo","waic"), file=here("fits","asc_multi_fo_bayes"))
asc_bayes_all_fo <- add_criterion(asc_bayes_all_fo, c("loo","waic"), file=here("fits","asc_bayes_all_fo"))

model_comp_asc <- loo_compare(asc_bayes_all$criteria$loo,
                              asc_multi_bayes$criteria$loo,
                              asc_multi_fo_bayes$criteria$loo,
                              asc_bayes_all_fo$criteria$loo)

top_row4 <- model_comp_asc[,3:4] %>% data.frame() %>% rownames_to_column(var="model_name") %>% ggplot(aes(x=model_name, y=elpd_loo, ymin=elpd_loo-se_elpd_loo, ymax=elpd_loo+se_elpd_loo))+
  geom_pointrange(shape=21)+
  coord_flip()+
  labs(x=NULL, y=NULL)+
  theme(axis.ticks.y = element_blank())+
  labs(y=expression(paste("elpd"[loo])))+
  scale_x_discrete(labels=c("asc_multi_fo_bayes"="multilevel first-order", "asc_multi_bayes"="multilevel nth-order", "asc_bayes_all_fo"="single-level first-order", "asc_bayes_all"="single-level nth-order"))

# Code for Figure 7:

pareto1 <- function(){plot(asc_multi_bayes$criteria$loo, main="multilevel nth-order")}
pareto2 <- function(){plot(asc_bayes_all$criteria$loo, main="single level nth-order")}
pareto3 <- function(){plot(asc_multi_fo_bayes$criteria$loo, main="multilevel first-order")}
pareto4 <- function(){plot(asc_bayes_all_fo$criteria$loo, main="single level first-order")}
top_row4 <- plot_grid(top_row4, labels=c("A"))
bottom_row4 <- plot_grid(pareto1,labels = c("B"))
bottom_row5 <- plot_grid(pareto3,labels = c("C"))
top_row4/bottom_row4/bottom_row5

# Code for Figure 8:

plot_me_predavg <- ggplot(asc70, aes(x = time, y = conc)) +
  geom_point(size=2, shape=21, stroke=1, fill="red")+
  geom_line(data=predavg2, aes(y=conc), colour="black", size=1.5) +
  geom_ribbon(data = predavg2, aes(ymin = lower, ymax = upper), fill = "blue", alpha = .5)+
  geom_ribbon(data = predavg1, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha =0.5)+
  geom_line(data=fitavg1, aes(y=conc), size=1.5, linetype=2, alpha=0.9) +
  coord_cartesian(ylim = c(0,7.5))+
  labs(x="time (h)", y="AA (mM)")

# asc_ind combines the data with predictions using the random effects with re_formula = NULL
newvary <- expand.grid(time=seq(from = 0, to =7, by=0.5),run=1:15)

asc_ind <- cbind(newvary, predict(asc_multi_bayes, newdata=newvary, re_formula = NULL)[,-2])

names(asc_ind) <- c("time", "run", "conc", "lower", "upper")

asc_ind$run=as.factor(asc_ind$run)

ind_fit3 <- ggplot(asc70, aes(x=time, y=conc))+
  geom_point(size=2, shape=21, stroke=1, fill="red")+
  facet_wrap(~run, ncol=5)+
  geom_line(data = asc_ind, aes(y = conc), size = 1, colour="blue") +
  geom_ribbon(data = asc_ind, aes(ymin = lower, ymax = upper), fill = "blue", alpha = .3)+
  geom_line(data=fitavg2, aes(y=conc), colour="red")+
  labs(x="time (h)", y="AA (mM)", subtitle="A")

#Figure 8B:
  
newdata <- expand_grid(time=1)
pred_1_avg <- posterior_epred(asc_avg_norm_bayes, newdata) %>% as.data.frame()
pred_1_cp <- posterior_epred(asc_bayes_all, newdata) %>% as.data.frame()
pred_1_pp <- posterior_epred(asc_multi_bayes, newdata, re_formula = NA) %>% as.data.frame()

colors <- c(V1="red", pred_1_cp="blue", pred_1_cp="turquoise")

pred_plot <- pred_1_avg %>% ggplot(aes(x=V1*5.6))+
  geom_density(fill="red",alpha=0.5)+
  geom_density(data=pred_1_cp, aes(x=V1, fill="blue"), fill="blue",alpha=0.5)+
  geom_density(data=pred_1_pp, aes(x=V1), fill="turquoise", alpha=0.7)+
  labs(x="global prediction (mM), t=1 h 70 ?C", y="density")+
  scale_colour_manual(values=c("blue", "red")) +
  theme(legend.position=c(0.9, 0.9))+
  annotate("text", label="averaged-normalized: red", x=2.8, y=6, color="red", size=2)+
  annotate("text", label="completely pooled: blue", x=2.77, y=5, color="blue", size=2)+
  annotate("text", label="partially pooled: turquoise", x=2.79, y=4, color="turquoise",size=2)
pred_plot #(Figure 8B)


(plot_me_predavg+plot_spacer())/ind_fit3+plot_annotation(tag_levels = 'A')+plot_layout(heights=c(1,2)) # Figure 8A, C

# Code for Figure 9:

asc_ind5 <-  asc_ind %>% filter(run==5)

names(asc_ind5) <- c("time", "run", "conc", "lower", "upper")

newavg <- data.frame(time = seq(from = 0, to = 7, by = 0.1))

fitind5 <- cbind(newavg, fitted(asc_model5, newdata = newavg, re_formula = NA)[,-2])

asc_ind5_plot <- ggplot(asc_run5, aes(x=time, y=conc))+
  geom_point()+
  geom_line(data=asc_ind5, aes(y = conc))+
  geom_line(data=fitavg2, aes(y=conc), colour="red") +
  geom_line(data=fitind5, aes(y=Estimate), lty=2, colour="blue")+
  labs(subtitle = "A: run 5")

asc_ind9 <-  asc_ind %>% dplyr::filter(run==9)

names(asc_ind9) <- c("time", "run", "conc", "lower", "upper")

fitind9 <- cbind(newavg, fitted(asc_model9, newdata = newavg, re_formula = NA)[,-2])

asc_ind9_plot <- ggplot(asc_run9, aes(x=time, y=conc))+
  geom_point()+
  geom_line(data=asc_ind9, aes(y = conc))+
  geom_line(data=fitavg2, aes(y=conc), colour="red") +
  geom_line(data=fitind9, aes(y=Estimate), lty=2, colour="blue")+
  labs(subtitle = "B: run 9")

(shrink_plot <- asc_ind5_plot+asc_ind9_plot)

# Code used for supplementary information:

#Figure S2
(asc70plot_sep <-  asc70 %>% ggplot(aes(x=time, y=conc, group=run))+
    geom_point()+
    geom_line()+
    facet_wrap(~run, scale="free", ncol=4)+
    labs(x="time (h)", y="ascorbic acid (mM)")+
    theme(strip.background = element_rect(fill="lightblue", size=0.2, color="darkblue")))

# Code for Figure S3:

#Figure S3A
asc_cor_fo1 <- dplyr::select(post_asc1_fo,b_c0_Intercept:sigma)
asc_cor_fo1 <- setNames(asc_cor_fo1, c(expression(paste("c"[0])), expression(paste("k"[r])), expression(sigma[e])))

asc_cor_fo1_plot <-asc_cor_fo1  %>%  ggpairs(diag=list(continuous="densityDiag"),
                                             mapping=aes(fill="red"),
                                             upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
                                             labeller=label_parsed)+ 
  theme(strip.text.x = element_text(size = 10, color = "red"),
        strip.text.y = element_text(size = 10, color = "red"))+
  theme(axis.text.x = element_text(angle=70, hjust=1))

#Figure S3B
fitplot_fo <- ggplot(asc70, aes(x=time, y=conc))+
  geom_point()+
  geom_ribbon(data=fit_unpooled_all_fo,aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),fill = "lightblue", alpha=0.5) +
  geom_line(data=fit_unpooled_all_fo, aes(x=time, y=Estimate))+
  facet_wrap(~run, ncol=8)+
  labs(x="time (h)", y="AA, mmol/L")

#Figure S3C

resid <- resid(asc_model1_fo)[, "Estimate"]
resid_run1 <- data.frame(resid) %>% bind_cols(time=asc_run1$time)
resid_run1$run <- as.factor(rep(1, length(nrow(time.seq))))

resid <-  resid(asc_model2_fo)[, "Estimate"]
resid_run2 <- data.frame(resid) %>% bind_cols(time=asc_run2$time)
resid_run2$run <- as.factor(rep(2, length(nrow(time.seq))))

resid = resid(asc_model3_fo)[, "Estimate"]
resid_run3 <- data.frame(resid) %>% bind_cols(time=asc_run3$time)
resid_run3$run <- as.factor(rep(3, length(nrow(time.seq))))

resid = resid(asc_model4_fo)[, "Estimate"]
resid_run4 <- data.frame(resid) %>% bind_cols(time=asc_run4$time)
resid_run4$run <- as.factor(rep(4, length(nrow(time.seq))))

resid = resid(asc_model5_fo)[, "Estimate"]
resid_run5 <- data.frame(resid) %>% bind_cols(time=asc_run5$time)
resid_run5$run <- as.factor(rep(5, length(nrow(time.seq))))

resid <-  resid(asc_model6_fo)[, "Estimate"]
resid_run6 <- data.frame(resid) %>% bind_cols(time=asc_run6$time)
resid_run6$run <- as.factor(rep(6, length(nrow(time.seq))))

resid <-  resid(asc_model7_fo)[, "Estimate"]
resid_run7 <- data.frame(resid) %>% bind_cols(time=asc_run7$time)
resid_run7$run <- as.factor(rep(7, length(nrow(time.seq))))

resid <- resid(asc_model8_fo)[, "Estimate"]
resid_run8 <- data.frame(resid) %>% bind_cols(time=asc_run8$time)
resid_run8$run <- as.factor(rep(8, length(nrow(time.seq))))

resid <-  resid(asc_model9_fo)[, "Estimate"]
resid_run9 <- data.frame(resid) %>% bind_cols(time=asc_run9$time)
resid_run9$run <- as.factor(rep(9, length(nrow(time.seq))))

resid <-  resid(asc_model10_fo)[, "Estimate"]
resid_run10 <- data.frame(resid) %>% bind_cols(time=asc_run10$time)
resid_run10$run <- as.factor(rep(10, length(nrow(time.seq))))

resid <-  resid(asc_model11_fo)[, "Estimate"]
resid_run11 <- data.frame(resid) %>% bind_cols(time=asc_run11$time)
resid_run11$run <- as.factor(rep(11, length(nrow(time.seq))))

resid <-  resid(asc_model12_fo)[, "Estimate"]
resid_run12 <- data.frame(resid) %>% bind_cols(time=asc_run12$time)
resid_run12$run <- as.factor(rep(12, length(nrow(time.seq))))

resid <-  resid(asc_model13_fo)[, "Estimate"]
resid_run13 <- data.frame(resid) %>% bind_cols(time=asc_run13$time)
resid_run13$run <- as.factor(rep(13, length(nrow(time.seq))))

resid <-  resid(asc_model14_fo)[, "Estimate"]
resid_run14 <- data.frame(resid) %>% bind_cols(time=asc_run14$time)
resid_run14$run <- as.factor(rep(14, length(nrow(time.seq))))

resid <-  resid(asc_model15_fo)[, "Estimate"]
resid_run15<- data.frame(resid) %>% bind_cols(time=asc_run15$time)
resid_run15$run <- as.factor(rep(15, length(nrow(time.seq))))

resid_all_fo <- rbind(resid_run1,resid_run2,resid_run3,resid_run4,resid_run5,resid_run6,resid_run7,resid_run8,resid_run9,resid_run10,resid_run11,resid_run12,resid_run13,resid_run14,resid_run15)

residplot_fo <-  ggplot(data = resid_all_fo, aes(x=time, y=resid)) + 
  geom_point()+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  xlim(0,7)+
  facet_wrap(~run, ncol = 8)+
  labs(x="time (h)", y="residual")

#Figure S3
top_row1 <- plot_grid(ggmatrix_gtable(asc_cor_fo1_plot),labels = c('A'))
middle_row1 <- plot_grid(fitplot_fo, labels=c('B'))
bottom_row1 <- plot_grid(residplot_fo,labels = c('C'))

combi_plot1 <- plot_grid(top_row1, middle_row1,bottom_row1, ncol=1, align="v")
combi_plot1

#Figure S4
asc_cor_nth <- dplyr::select(post_asc1,b_c0_Intercept:sigma)
asc_cor_nth <- setNames(asc_cor_nth, c(expression(paste("c"[0])), expression(paste("n"[t])),expression(paste("k"[r])), expression(sigma[e])))

(asc_corplot_nth <-asc_cor_nth  %>%  ggpairs(diag=list(continuous="densityDiag"),
                                             mapping=aes(fill="red"),
                                             upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
                                             labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red"))+
    theme(axis.text.x = element_text(angle=70, hjust=1))
)

# brms regression of single-level first-order model:
nlform<-bf(conc ~ c0*exp(-k*time), 
           c0~1, k~1,nl=TRUE)

nlprior<-c(prior(normal(6,1), nlpar = "c0"),
           prior(normal(0.3,0.3), nlpar="k"),
           prior(cauchy(0,10), class="sigma")
)           

asc_bayes_all_fo <- brm(formula=nlform, data=asc70, family = gaussian(), prior = nlprior, warmup=2000, iter=4000, chains=4, control = list(adapt_delta = 0.95), file=here("fits","asc_bayes_all_fo"))

asc_bayes_all_fo_post <- as_draws_df(asc_bayes_all_fo)

#Figure S4
asc_cor_fo <- dplyr::select(asc_bayes_all_fo_post,b_c0_Intercept:sigma)
asc_cor_fo <- setNames(asc_cor_fo, c(expression(paste("c"[0])), expression(paste("k"[r])), expression(sigma[e])))

asc_corplot_fo <-asc_cor_fo  %>%     ggpairs(diag=list(continuous="densityDiag"),
                                             mapping=aes(fill="red"),
                                             upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
                                             labeller=label_parsed)+ 
  theme(strip.text.x = element_text(size = 10, color = "red"),
        strip.text.y = element_text(size = 10, color = "red"))+
  theme(axis.text.x = element_text(angle=70, hjust=1))

# Figure S5:

newavg <- data.frame(time = seq(from = 0, to = 7, by = 0.1))

fitavg1_fo <- cbind(newavg, fitted(asc_bayes_all_fo, newdata = newavg, re_formula = NA)[,-2])
#[,-2] means that the second column is omitted from fitted, containing the est.error and leaving estimate and Q2.5 and Q97.5, re_formula=NA means that group-level effects are omitted
names(fitavg1_fo) <- c("time", "conc", "lower", "upper")
#same for prediction intervals
predavg1_fo <- cbind(newavg, predict(asc_bayes_all_fo, newdata=newavg, re_formula=NA)[,-2])
names(predavg1_fo) <- c("time", "conc", "lower", "upper")

plot_fo_asc <- ggplot(asc70, aes(x = time, y = conc)) +
  geom_point(size=2, shape=21, stroke=1, fill="red") +
  geom_line(data=fitavg1_fo, aes(y=conc)) +
  geom_ribbon(data = fitavg1_fo, aes(ymin = lower, ymax = upper), fill = "blue", alpha = .3)+
  geom_ribbon(data = predavg1_fo, aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = .6)+
  coord_cartesian(ylim = c(0,7.5))+
  labs(x="time (h)", y= "[AA] (mM)")

#residuals, QQ and lag plots
resid_asc1_fo = resid(asc_bayes_all_fo)[, "Estimate"]

residplot_asc1_fo <-  ggplot() + geom_point(data = NULL, size = 2, aes(x = asc70$time, y = resid_asc1_fo))+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  labs(x = "time (h)", y = "residual")

# calculating the QQ plot for residuals:
res_asc1_fo <- as.data.frame(residuals(asc_bayes_all_fo, summary=T))

qqplot_asc1_fo <-  ggplot(res_asc1_fo, aes(sample = Estimate)) +
  stat_qq() +
  stat_qq_line() 

#lagplot
lagx <- resid_asc1_fo[-1]
lagy <- resid_asc1_fo[-length((resid_asc1_fo))]
lag <- data.frame(lagx,lagy)
lagplot_asc1_fo <- ggplot(lag, aes(x=lagx, y=lagy))+
  geom_point(size=1)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(x="residual(i-1)" , y="residual(i)")

top_row2 <- plot_grid(ggmatrix_gtable(asc_corplot_fo),plot_fo_asc,labels = c('A', 'B'))

bottom_row2 <- plot_grid(residplot_asc1_fo,qqplot_asc1_fo,lagplot_asc1_fo, ncol = 3,labels = c('C', 'D','E'))

combi_plot2 <- plot_grid(top_row2, bottom_row2, nrow=2, align="v")
combi_plot2

# brms multilevel regression of first-order model:

nlform<-bf(conc ~ c0*exp(-kr*time), 
           c0 ~1+(1|ID|run),
           kr ~ 1+(1|ID|run), 
           nl=TRUE)

nlprior<-c(prior(normal(6.0,1.0), nlpar = "c0", lb=0),
           prior(normal(0.5,0.3), nlpar="kr", lb=0),
           prior(cauchy(0,10), class="sigma"),
           prior(cauchy(0,10), class="sd", nlpar="c0"),
           prior(cauchy(0,10), class="sd", nlpar="kr"),
           prior(lkj(2), class="cor")
)           

asc_multi_fo_bayes <- brm(formula=nlform, data=asc70, family=gaussian(),
                          prior=nlprior, iter=8000, warmup=4000, seed=1234, chains=4, control = list(adapt_delta=0.999, max_treedepth=15),
                          file=here("fits","asc_multi_fo_bayes"))

asc_multi_fo_bayes_post <- as_draws_df(asc_multi_fo_bayes)

#Figure S6
asc_cor_multi_fo <- dplyr::select(asc_multi_fo_bayes_post,b_c0_Intercept:sigma)

asc_cor_multi_fo <- setNames(asc_cor_multi_fo, c(expression(paste("c"[0])), expression(paste("k"[r])), expression(sigma[c[0]]),expression(sigma[k[r]]),expression(rho[c[0]-k[r]]), expression(sigma[e])))

asc_corplot_multi_fo <-asc_cor_multi_fo  %>% ggpairs(diag=list(continuous="densityDiag"),upper = list(continuous = wrap("cor", size = 2, stars=FALSE, title="corr. coef")), mapping=aes(fill="red"), labeller=label_parsed)+ 
  theme(strip.text.x = element_text(size = 8, color = "red"),
        strip.text.y = element_text(size = 8, color = "red"))+
  theme(axis.text.x = element_text(size=3,angle=70, hjust=1))+
  theme(axis.text.y = element_text(size=3))  

intercept_run <- asc_multi_fo_bayes %>% spread_draws(b_c0_Intercept, r_run__c0[run,term]) %>% 
  mutate(b_c0_Intercept = r_run__c0 + b_c0_Intercept) %>% filter(term=="Intercept")

intercept_run$run <- as.character(intercept_run$run)

# Average effect: the grand mean
intercept_average <- spread_draws(asc_multi_fo_bayes, b_c0_Intercept) %>% 
  mutate(run = "Grand mean")

# Combine average and study-specific effects data frames

intercept_all <- bind_rows(intercept_run, intercept_average) %>% 
  ungroup() %>% mutate(run=fct_relevel(run,"Grand mean","1","2","3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" )) 
# Data frame of summary numbers
intercept_all_sum <- group_by(intercept_all, run) %>% 
  mean_qi(b_c0_Intercept)

forest_c0_fo <- intercept_all %>% ggplot(aes(b_c0_Intercept, run))+
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1,fill="Lightblue"
  ) +
  geom_pointintervalh(
    data = intercept_all_sum, size = 1
  ) +
  labs(x=expression(paste("c"[0])), y="run")

#> Warning: unnest() has a new interface. See ?unnest for details.
#> Try `cols = c(.lower, .upper)`, with `mutate()` needed

intercept_run <- asc_multi_fo_bayes %>% spread_draws(b_kr_Intercept, r_run__kr[run,term]) %>% 
  mutate(b_kr_Intercept = r_run__kr + b_kr_Intercept) %>% filter(term=="Intercept")

intercept_run$run <- as.character(intercept_run$run)

# Average effect: the grand mean
intercept_average <- spread_draws(asc_multi_fo_bayes, b_kr_Intercept) %>% 
  mutate(run = "Grand mean")

# Combine average and study-specific effects data frames

intercept_all <- bind_rows(intercept_run, intercept_average) %>% 
  ungroup() %>% mutate(run=fct_relevel(run,"Grand mean","1","2","3","4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15" )) 

# Data frame of summary numbers
intercept_all_sum <- group_by(intercept_all, run) %>% 
  mean_qi(b_kr_Intercept)

forest_kr_fo <- intercept_all %>% ggplot(aes(b_kr_Intercept, run))+
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1,fill="Lightblue"
  ) +
  geom_pointintervalh(
    data = intercept_all_sum, size = 1
  ) +
  labs(x=expression(paste("k"[r])), y="run")

multi_fo_forest <- forest_c0_fo + forest_kr_fo

# asc_ind combines the data with predictions using the random effects with re_formula = NULL, fo=first-order
newvary <- expand.grid(time=seq(from = 0, to =7, by=0.5),run=1:15)

asc_ind_fo <- cbind(newvary, predict(asc_multi_fo_bayes, newdata=newvary, re_formula = NULL)[,-2])

names(asc_ind_fo) <- c("time", "run", "conc", "lower", "upper")

asc_ind_fo$run=as.factor(asc_ind_fo$run)

ind_fit_fo <- ggplot(asc70, aes(x=time, y=conc))+
  geom_point(size=2, shape=21, stroke=1, fill="red")+
  facet_wrap(~run, ncol=5)+
  geom_line(data = asc_ind_fo, aes(y = conc), size = 1, colour="blue") +
  geom_ribbon(data = asc_ind_fo, aes(ymin = lower, ymax = upper), fill = "blue", alpha = .3)+
  labs(x="time (h)", y="AA (mM)")

top_row_fo <- plot_grid(ggmatrix_gtable(asc_corplot_multi_fo),labels = c('A'))
middle_row_fo <- plot_grid(ind_fit_fo,labels = c('B'))
bottom_row_fo <- plot_grid(multi_fo_forest, labels = c('C'))
combi_plot_fo <- plot_grid(top_row_fo, middle_row_fo,bottom_row_fo, ncol=1)
combi_plot_fo



