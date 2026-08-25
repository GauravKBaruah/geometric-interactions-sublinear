
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(nleqslv) 
source("functions_sublinear.R")

alpha_from_traits <- function(mu, varA, w) {
  dmA <- outer(mu, mu, "-")
  svA <- outer(varA, varA, "+")
  exp(-dmA^2 / (2*svA + w^2)) * w / sqrt(2*svA + w^2)
}

w     <- 0.1
theta <- 0.75
tmax  <- 10000
kappa_seq <- c(0.2, 0.4, 0.8, 1, 2, 4, 8, 10)
fact <- expand.grid(S= c(30, 40, 50), model= c("gamma=1", "gamma=0", "LV"),variation = c("high","low"), kappa= kappa_seq,rep= 1:20)
fact <- fact[!(fact$model == "LV" & fact$kappa != fact$kappa[1]), ]
rownames(fact) <- NULL

fact$eta_scale <- fact$eta_ii <- fact$richness <- fact$frac_richness <-
  fact$inv_simpson <- fact$traitCV <- NA_real_

for (r in 1:nrow(fact)) {
  
  S <- fact$S[r]
  kappa <- fact$kappa[r]
  
  #same community for a given across model and kappa
  set.seed(1e5 * as.integer(fact$variation[r]) + 1e3 * S + fact$rep[r])
  
  d <- runif(S, 0.2, 0.25)
  if (fact$variation[r] == "high") {
    sigma <-runif(S, 0.01, 0.1)
  } else {
    sigma <-runif(S, 0.001, 0.005)
  }
  h2 <-runif(S, 0.1, 0.15)
  mu_max <- sqrt(pmax(theta^2 * (1 - d) - sigma, 0))
  m <- runif(S, -0.95 * mu_max, 0.95 * mu_max)
  
  if (fact$model[r] == "LV") {
    c <- 0 }else if (fact$model[r] == "gamma=1"){
    c <- 1; gamma <- 0
    gamma <- 1
    }else{
    c <- 1;gamma <- 0.5
  } 
  
 
  a_mat <- alpha_from_traits(m, sigma, w)
  eta_scale <- mean(diag(a_mat))
  eta <- eta_scale * matrix(runif(S*S, 0, 2), S, S)  # mean= eta_scale
  diag(eta) <- kappa * eta_scale                            #kappa = eta_ii/mean(eta_ij)
  
  init<- rep(1, S)
  inmu<- m
  ic<-c(init, m)
  
  params2 <- list(S = S, m = m, w = w, sigma = sigma, h2 = h2,
                  theta = theta, d = d, ic = ic, c = c, eta = eta,gamma=gamma)
  
  OUT <- ode(func = eqs, y = ic, parms = params2,
             times = seq(0, tmax, by = 100)) %>%
    organize_results_ode(pars = params2)
  #print(OUT %>% plot_all(moment=tmax))
  fin<- OUT %>% filter(time == tmax, n >= 1e-5)
  Ns<- fin$n
  mus<- fin$m
  p <- Ns / sum(Ns)
  
  fact$eta_scale[r] <- eta_scale
  fact$eta_ii[r]<- kappa * eta_scale
  fact$richness[r]<- length(Ns)
  fact$frac_richness[r] <- length(Ns) / S
  fact$inv_simpson[r]<- 1 / sum(p^2)
  fact$traitCV[r] <- traitCV(mus = mus)
  
  print(r)
}


library(ggthemes)
library(cowplot)
library(plyr)
kappa_levels <- sort(unique(fact$kappa[fact$model == "gamma=0"]))


lv_expanded <- fact %>%
  filter(model == "LV") %>%
  select(-kappa) %>%
  tidyr::crossing(kappa = kappa_levels)

plotdat <- bind_rows(fact %>% filter(model != "Sublinear"), lv_expanded) %>%
  mutate(kappa_lab = factor(round(kappa, 3), levels = round(kappa_levels, 3)),
         variation = factor(variation, levels = c("low", "high"),
                            labels = c("Low trait variation",
                                       "High trait variation")))

plotdat$model <- factor(plotdat$model,
                        levels = c("LV", "gamma=0", "gamma=1"),
                        labels = c("LV", "gamma == 0", "gamma == 1"))

plotdat$model <- revalue(plotdat$model, c("gamma == 0"="gamma == 0.5", "gamma == 1" = "gamma == 1",
                                          "LV"="LV"))

p1 <- plotdat %>%
  ggplot(aes(x = kappa_lab, y = inv_simpson, colour = factor(model))) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.25),
              alpha = 0.4, size = 2, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.8,
               width = 0.6) +
  geom_vline(xintercept = which(round(kappa_levels, 3) ==
                                  round(kappa_levels[which.min(abs(kappa_levels - 1))], 3)),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(
    values = c("LV"           = "#1b9e77",
               "gamma == 0.5" = "#7570b3",
               "gamma == 1"   = "#d95f02"),
    labels = parse(text = levels(plotdat$model))
  )+
  labs(x = expression(psi == eta[ii] / bar(eta)[ij] *
                        " scaled intra-inteference"),
       y = "Inverse simpson's diversity",
       colour = "Model") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  facet_grid(S ~ variation, scales = "free_y")

p1


plotdat %>%
  ggplot(aes(x = kappa_lab, y = richness, colour = factor(model))) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.25),
              alpha = 0.4, size = 2, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.8,
               width = 0.6) +
  geom_vline(xintercept = which(round(kappa_levels, 3) ==
                                  round(kappa_levels[which.min(abs(kappa_levels - 1))], 3)),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(
    values = c("LV"           = "#1b9e77",
               "gamma == 0.5" = "#7570b3",
               "gamma == 1"   = "#d95f02"),
    labels = parse(text = levels(plotdat$model))
  )+
  labs(x = expression(psi == eta[ii] / bar(eta)[ij] *
                        " scaled intra-inteference"),
       y = "Richness",
       colour = "Model") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  facet_grid(S ~ variation, scales = "free_y")

plotdat %>%
  ggplot(aes(x = kappa_lab, y = traitCV, colour = factor(model))) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.25),
              alpha = 0.4, size = 2, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.8,
               width = 0.6) +
  geom_vline(xintercept = which(round(kappa_levels, 3) ==
                                  round(kappa_levels[which.min(abs(kappa_levels - 1))], 3)),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(
    values = c("LV"           = "#1b9e77",
               "gamma == 0.5" = "#7570b3",
               "gamma == 1"   = "#d95f02"),
    labels = parse(text = levels(plotdat$model))
  )+
  labs(x = expression(psi == eta[ii] / bar(eta)[ij] *
                        " scaled intra-inteference"),
       y = "Trait clustering",
       colour = "Model") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  facet_grid(S ~ variation, scales = "free_y")
