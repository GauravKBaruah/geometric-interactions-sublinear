#R script for figure 5

library(deSolve)
library(dplyr)
library(ggplot2) 
library(tidyr)

#niche overlap and fitness ratio 
#mu1- trait value of species 1
#mu2 - trait value of species 2
#s2_sq - trait variance species 2
#s1_sq - trait variance species 1
#w - width of competition kernel
#K - equilibrium density
#theta - optimum
#returns niche overlap and fitness ratio
mct_metrics <- function(mu1, mu2, s1_sq, s2_sq, w, K, theta) {
  rho <- exp(-(mu1 - mu2)^2 / (w^2 + 2*s1_sq + 2*s2_sq))*((w^2 + 4*s1_sq)*(w^2 + 4*s2_sq)/(w^2 + 2*s1_sq + 2*s2_sq)^2)^(1/4)
  num <- K[1]*theta^2 - mu1^2 - s1_sq
  den <- K[2]*theta^2 - mu2^2 - s2_sq
  k_ratio <- (num / den) * ((w^2 + 4*s1_sq) / (w^2 + 4*s2_sq))^(1/4)
  list(rho = rho, niche_diff = 1 - rho, k_ratio = k_ratio)
}

#simulate for two species the dynamical equations:
#sigma_val- vector of trait variances
#c_val - vector of c values
#gamma - 0.5 to 1
#returns- niche difference and fitness ratio
simulate_one <- function(sigma_val, c_val, gamma, label, tmax = 100000) {
  S      <- 2
  m_init <- c(0.21, 0.20)
  params <- list(S= S,m = m_init,w = 0.1,sigma = rep(sigma_val, S),h2= c(0.1, 0.15),theta = 0.75, d  = rep(0.25, S),c = rep(c_val, S),gamma=gamma,eta =  matrix(c(1, .01,.01, 1), 
                  nrow = 2, byrow = FALSE),points = 10)
  ic <- c(rep(1, S), m_init)
  out <- ode(func = eqs, y = ic, parms = params,times = seq(0, tmax, by = tmax/200), method = "lsoda", rtol = 1e-8, atol = 1e-10)
  out_df <- as.data.frame(out)
  names(out_df) <- c("time", "N1", "N2", "mu1", "mu2")
   K <- 1 - params$d                  # K_i = b_max - m0 = 1 - d_i
  mct_list <- mapply(mct_metrics,mu1 = out_df$mu1, mu2 = out_df$mu2,MoreArgs = list(s1_sq = sigma_val, s2_sq = sigma_val,
                    w = params$w, K = K, theta = params$theta),SIMPLIFY = FALSE)
  
  out_df$rho<- sapply(mct_list, `[[`, "rho")
  out_df$niche_diff<- sapply(mct_list, `[[`, "niche_diff")
  out_df$k_ratio<- sapply(mct_list, `[[`, "k_ratio")
  out_df$scenario <- label
  out_df$sigma_val<- sigma_val
  out_df$c_val<- c_val
  out_df$gamma  <- gamma
  out_df
}


#different scenario lists tested
scenarios <- list( list(sigma = 0.01, c = 0, gamma=0, label = "Low sigma, LV"),
  list(sigma = 0.01, c = 1, gamma=0.5, label = "Low sigma, sublinear, gamma =0.5 "),
  list(sigma = 0.01, c = 1, gamma=1, label = "Low sigma, sublinear, gamma =1"),
  list(sigma = 0.10, c = 0, gamma=0, label = "High sigma, LV"),
  list(sigma = 0.10, c = 1, gamma =0.5, label = "High sigma, sublinear, gamma =0.5"),
  list(sigma = 0.10, c = 1, gamma =1, label = "High sigma, sublinear, gamma =1"))


#trajectories 
trajs <- bind_rows(lapply(scenarios, function(s) simulate_one(s$sigma, s$c,s$gamma, s$label)))

x_lims <- c(0, 1)
y_lims <- c(0.3, 5)
rho_grid <- seq(0.001, 1, length.out = 1000)
wedge <- data.frame(rho = rho_grid, upper = 1/rho_grid,lower = rho_grid)

trajs <- trajs %>% mutate(sigma_lab = ifelse(sigma_val < 0.05, "Low sigma", "High sigma"),
         gamma_lab = factor(gamma,levels = c(0, 0.5, 1),labels = c("LV (c = 0)","gamma == 0.5", "gamma == 1")))

start_end <- trajs %>%group_by(scenario) %>%filter(time == min(time) | time == max(time)) %>%mutate(point_type = ifelse(time == min(time), "start", "end")) %>%
  ungroup()

mct1 <- ggplot() +geom_ribbon(data = wedge,aes(x = rho, ymin = lower, ymax = upper),fill = "grey88", alpha = 0.5
  ) +geom_line(data = wedge,aes(x = rho, y = upper),linetype = "dashed", color = "grey30") +
  geom_line( data = wedge,aes(x = rho, y = lower),linetype = "dashed", color = "grey30") +
  geom_path(data = trajs, aes(x = rho, y = k_ratio, color = sigma_lab),linewidth = 0.7, alpha = 0.8, linetype="dashed")  +
  geom_point(data = start_end %>% filter(point_type == "start"),aes(x = rho, y = k_ratio, color = sigma_lab), shape = 1, size = 3) +
  geom_point(data = start_end %>% filter(point_type == "end"),aes(x = rho, y = k_ratio, color = sigma_lab), shape = 16, size = 3)  +
  scale_y_log10() +
  scale_x_continuous() +
  coord_cartesian(xlim = x_lims, ylim = y_lims, expand = FALSE) +
  scale_color_manual(values = c("Low sigma"  = "#1b9e77","High sigma" = "#d95f02")) +
  facet_wrap(~ gamma_lab, nrow = 1, labeller = label_parsed) +
  labs( x = expression("Niche overlap  " * (rho)),y = expression("Fitness ratio  " * kappa[1] / kappa[2]),
    color = expression(sigma),title = "") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = expression("Niche overlap  " * ( rho)),
    y = expression("Fitness ratio  " * kappa[1] / kappa[2]),
    color = "Scenario",linetype = NULL, title = "") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

mct1

library(patchwork)

p_full <-mct1

# Zoomed inset around the apex
p_zoom <- p_full +
  coord_cartesian(xlim = c(0.85, 1.0), ylim = c(0.85, 1.05)) +labs(title = NULL, subtitle = NULL,
       x = NULL, y = NULL) +theme(legend.position = "none",
        plot.background = element_rect(fill = "white",color = "grey50",linewidth = 0.5), plot.margin = margin(4, 4, 4, 4))
p_zoom
p_full_annotated <- p_full +annotate("rect", xmin = 0.85, xmax = 1, ymin = 0.5, ymax = 2,fill = NA, color = "grey30", linewidth = 0.5,linetype = "dotted")

p_final <- p_full_annotated + inset_element(p_zoom,left   = 0.55, right = 0.99,bottom = 0.5, top   = 1.0)

print(p_final)




