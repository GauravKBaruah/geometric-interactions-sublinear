, 


########### FIGURE 4 #######################



# list of parameters
S  <- 2 #no of species
m0 <- rep(0.25, S) #mortality rate constant
b  <- c(1.0, 1.0) # growth rate constant
a11 <- 1 #intra comp
a22 <- 1 #intra competition
alpha_21_fixed <- 0.4
eta0  <- 0.5
c_F   <- c(1, 1)

# rho - niche overlap
rho_seq <- seq(0.05, 0.98, length.out = 500)

# gLV exclusion boundary in niche-difference
lv_nd_boundary <- 1 - sqrt(alpha_21_fixed) 

# eta_12 is the interference exerted by species 2
#on species 1. k_asym = eta_12 / eta_21

make_eta_asym <- function(k_asym) {
  
  eta_11 <- eta0
  eta_12 <- eta0 * k_asym
  eta_21 <- eta0
  eta_22 <- eta0
  
  matrix(c(
      eta_11, eta_12,
      eta_21, eta_22),  nrow = 2, byrow = TRUE)
}


# eta_22 is self-interference within species 2.
# k_self = eta_22 / eta_11

make_eta_self <- function(k_self) {
  
  eta_11 <- eta0
  eta_12 <- eta0
  eta_21 <- eta0
  eta_22 <- eta0 * k_self
  
  matrix( c(
      eta_11, eta_12,
      eta_21, eta_22),nrow = 2,  byrow = TRUE)
}

# eta_21 is the interference exerted by species 1
# on species 2.
#
# k_cross = eta_21 / eta_11
make_eta_cross <- function(k_cross) {
  
  eta_11 <- eta0
  eta_12 <- eta0
  eta_21 <- eta0 * k_cross
  eta_22 <- eta0
  
  matrix(c(eta_11, eta_12,
      eta_21, eta_22),nrow = 2,byrow = TRUE)
}

run_interference_sweep <- function(knob_seq,rho_seq,gamma_seq,make_eta) {
  
  grid <- expand.grid(rho = rho_seq,knob = knob_seq,g = gamma_seq)

  grid$niche_diff <- 1 - grid$rho
  
  for (i in seq_len(nrow(grid))) {
    
    rho_i<- grid$rho[i]
    knob_i<- grid$knob[i]
    gamma_i<- grid$g[i]
      alpha_12 <-(rho_i^2 *a11*a22)/alpha_21_fixed
    
    a_mat <- matrix(c(a11, alpha_12, alpha_21_fixed, a22),nrow = 2,byrow = TRUE)
    eta_mat <- make_eta(knob_i)
    output <- run_to_eq(eta_mat = eta_mat,a_mat = a_mat,b_vec = b,m0_vec = m0,c_vec = c_F,gamma = gamma_i)
    
    grid$N1[i] <- output[1]
    grid$N2[i] <- output[2]
    grid$outcome[i] <- as.character(classify(grid$N1[i],grid$N2[i]))
  }
  return(grid)
}

k_asym_seq<- seq(0.01,5.0,length.out = 50)
k_self_seq  <- seq(0.01,5.0,length.out = 50)
k_cross_seq<- seq(0.01,5.0,length.out = 50)
gamma <- c(0.5,1)

# Parameter sequences
k_asym_seq <- seq(0.01,5,length.out = 50)
k_self_seq <- seq(0.01,5,length.out = 50)
k_cross_seq <- seq(0.01,5,length.out = 50)
gamma_seq <- c(0.5, 1)

g1 <- run_interference_sweep(knob_seq = k_asym_seq, rho_seq = rho_seq,  
                             gamma_seq = gamma_seq,make_eta = make_eta_asym)

g2 <- run_interference_sweep(  knob_seq = k_self_seq,
  rho_seq = rho_seq,gamma_seq = gamma_seq,make_eta = make_eta_self)

g3 <- run_interference_sweep(knob_seq = k_cross_seq,rho_seq = rho_seq,
                             gamma_seq = gamma_seq,make_eta = make_eta_cross)

#make_plot <- function(g, ylab_expr, ttl) {
  a4<-ggplot(g1, aes(x = niche_diff, y = knob, fill = outcome)) +
    geom_raster(interpolate = FALSE) +
      geom_vline(xintercept = lv_nd_boundary,
               linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dotted",
               color = "grey85", linewidth = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = outcome_colors) +
    labs(x = expression("Niche difference  " * (1 - rho)),
         y =  expression(k[asym] * " = " * eta[12] / eta[21]), fill = "Outcome",
         title = "") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),legend.position = "NA") + facet_wrap(.~g)

a5<-  ggplot(g2, aes(x = niche_diff, y = knob, fill = outcome)) +
    geom_raster(interpolate = FALSE) +
    geom_vline(xintercept = lv_nd_boundary,
               linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dotted",
               color = "grey85", linewidth = 0.5) +
    scale_fill_manual(values = outcome_colors) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = expression("Niche difference  " * (1 - rho)),
         y =  expression(k[self] * " = " * eta[22] / eta[11]), fill = "Outcome",
         title = "") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "NA") + facet_wrap(.~g)
  
a6<-  ggplot(g3, aes(x = niche_diff, y = knob, fill = outcome)) +
    geom_raster(interpolate = FALSE) +
    geom_vline(xintercept = lv_nd_boundary,
               linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dotted",
               color = "grey85", linewidth = 0.5) +
    scale_fill_manual(values = outcome_colors) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = expression("Niche difference  " * (1 - rho)),
         y =  expression(k[cross] * " = " * eta[21] / eta[11]), fill = "Outcome",
         title = "F)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),legend.position = "NA") + facet_wrap(.~g)
  
  
  

final<-ggpubr::ggarrange(a3,a2,a1,a4,a5, nrow=2, ncol=3,
                         labels=c("A","B","C","D", "E"))
final

