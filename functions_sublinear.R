
library(cowplot)
library(statmod)

erf <- function(x) 2*pnorm(x*sqrt(2))-1

cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)
traitCV <- function(mus){
  
  v <- sort(mus[is.na(mus)==FALSE])
  d <- diff(v)
  cvs<-sd(d,na.rm=T)/mean(d,na.rm=T)
  cvs[is.na(cvs)]<-0
  
  return(cvs)
}


classify <- function(N1, N2, thr = 1e-6) {
  dplyr::case_when(
    N1 > thr & N2 > thr   ~ "Coexistence",
    N1 > thr & N2 <= thr  ~ "Sp 1 survives",
    N1 <= thr & N2 > thr  ~ "Sp 2 survives",
    TRUE                  ~ "Both extinct"
  )
}

outcome_colors <- c("Coexistence"  = "#2ca25f",
                    "Sp 1 survives"    = "grey85",
                    "Sp 2 survives"    = "grey85",
                    "Both extinct" = "grey85")

#subliner dynamics
# sublinear dynamics with Q_i * Q_j competition structure
sublinear_eqs <- function(t, state, pars) {
  N <- pmax(state, 0)
  N[N < 1e-10] <- 0
  crowding_effect <- as.vector(pars$eta %*% N)
  Q <- 1/(1+pars$c*(crowding_effect)^pars$gamma)
  competition <- as.vector(pars$a %*% (Q*N))
  dndt <- N*(pars$b*Q - pars$m_0 - Q*competition)
  list(c(dndt))
}


eqs <- function(time, state, pars) {
  S <- length(pars$sigma) ## number of species
  n <- state[1:S] ## species densities
  m <- state[(S+1):(2*S)] ## species trait means
  varA <- pars$sigma ## species' trait standard deviations
  w<-pars$w
  theta<-pars$theta
  d<- pars$d#(m ^2 + varA^2) / theta^2 #
  crowding_effect <- as.vector(pars$eta %*% n)
  Q<- 1 / (1 + pars$c * (crowding_effect)^pars$gamma)
  b<- 1- (m ^2 + varA)/theta^2    #0.5*(erf((theta-m)/sqrt(2*varA))+erf((theta+m)/sqrt(2*varA)))
  g<-  -2 * m * varA / theta^2   #sqrt(varA)/sqrt(2*pi)*exp(-(theta+m)^2/(2*varA))*(1-exp(2*theta*m/varA))
  dmA <- outer(m, m, FUN="-") ## difference matrix of trait means
  svA <- outer(varA, varA, FUN="+") ## sum matrix of trait variances
  alpha <- exp(-dmA^2/(2*svA+pars$w^2))*pars$w/sqrt(2*svA+pars$w^2) ## alpha matrix
   beta<- alpha*2*varA*(-dmA)/(2*svA+pars$w^2)^1.5
  
  alpha_eff <- as.vector(alpha %*% n)
  beta_eff  <- as.vector(beta %*% n) 
  
  dndt <- n*(b*Q-d- Q^2*alpha_eff )*cutoff(n/(1e-4)) ##  ddensity eqs
  alive_threshold <- ifelse(n > 1e-5, 1, 0)
  dmdt <- pars$h2*(g*Q -Q^2*beta_eff)*alive_threshold  #trait equation
  
  return(list(c(dndt, dmdt)))
}



# trait-averaged Q quantities for species i 
# Computes <Q>_i, <Q^2>_i, <b(z) Q(z,N)>_i, <g(z) Q(z,N)>_i
# where b(z) = 1 - z^2/theta^2  and  g(z) = (z - u_i) * b(z) restoring term
# is folded into the trait equation through a separate quadrature.

trait_Q_quantities <- function(i_focal, m_vec, var_vec, N, c_i,
                               eta_w, theta, points = 9) {
  S <- length(N)
  
  outer_q <- gauss.quad.prob(points, dist = "normal",
                             mu    = m_vec[i_focal],
                             sigma = sqrt(var_vec[i_focal]))
  z1 <- outer_q$nodes
  w1 <- outer_q$weights
  
  z2_list <- lapply(1:S, function(j)
    gauss.quad.prob(points, dist = "normal",
                    mu    = m_vec[j],
                    sigma = sqrt(var_vec[j])))
  
  Q_at_z1  <- numeric(points)
  
  for (k in seq_len(points)) {
    z_focal <- z1[k]
    crowding_sum <- 0
    for (j in 1:S) {
      eta_kernel <- exp(-(z_focal - z2_list[[j]]$nodes)^2 / eta_w^2)
      crowding_sum <- crowding_sum +
        N[j] * sum(eta_kernel * z2_list[[j]]$weights)
    }
    Q_at_z1[k] <- 1 / (1 + c_i * sqrt(pmax(crowding_sum, 0)))
  }
  
  b_at_z1 <- 1 - z1^2 / theta^2
  # selection gradient piece for trait equation: (z - u_i) * b(z) * Q(z,N)
  trait_grad_at_z1 <- (z1 - m_vec[i_focal]) * b_at_z1 * Q_at_z1
  
  list(
    Q_avg          = sum(Q_at_z1        * w1),
    Q2_avg         = sum(Q_at_z1^2      * w1),
    bQ_avg         = sum(b_at_z1 * Q_at_z1 * w1),
    b_avg          = sum(b_at_z1 * w1),
    trait_grad_avg = sum(trait_grad_at_z1 * w1)
  )
}

# main eco-evolutionary equations

eqs_trait_based_crowding <- function(time, state, pars) {
  S      <- length(pars$sigma)
  n      <- state[1:S]
  m_vec  <- state[(S+1):(2*S)]
  varA   <- pars$sigma
  w      <- pars$w
  theta  <- pars$theta
  eta_w  <- pars$eta_w
  c_vec  <- pars$c
  pts    <- pars$points %||% 9
  

  dmA <- outer(m_vec, m_vec, FUN = "-")
  svA <- outer(varA, varA, FUN = "+")
  alpha <- exp(-dmA^2 / (2*svA + w^2)) * w / sqrt(2*svA + w^2)
  beta  <- alpha * 2 * varA * (-dmA) / (2*svA + w^2)^1.5
  
  
  Q_avg  <- numeric(S)
  Q2_avg <- numeric(S)
  bQ_avg <- numeric(S)
  trait_grad_Q <- numeric(S)   # <(z - u_i) b(z) Q(z,N)>_i
  if (all(c_vec == 0)) {
    Q_avg  <- rep(1, S)
    Q2_avg <- rep(1, S)
    bQ_avg <- 1 - (m_vec^2 + varA) / theta^2
    trait_grad_Q <- -2 * m_vec * varA / theta^2
  } else {
    for (i in 1:S) {
      q <- trait_Q_quantities(i, m_vec, varA, n,
                              c_i   = c_vec[i],
                              eta_w = eta_w,
                              theta = theta,
                              points = pts)
      Q_avg[i]        <- q$Q_avg
      Q2_avg[i]       <- q$Q2_avg
      bQ_avg[i]       <- q$bQ_avg
      trait_grad_Q[i] <- q$trait_grad_avg
    }
  }
  
  
  
  d <- pars$d
  
  
  # The trait equation has <- Q^2> * (beta %*% n) for the competition pull,
  # using the species-averaged <Q^2>_i.
  competition_trait_pull <- Q2_avg * as.vector(beta %*% n)
  
 
  # dn_i/dt = n_i * ( <b Q>_i  -  d_i  -  <Q^2>_i * (alpha %*% n)_i )
  dndt <- n * (bQ_avg - d - Q2_avg * as.vector(alpha %*% n)) *
    cutoff(n / 1e-6)
  
  # 
  # du_i/dt = h^2 * ( <(z - u_i) b(z) Q(z,N)>_i  -  <Q^2>_i * (beta %*% n)_i )
  dmdt <- pars$h2 * (trait_grad_Q - competition_trait_pull)
  
  return(list(c(dndt, dmdt)))
}

# null-coalescing operator for default arg
`%||%` <- function(a, b) if (!is.null(a)) a else b




organize_results_ode <- function(sol, pars) {
  S <- length(pars$sigma) ## number of species
  dat <- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  names(dat)[1] <- "time" ## name the first column "time"
  names(dat)[2:(S+1)] <- paste0("n_", 1:S) ## name abundance columns (n_k)
  names(dat)[(S+2):(2*S+1)] <- paste0("m_", 1:S) ## name trait mean columns
  dat <- dat %>%
    gather("variable", "v", 2:ncol(dat)) %>% ## normalize the data
    separate(variable, c("type", "species"), sep="_") %>%
    spread(type, v) %>% ## separate columns for densities n and trait means m
    dplyr::select(time, species, n, m) %>% ## rearrange columns
    mutate(species=as.integer(species), sigma=pars$sigma[species], w=pars$w,theta=pars$theta) ## add parameters
  return(as_tibble(dat))
}


# when triggered, snap below-threshold species to exactly zero
extinction_event <- function(t, state, pars) {
  state[state < 1e-8] <- 0
  state
}

#zero-crossing detector — triggers when any species hits 1e-8
extinction_root <- function(t, state, pars) {
  state - 1e-8
}

run_to_eq <- function(eta_mat, a_mat, b_vec, m0_vec, c_vec, gamma,
                      ic = c(1, 1), tmax = 5e3) {
  pars <- list(c = c_vec, eta = eta_mat, a = a_mat,
               m_0 = m0_vec, b = b_vec, gamma=gamma)
  out <- tryCatch(
    ode(func = sublinear_eqs, y = ic, parms = pars,
        times = c(0, tmax),
        method = "lsoda", rtol = 1e-8, atol = 1e-10,
        events = list(func = extinction_event, root = TRUE),
        rootfun = extinction_root),
    error = function(e) NULL
  )
  if (is.null(out)) return(c(N1 = NA, N2 = NA))
  final <- out[nrow(out), -1]
  c(N1 = as.numeric(final[1]), N2 = as.numeric(final[2]))
}

cclassify <- function(N1, N2, thr = 1e-6) {
  dplyr::case_when(
    N1 > thr & N2 > thr   ~ "Coexistence",
    N1 > thr & N2 <= thr  ~ "Sp 1 survives only",
    N1 <= thr & N2 > thr  ~ "Sp 2 survives only",
    TRUE                  ~ "Both extinct"
  )
}
outcome_colors <- c("Coexistence"  = "#2ca25f",
                    "Sp 1 survives only"    = "grey70",# "#3182bd",
                    "Sp 2 survives only"    ="grey70",# "#de2d26",
                    "Both extinct" = "grey70")#"grey70")

plot_outcome <- function(grid, xvar, yvar, xlab, ylab, title, subtitle) {
  ggplot(grid, aes(x = .data[[xvar]], y = .data[[yvar]], fill = outcome)) +
    geom_raster(interpolate = FALSE) +
    scale_fill_manual(values = c("Coexistence"  = "#2ca25f",
                                 "Sp 1 survives only"    =  "grey70", #"#3182bd",
                                 "Sp 2 survives only"    = "grey70", #, "#de2d26",
                                 "Both extinct" = "grey70")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = xlab, y = ylab, fill = "Outcome",
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
}



plot_snapshot <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
  S <- dat %>% pull(species) %>% max ## number of species
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  snap <- dat %>% filter(time==moment) %>% dplyr::select(-time) ## time = moment
  traits <- expand.grid(species=1:S, trait=traitaxis) %>% as_tibble ## trait table
  traits["density"] <- 0 ## add column for population densities
  for (i in 1:S) {
    v <- snap %>% filter(species==i) %>% dplyr::select(n, m, sigma)
    traits$density[(traits$species==i)] <- v$n * ## trait normally distributed,
      dnorm(traits$trait[(traits$species==i)], v$m, v$sigma) ## times density
  }
  landscape <- tibble(trait=traitaxis, r=0) ## for plotting intrinsic rates
  landscape$r <- 1 - (landscape$trait^2+0)/dat$theta[1]^2 #   1-( u[t,k,]^2+sigma_i[k,]^2)/theta[1]^2 
  #(sign(dat$theta[1]+landscape$trait)+1)*  ## phenotype-specific
  #(sign(dat$theta[1]-landscape$trait)+1)/4              ## growth rates
  # scale growth rates to show nicely on plot:
   landscape$r <- landscape$r*max(traits$density,na.rm=TRUE)
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15) +
    geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
             colour="darkred", alpha=0.9, na.rm=TRUE) +
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    # scale_color_viridis_D(alpha = 1)+
    scale_fill_viridis_d()+
    theme_bw()+
    ggtitle(label = bquote(t == .(moment)))+
    theme(legend.position="none",text = element_text(size=12)) %>%
    return
}

## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
 plot_grid(plot_density(dat), plot_snapshot(dat, moment=0, limits, res),
           plot_snapshot(dat, moment, limits, res), ncol=1, align="hv") %>%
    
    #plot_grid( plot_snapshot(dat, moment=0, limits, res),
     #         plot_snapshot(dat, moment, limits, res), ncol=1, align="hv") %>%
    return
}

## Output:
## - a ggplot2 plot
plot_density <- function(dat) {
  dat %>%
    ggplot +
    theme_bw()+
    geom_line(aes(x=time, y=n, colour=as.factor(species)),linetype = "dashed",size=1.1) +
    scale_y_continuous(name="population \n density", limits=c(0, NA)) +
    theme(legend.position="none") %>%
    return
}


