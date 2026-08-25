library(cowplot)
library(statmod)

#error function
erf <- function(x) 2*pnorm(x*sqrt(2))-1
#cutoff function for smoothing extinctions
cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)

#trait Cv - takes equilibrium trait values and returns the CV - high CV high trait clustering
traitCV <- function(mus){
  
  v <- sort(mus[is.na(mus)==FALSE])
  d <- diff(v)
  cvs<-sd(d,na.rm=T)/mean(d,na.rm=T)
  cvs[is.na(cvs)]<-0
  
  return(cvs)
}


#function to classify which species survive, go extinct or coexist
#N1 - density of species 1
#N2 - density of species 2
#thr- extinction threshold
classify <- function(N1, N2, thr = 1e-6) {
  dplyr::case_when(
    N1 > thr & N2 > thr   ~ "Coexistence",
    N1 > thr & N2 <= thr  ~ "Sp 1 survives",
    N1 <= thr & N2 > thr  ~ "Sp 2 survives",
    TRUE                  ~ "Both extinct"
  )
}

#color scheme used to produce MCT plots

outcome_colors <- c("Coexistence"  = "#2ca25f",
                    "Sp 1 survives only"    = "grey70",# "#3182bd",
                    "Sp 2 survives only"    ="grey70",# "#de2d26",
                    "Both extinct" = "grey70")#"grey70")

#subliner dynamical function
# sublinear dynamics withQ_i*Q_j competition structure
sublinear_eqs <- function(t, state, pars) {
  N <- pmax(state, 0)
  N[N < 1e-10] <- 0
  crowding_effect <- as.vector(pars$eta %*% N)
  Q <- 1/(1+pars$c*(crowding_effect)^pars$gamma)
  competition <- as.vector(pars$a %*% (Q*N))
  dndt <- N*(pars$b*Q - pars$m_0 - Q*competition)
  list(c(dndt))
}


#dynamical equations and function that goes into the ODE solver that either has sublinear c>0, or LV, c=0
#time - time 
#state - states of the community, N, and mus (trat values)
#pars- parameter list
eqs <- function(time, state, pars) {
  S <- length(pars$sigma) ## number of species
  n <- state[1:S] ## species densities
  m <- state[(S+1):(2*S)] ## species trait means
  varA <- pars$sigma ## species' trait standard deviations
  w<-pars$w
  theta<-pars$theta
  d<- pars$d #mortality rate
  crowding_effect <- as.vector(pars$eta %*% n)
  Q<- 1 / (1 + pars$c * (crowding_effect)^pars$gamma)
  b<- 1- (m ^2 + varA)/theta^2    
  g<-  -2 * m * varA / theta^2   
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




`%||%` <- function(a, b) if (!is.null(a)) a else b


#for ploting this function is used to organise the ODE output into a nice tibble
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


#when happens, it snaps below-threshold species to exactly zero
extinction_event <- function(t, state, pars) {
  state[state < 1e-8] <- 0
  state
}
extinction_root <- function(t, state, pars) {
  state - 1e-8
}


#function used to simulate dynamics of LV and sublinear regime (c>0) and used to produce the MCT plots
#eta_mat - interference matrix
#a_mat - competition matrix
#b_vec - vector of growht rates
#m0_vec- vector of mortality rates
#c_vec - vector of species interference sensitivity , c=0, or c = 1(sublinear)
#gamma - interference exponent, either 0.5, or 1
#tmax - max time steps
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


#plotting functions
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


#plotting functions :NOT used in the MS but could be useful
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


