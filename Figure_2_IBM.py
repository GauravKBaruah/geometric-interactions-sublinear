#Individual-based model simulating feeding/foraging phase with non-consumptive interference
#Figure 02
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from joblib import Parallel, delayed
from tqdm import tqdm


def simulate_feeding_phase(
    N, xi, arena_radius=10.0, larva_radius=0.45, t_max=500,
    feed_rate=0.1, move_speed=0.5, penalty_time=8,
    m0=5, conversion_efficiency=0.2,conversion_efficiency_sd=0.01, total_resource=2500.0, r_int=1.0
):
    """
    Simulates consumer feeding in a bounded arena with social interaction rules.
    Returns the realized per-capita growth rate and the mean feeding time fraction Q(N).
    """
    r = np.sqrt(np.random.uniform(0, 1, N)) * (arena_radius - larva_radius)
    theta = np.random.uniform(0, 2 * np.pi, N)
    pos = np.column_stack((r * np.cos(theta), r * np.sin(theta)))

    energy = np.zeros(N)
    timeout = np.zeros(N, dtype=int)
    feeding_time = np.zeros(N, dtype=int)
    speeds = np.abs(np.random.normal(move_speed, 1, N))
    individual_conversion_efficiency = np.random.normal(conversion_efficiency, conversion_efficiency_sd, N)
    individual_conversion_efficiency = np.clip(individual_conversion_efficiency, 0, 1.0)
    
    for t in range(t_max):
        timeout[timeout > 0] -= 1
        active = timeout == 0

        if np.any(active):
            active_idx = np.where(active)[0]
            num_active = len(active_idx)
            
            angles = np.random.uniform(0, 2 * np.pi, num_active)
            active_pos = pos[active_idx]
            
            if xi != 0:
                diff_all = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
                dists_all = np.linalg.norm(diff_all, axis=2)
                np.fill_diagonal(dists_all, np.inf)

                in_range = dists_all < (r_int + larva_radius)
                
                for i in range(num_active):
                    global_i = active_idx[i]
                    dists_from_i = dists_all[global_i]
                    nearest_neighbour = np.argmin(dists_from_i)
                    vector_to_neighbor = pos[nearest_neighbour] - active_pos[i]
                    norm = np.linalg.norm(vector_to_neighbor)
                    if norm > 0:
                        if xi < 0:
                            target_vector = -vector_to_neighbor
                        else:
                            if norm < (2 * larva_radius):
                                target_vector = -vector_to_neighbor
                            else:
                                target_vector = vector_to_neighbor
                                
                        social_theta = np.arctan2(target_vector[1], target_vector[0])
                        current_theta = angles[i]
                        
                        delta_theta = np.arctan2(
                            np.sin(social_theta - current_theta), 
                            np.cos(social_theta - current_theta)
                        )
                        turn_amount = abs(xi) * delta_theta
                        turn_amount = np.clip(turn_amount, -np.abs(delta_theta), np.abs(delta_theta))
                        
                        angles[i] = current_theta + turn_amount
                
            new_x = active_pos[:, 0] + speeds[active_idx] * np.cos(angles)
            new_y = active_pos[:, 1] + speeds[active_idx] * np.sin(angles)
            new_pos = np.column_stack((new_x, new_y))

            dist_to_center = np.linalg.norm(new_pos, axis=1)
            valid_move = dist_to_center <= (arena_radius - larva_radius)
            invalid_move = ~valid_move
            angles[invalid_move] += np.pi  
            pos[active_idx[valid_move]] = new_pos[valid_move]

        diff_all = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
        dists_all = np.linalg.norm(diff_all, axis=2)
        np.fill_diagonal(dists_all, np.inf)

        bumped = np.any(dists_all < (2 * larva_radius), axis=1)
        
        newly_interrupted = active & bumped
        timeout[newly_interrupted] = penalty_time

        successfully_fed = active & ~bumped
        num_feeding = np.sum(successfully_fed)
        feeding_time[successfully_fed] += 1

        if num_feeding > 0 and total_resource > 0:
            total_demand = num_feeding * feed_rate
            
            if total_resource >= total_demand:
                energy[successfully_fed] += feed_rate
                total_resource -= total_demand
            else:
                fair_share = total_resource / num_feeding
                energy[successfully_fed] += fair_share
                total_resource = 0

    net_energy = np.maximum(energy - m0, 0)
    mean_offspring = net_energy * conversion_efficiency
    
    offspring = np.random.poisson(mean_offspring)
    per_capita_growth = np.mean(offspring)
    mean_feeding_time = np.mean(feeding_time)
    Q = mean_feeding_time / t_max
    
    return per_capita_growth, Q

# Analytical Functions
def interference_function(N, c, gamma):
    return 1.0 / (1.0 + c * (N ** gamma))

def analytical_per_capita_growth(N, c, gamma, b=1.0, m0=0.1, a=0.006):
    """
    Derived from Eq 5: (1/N)(dN/dt) = b*Q(N) - m0 - a*Q(N)^2*N
    """
    Q = interference_function(N, c, gamma)
    return b * Q - m0 - a * (Q ** 2) * N
    
#Simulate with and without interaction

def simulate_and_plot_xi_range(xi_values):
    N_values = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150])
    replicates = 50
    
    all_raw_growths = []
    all_raw_Qs = []
    
    all_gamma_estimates = []
    all_mean_growths = []
    all_std_growths = [] 
        
    for xi in tqdm(xi_values, desc="Overall Progress (xi)", position=0):
        
        results = Parallel(n_jobs=-1)(
            delayed(simulate_feeding_phase)(N, xi=xi) 
            for N in N_values for _ in range(replicates)
        )
        
        growth_rates = np.array([res[0] for res in results]).reshape(len(N_values), replicates)
        Q_values = np.array([res[1] for res in results]).reshape(len(N_values), replicates)
        
        all_raw_growths.append(growth_rates)
        all_raw_Qs.append(Q_values)
        
        all_mean_growths.append(np.mean(growth_rates, axis=1))
        all_std_growths.append(np.std(growth_rates, axis=1))
        
        gamma_estimates = np.full(replicates, np.nan)
        
        for rep in range(replicates):
            Q_curve = Q_values[:, rep] 
            
            try:
                popt, _ = curve_fit(
                    interference_function, 
                    N_values, 
                    Q_curve, 
                    p0=[0.1, 0.5], 
                    bounds=([0.0, 0.0], [np.inf, np.inf]),
                    maxfev=2000
                )
                gamma_estimates[rep] = popt[1]  # popt is [c_fit, gamma_fit]
            except RuntimeError:
                pass 
                
        all_gamma_estimates.append(gamma_estimates)
        successful_fits = np.sum(~np.isnan(gamma_estimates))
    
    # Execute IBM with No Interference (penalty_time = 0)
    print("Running IBM No Interference (Penalty Time = 0)...")
    results_no_int = Parallel(n_jobs=-1)(
        delayed(simulate_feeding_phase)(N, xi=0.0, penalty_time=0) 
        for N in N_values for _ in range(replicates)
    )
    
    growth_rates_no_int = np.array([res[0] for res in results_no_int]).reshape(len(N_values), replicates)
    Q_values_no_int = np.array([res[1] for res in results_no_int]).reshape(len(N_values), replicates)

    # Raw observations and linear fit for the penalty_time = 0 panel.
    # Jitter is used only when plotting; the model is fitted to the true N values.
    raw_N_no_int = np.repeat(N_values, replicates)
    raw_growth_no_int = growth_rates_no_int.reshape(-1)
    linear_slope, linear_intercept = np.polyfit(
        raw_N_no_int,
        raw_growth_no_int,
        deg=1
    )
    fitted_raw_growth = linear_intercept + linear_slope * raw_N_no_int
    residual_sum_squares = np.sum(
        (raw_growth_no_int - fitted_raw_growth) ** 2
    )
    total_sum_squares = np.sum(
        (raw_growth_no_int - np.mean(raw_growth_no_int)) ** 2
    )
    linear_r_squared = (
        1.0 - residual_sum_squares / total_sum_squares
        if total_sum_squares > 0
        else np.nan
    )
    linear_N = np.linspace(N_values.min(), N_values.max(), 200)
    linear_growth = linear_intercept + linear_slope * linear_N

    jitter_rng = np.random.default_rng(12345)
    jittered_N_no_int = raw_N_no_int + jitter_rng.uniform(
        -0.8,
        0.8,
        size=raw_N_no_int.size
    )


    # Plotting 
    plt.style.use('ggplot') 
    fig, axs = plt.subplots(2, 2, figsize=(16, 12))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(xi_values)))
    N_continuous = np.linspace(0, 150, 200)

    # 1. Top-Left (0, 0): Analytical Growth Function
    axs[0, 0].plot(N_continuous, analytical_per_capita_growth(N_continuous, c=0.0, gamma=1.0), 
                   label='Logistic ($c=0$)', color='crimson', linewidth=2)
    axs[0, 0].plot(N_continuous, analytical_per_capita_growth(N_continuous, c=0.15, gamma=0.5), 
                   label=r'Sublinear ($\gamma=0.5$)', color='darkblue', linewidth=2)
    axs[0, 0].plot(N_continuous, analytical_per_capita_growth(N_continuous, c=0.15, gamma=1.0), 
                   label=r'Sublinear ($\gamma=1.0$)', color='steelblue', linewidth=2)
    axs[0, 0].set_xlabel('Population Density (N)', fontsize=12)
    axs[0, 0].set_ylabel('per-capita growth rate', fontsize=12)
    axs[0, 0].legend(loc='upper right')

    # 2. Top-Right (0, 1): Raw IBM observations and fitted linear model
    # for no persistent interference (penalty_time = 0).
    axs[0, 1].scatter(
        jittered_N_no_int,
        raw_growth_no_int,
        color='crimson',
        s=18,
        alpha=0.25,
        edgecolors='none',
        label='IBM replicates'
    )
    axs[0, 1].plot(
        linear_N,
        linear_growth,
        color='black',
        linewidth=2.2,
        label=rf'Linear fit ($R^2={linear_r_squared:.3f}$)'
    )
    axs[0, 1].set_xlabel('Population Density (N)', fontsize=12)
    axs[0, 1].set_ylabel('Realized per-capita growth rate', fontsize=12)
    axs[0, 1].legend(loc='upper right')

    # 3. Bottom-Left (1, 0): Per-Capita Growth Rates across N for each xi
    for i, xi in enumerate(xi_values): 
        axs[1, 0].errorbar(N_values, all_mean_growths[i], yerr=all_std_growths[i], fmt='o', 
                           color=colors[i], label=rf'$\xi = {xi}$', 
                           markersize=4, capsize=3, alpha=0.8)
    
    axs[1, 0].set_xlabel('Population Density (N)', fontsize=12)
    axs[1, 0].set_ylabel('Realized per-capita growth rate', fontsize=12)
    axs[1, 0].legend(loc='upper right', fontsize=10)
    
    # 4. Bottom-Right (1, 1): Boxplot of Gamma estimates
    clean_gamma_estimates = [gamma[~np.isnan(gamma)] for gamma in all_gamma_estimates]
    
    box = axs[1, 1].boxplot(clean_gamma_estimates, orientation='vertical', patch_artist=True, tick_labels=xi_values)
    
    for patch in box['boxes']:
        patch.set_facecolor('steelblue')
        patch.set_alpha(0.8)
    for median in box['medians']:
        median.set(color='darkorange', linewidth=2.5)
        
    axs[1, 1].axhline(0.5, color='grey', linestyle='--', linewidth=1.5, label=r'Nearest-Neighbour Attraction($\gamma = 0.5$)')
    axs[1, 1].axhline(1.0, color='black', linestyle=':', linewidth=1.5, label=r'Well-Mixed ($\gamma = 1.0$)')
        
    axs[1, 1].set_xlabel(r'Interaction value ($\xi$)', fontsize=12)
    axs[1, 1].set_ylabel(r'Estimated ($\gamma$)', fontsize=12)
    axs[1, 1].legend(loc='best')
    
    # Aesthetics
    for ax in axs.flat:
        ax.grid(False) 
        ax.set_facecolor('white') 
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        if ax != axs[1, 1]: 
            ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)

    plt.tight_layout()
    plt.savefig("Figure_02.pdf")
    plt.show()

# Execute
if __name__ == "__main__":
    xi_values = [-1.0,-0.75,-0.5, 0.0,0.5 ,0.75,1.0]
    simulate_and_plot_xi_range(xi_values)
