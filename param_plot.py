import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.colors as mcolors
import os
from scipy.interpolate import griddata
from matplotlib.cm import get_cmap
import seaborn as sns
from numpy.polynomial.polynomial import Polynomial
#from sklearn.model_selection import train_test_split
#from sklearn.preprocessing import StandardScaler
#from tensorflow.keras.models import Sequential
#from tensorflow.keras.layers import Dense
#from tensorflow.keras.callbacks import EarlyStopping
from scipy.stats import binned_statistic
from scipy.optimize import curve_fit

#from sklearn.gaussian_process import GaussianProcessRegressor
#from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

# Load data from file
def load_data(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            values = list(map(float, line.strip().split('\t')))
            data.append(values)
    return np.array(data)

def sigmoid_fixed_amp(x, A, B):
    """
    Sigmoid with fixed amplitude 14:
        f(x) = 14 / (1 + exp(-A (x - B)))
    """
    return 14.0 / (1.0 + np.exp(-A * (x - B)))

def sigmoid_14_asym(x, a, b, c):
    """
    Sigmoid that asymptotes to 14 as x -> +inf:
        f(x) = c / (1 + exp(-a (x - b))) + 14 - c

    As x -> -inf: f -> 14 - c
    As x -> +inf: f -> 14
    """
    return c / (1.0 + np.exp(-a * (x - b))) + 14.0 - c


# Train a neural network to model decay time
def train_decay_time_model(data, feature_cols, target_col):
    X = data[:, feature_cols]
    y = data[:, target_col]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    model = Sequential([
        Dense(64, activation='relu', input_shape=(X.shape[1],)),
        Dense(64, activation='relu'),
        Dense(1)
    ])

    model.compile(optimizer='adam', loss='mse')

    early_stop = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

    history = model.fit(X_train, y_train, validation_split=0.2, epochs=100, batch_size=16, callbacks=[early_stop], verbose=1)

    loss = model.evaluate(X_test, y_test)
    print(f"Test MSE: {loss:.4f}")

    return model, scaler

def train_gp_model(data, feature_cols, target_col):
    X = data[:, feature_cols]
    y = data[:, target_col]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    kernel = C(1.0, (1e-3, 1e3)) * RBF(length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-2, 1e2))
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=1e0, normalize_y=True)

    gp.fit(X_train, y_train)
    score = gp.score(X_test, y_test)
    print(f"Test R^2 Score: {score:.4f}")

    return gp, scaler

def plot_decay_vs_angle_kernel_bands_sig(
    data,
    columns,
    fixed_param='fg',
    fixed_value=0.5,
    tolerance=1e-2,
    bandwidth=5.0
):
    """
    Plots decay time vs. angle with a sigmoid fit of the form

        f(x) = c / (1 + exp(-a (x - b))) + 14 - c

    and adaptive ±1σ, ±2σ bands using kernel-smoothed residual std.
    """

    # Filter data at fixed parameter
    col_index = columns[fixed_param]
    mask = np.abs(data[:, col_index] - fixed_value) < tolerance
    filtered_data = data[mask]

    angles = filtered_data[:, columns['angle']]
    decay_times = filtered_data[:, columns['decay_time']]

    angle_range = np.linspace(np.min(angles), np.max(angles), 300)

    # ---------- Sigmoid Fit (asymptote at 14) ----------

    amin, amax = np.min(angles), np.max(angles)
    span = amax - amin if amax > amin else 1.0

    # Initial guesses:
    # lower asymptote ~ min(decay_times) ≈ 14 - c  =>  c ≈ 14 - min
    dt_min = float(np.min(decay_times))
    dt_max = float(np.max(decay_times))

    c0 = 14.0 - dt_min
    # if that does something crazy, you can clamp it later if needed

    b0 = 0.5 * (amin + amax)        # midpoint of angle range
    a0 = 4.0 / span                 # steepness guess

    p0 = [a0, b0, c0]

    popt, pcov = curve_fit(
        sigmoid_14_asym,
        angles,
        decay_times,
        p0=p0,
        maxfev=10000
    )
    a_fit, b_fit, c_fit = popt

    # --- Print coefficients ---
    print("\nSigmoid fit coefficients for f(x) = c / (1 + exp(-a (x - b))) + 14 - c:")
    print(f"  a (steepness)      = {a_fit:.6g}")
    print(f"  b (midpoint angle) = {b_fit:.6g}")
    print(f"  c (amplitude)      = {c_fit:.6g}")
    print(f"  Lower asymptote    ≈ {14.0 - c_fit:.6g}")
    print(f"  Upper asymptote    → 14\n")

    # Compute fit values
    fit_vals = sigmoid_14_asym(angle_range, a_fit, b_fit, c_fit)

    # Residuals
    residuals = decay_times - sigmoid_14_asym(angles, a_fit, b_fit, c_fit)

    # Kernel smoothed standard deviation
    local_std = kernel_smooth_std(angles, residuals, angle_range, bandwidth)

    # Optional LaTeX equation string
    equation_str = (
        r"$t_{\rm p}(i_0) = \frac{" + f"{c_fit:.3g}" +
        r"}{1 + \exp\left[-" + f"{a_fit:.3g}" + r"(i_0 - " +
        f"{b_fit:.3g}" + r")\right]} + 14 - " + f"{c_fit:.3g}$"
    )
    # print(equation_str)

    # ---------- Plot ----------
    plt.figure(figsize=(9, 5))

    # Data
    plt.scatter(
        angles,
        decay_times,
        color='blue',
        alpha=0.7,
        edgecolor='k',
        s=30,
        label='Data'
    )

    # Sigmoid fit
    plt.plot(
        angle_range,
        fit_vals,
        color='darkred',
        linewidth=2.5,
        label='Sigmoid fit'
    )

    # ±1σ and ±2σ bands
    for n, color, alpha in zip([1, 2], ['orange', 'gold'], [0.3, 0.2]):
        upper = fit_vals + n * local_std
        lower = fit_vals - n * local_std
        plt.fill_between(
            angle_range,
            lower,
            upper,
            color=color,
            alpha=alpha,
            label=fr'$\pm{n}\sigma$' if n == 1 else None
        )

    plt.xlabel(r'$i_{\rm 0}$ (deg)', fontsize=16)
    plt.ylabel(r'$t_{\rm p}$ (Gyr)', fontsize=16)
    plt.ylim([0, 14.7])

    # If you want the equation on the plot:
    # plt.text(
    #     0.95, 0.05, equation_str, transform=plt.gca().transAxes,
    #     fontsize=13, verticalalignment='bottom', horizontalalignment='right',
    #     bbox=dict(facecolor='white', alpha=0.7)
    # )

    plt.grid(True)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    #plt.legend(fontsize=12)
    plt.savefig(
        '/Users/sghob/OneDrive/Desktop/Res_code/108.png',
        format='png',
        bbox_inches='tight'
    )
    plt.show()
        
# Plot 3D scatter and neural net prediction surface
def plot_3d_prediction_surface_nn(data, model, scaler, columns):
    # Choose x_coord and angle as independent variables, fix the rest
    fixed_vals = {
        'mass': 1e6,
        'vg': -0.7,
        'ngd': 300.0
    }

    x_vals = np.unique(data[:, columns['fg']])
    angle_vals = np.unique(data[:, columns['angle']])
    xx, aa = np.meshgrid(x_vals, angle_vals)

    input_grid = []
    for x, a in zip(np.ravel(xx), np.ravel(aa)):
        sample = [
            fixed_vals['mass'],
            x,
            fixed_vals['vg'],
            fixed_vals['ngd'],
            a
        ]
        input_grid.append(sample)

    input_grid = np.array(input_grid)
    input_scaled = scaler.transform(input_grid)
    predicted = model.predict(input_scaled).flatten()
    predicted_grid = predicted.reshape(xx.shape)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the original scatter data
    mask = (data[:, columns['mass']] == fixed_vals['mass']) & \
           (data[:, columns['vg']] == fixed_vals['vg']) & \
           (data[:, columns['ngd']] == fixed_vals['ngd'])

    scatter = data[mask]
    ax.scatter(scatter[:, columns['fg']], scatter[:, columns['angle']], scatter[:, columns['decay_time']], color='red', label='Actual Data')

    # Plot the surface
    ax.plot_surface(xx, aa, predicted_grid, cmap='viridis', alpha=0.6, edgecolor='none')
    ax.set_xlabel('fg')
    ax.set_ylabel('angle')
    ax.set_zlabel('Decay Time')
    ax.set_title('Neural Net Prediction Surface vs Actual Data')
    ax.legend()
    plt.tight_layout()
    plt.show()

def plot_3d_prediction_surface(data, model, scaler, columns):
    fixed_vals = {
        'mass': 1E6,
        'vg': -0.7,
        'ngd': 300
    }

    x_vals = np.unique(data[:, columns['fg']])
    angle_vals = np.unique(data[:, columns['angle']])
    xx, aa = np.meshgrid(x_vals, angle_vals)

    input_grid = []
    for x, a in zip(np.ravel(xx), np.ravel(aa)):
        sample = [
            fixed_vals['mass'],
            x,
            fixed_vals['vg'],
            fixed_vals['ngd'],
            a
        ]
        input_grid.append(sample)

    input_grid = np.array(input_grid)
    input_scaled = scaler.transform(input_grid)
    predicted = model.predict(input_scaled).flatten()
    predicted_grid = predicted.reshape(xx.shape)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    mask = (data[:, columns['mass']] == fixed_vals['mass']) & \
           (data[:, columns['vg']] == fixed_vals['vg']) & \
           (data[:, columns['ngd']] == fixed_vals['ngd'])

    scatter = data[mask]
    ax.scatter(scatter[:, columns['fg']], scatter[:, columns['angle']], scatter[:, columns['decay_time']], color='red', label='Actual Data')

    ax.plot_surface(xx, aa, predicted_grid, cmap='viridis', alpha=0.6, edgecolor='none')
    ax.set_xlabel('fg')
    ax.set_ylabel('angle')
    ax.set_zlabel('Decay Time')
    ax.set_title('Gaussian Process Prediction Surface vs Actual Data')
    ax.legend()
    plt.tight_layout()
    plt.show()

# Plot decay time against two chosen parameters in 3D scatter plot
def plot_3d_scatter(data, x_col, y_col, z_col, fixed_params):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)
    
    filtered_data = data[mask]
    #x = np.log10(filtered_data[:, x_col])
    x = filtered_data[:, x_col]
    y = filtered_data[:, y_col]
    z = filtered_data[:, z_col]
    
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=z, cmap='viridis', marker='o')
    
    #ax.set_xlabel(r'$\log{\frac{M}{M_{\odot}}}$')
    ax.set_xlabel(r'$n_{gd}$')
    ax.set_ylabel(r'Inclination ($\degree$)')
    ax.set_zlabel('Decay Time (Gyr)')
    ax.set_title('3D Decay Time Visualization')
    
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=10)
    cbar.set_label('Decay Time (Gyr)')
    
    plt.show()

# Plot decay time using a 2D color representation
def plot_2d_heatmap(data, x_col, y_col, z_col, fixed_params):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)
    
    filtered_data = data[mask]
    #x = np.log10(filtered_data[:, x_col])
    x = filtered_data[:, x_col]
    y = filtered_data[:, y_col]
    z = filtered_data[:, z_col]
    
    grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')
    
    plt.figure(figsize=(8, 6))
    plt.contourf(grid_x, grid_y, grid_z, levels=100, cmap='viridis')
    plt.colorbar(label='Decay Time (Gyr)')
    #plt.xlabel(r'$\log{\frac{M}{M_{\odot}}}$')
    plt.xlabel(r'$n_{gd}$')
    plt.ylabel(r'Inclination ($\degree$)')
    plt.show()

# Plot decay time using a 2D discrete color representation
def plot_2d_discrete_heatmap(data, x_col, y_col, z_col, fixed_params, num_levels=10):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)
    
    filtered_data = data[mask]
    #x = np.log10(filtered_data[:, x_col])  # Apply log to mass
    x = filtered_data[:, x_col]
    y = filtered_data[:, y_col]
    z = filtered_data[:, z_col]
    
    grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')
    
    plt.figure(figsize=(8, 6))
    levels = np.linspace(np.min(z), np.max(z), num_levels)  # Discrete levels
    plt.contourf(grid_x, grid_y, grid_z, levels=levels, cmap='viridis')
    plt.colorbar(label='Decay Time (Gyr)')
    #plt.xlabel(r'$\log{\frac{M}{M_{\odot}}}$')
    plt.xlabel(r'$n_{gd}$')
    plt.ylabel(r'Inclination ($\degree$)')
    plt.show()


# Plot decay time using a binned scatter plot with a continuous colorbar
def plot_binned_scatter(data, x_col, y_col, z_col, fixed_params):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)
    
    filtered_data = data[mask]
    #x = np.log10(filtered_data[:, x_col])  # Apply log to mass
    x = filtered_data[:, x_col]
    y = filtered_data[:, y_col]
    z = filtered_data[:, z_col]
    
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(x, y, c=z, cmap='plasma', marker='o', edgecolors='k', s=100)
    cbar = plt.colorbar(scatter)
    cbar.set_label('Decay Time (Gyr)')
    
    plt.xlabel(r'$f_{gd}$')
    plt.ylabel(r'Inclination ($\degree$)')
    plt.show()


# Plot decay time as a colormap (grid of color-coded cells)
def plot_colormap(data, x_col, y_col, z_col, fixed_params):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)
    
    filtered_data = data[mask]
    #x = np.log10(filtered_data[:, x_col])  # Apply log to mass
    x = filtered_data[:, x_col]
    y = filtered_data[:, y_col]
    z = filtered_data[:, z_col]
    
    xi = np.linspace(min(x), max(x), 50)
    yi = np.linspace(min(y), max(y), 50)
    grid_x, grid_y = np.meshgrid(xi, yi)
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')
    
    plt.figure(figsize=(10, 6))
    plt.imshow(grid_z, extent=[min(x), max(x), min(y), max(y)], origin='lower', aspect='auto', cmap='plasma_r')
    cbar = plt.colorbar()
    cbar.set_label('Decay Time (Gyr)')
    
    #plt.xlabel(r'$\log(\frac{M}{M_{\odot}})$')
    plt.xlabel(r'$f_g$')
    plt.ylabel(r'Inclination ($\degree$)')
    #plt.title('Deca')
    plt.show()

# Plot decay time as a colormap with white space between columns to reflect discrete x values

# Plot decay time as a colormap with white space between columns to reflect discrete x values

# Plot decay time as a colormap with white space between columns to reflect discrete x values
# Plot decay time as a colormap with white space between columns to reflect discrete x values
def plot_colormap_space(data, x_col, y_col, z_col, fixed_params, log_x=True, spacing=0.3):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)

    filtered_data = data[mask]
    x_data = np.log10(filtered_data[:, x_col]) if log_x else filtered_data[:, x_col]
    x_vals = np.unique(x_data)
    y_vals = np.unique(filtered_data[:, y_col])

    # Create a grid to hold z values
    z_grid = np.full((len(y_vals), len(x_vals)), np.nan)

    for i, xv in enumerate(x_vals):
        for j, yv in enumerate(y_vals):
            mask = (np.isclose(x_data, xv) & np.isclose(filtered_data[:, y_col], yv))
            if np.any(mask):
                z_grid[j, i] = np.mean(filtered_data[mask][:, z_col])

    # Define new width to simulate spacing
    bar_width = 1
    total_width = bar_width + spacing

    fig, ax = plt.subplots(figsize=(10, 6))

    for i in range(len(x_vals)):
        x_start = i * total_width
        for j in range(len(y_vals)):
            value = z_grid[j, i]
            if not np.isnan(value):
                rect = plt.Rectangle((x_start, j), bar_width, 1, color=plt.cm.plasma((value - np.nanmin(z_grid)) / (np.nanmax(z_grid) - np.nanmin(z_grid))))
                ax.add_patch(rect)

    ax.set_xlim(-0.3, len(x_vals) * total_width)
    ax.set_ylim(0, len(y_vals))
    ax.set_xticks([i * total_width + bar_width / 2 for i in range(len(x_vals))])
    ax.set_xticklabels([f"{val:.2f}" for val in x_vals])
    ax.set_yticks(range(len(y_vals)))
    ax.set_yticklabels([f"{val:.0f}" for val in y_vals])

    #ax.set_xlabel(f'{"Log10(" if log_x else ""}Parameter {x_col}{")" if log_x else ""} (Discrete)')
    #ax.set_ylabel(f'Parameter {y_col}')
    #ax.set_title('Colormap of Decay Time with Spaced Columns')
    
    ax.set_xlabel(r'$f$')
    ax.set_ylabel(r'$\dot{\iota}$ ($\degree$)')
    #ax.set_title('Colormap of Decay Time with Spaced Columns')

    sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=np.nanmin(z_grid), vmax=np.nanmax(z_grid)))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('Decay Time (Gyr)')

    plt.tight_layout()
    plt.show()
# Plot grouped bar chart of decay time
def plot_grouped_bar(data, x_col, y_col, z_col, fixed_params, log_x=True):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)

    filtered_data = data[mask]
    x_data = np.log10(filtered_data[:, x_col]) if log_x else filtered_data[:, x_col]
    y_data = filtered_data[:, y_col]
    z_data = filtered_data[:, z_col]

    x_vals = np.unique(x_data)
    y_vals = np.unique(y_data)
    width = 0.8 / len(y_vals)  # width of each bar within a group

    cmap = get_cmap('plasma')
    colors = [cmap(i / (len(y_vals) - 1)) for i in range(len(y_vals))]

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (yv, color) in enumerate(zip(y_vals, colors)):
        bar_heights = []
        for xv in x_vals:
            mask = (np.isclose(x_data, xv) & np.isclose(y_data, yv))
            if np.any(mask):
                bar_heights.append(np.mean(z_data[mask]))
            else:
                bar_heights.append(0)
        x_positions = np.arange(len(x_vals)) + i * width - 0.4 + width/2
        ax.bar(x_positions, bar_heights, width=width, label=f'{yv:.1f}', color=color)

    ax.set_xticks(np.arange(len(x_vals)))
    ax.set_xticklabels([f'{val:.2f}' for val in x_vals])
    ax.set_xlabel(r'$f_{\rm g}$')
    ax.set_ylabel('Decay Time (Gyr)')
    #ax.set_title('Grouped Bar Plot of Decay Time')
    ax.legend(title=r'$\dot{\iota}$ $(\degree)$')
    plt.tight_layout()
    plt.show()

def plot_grouped_bar_1(data, x_col, y_col, z_col, fixed_params, log_x=True):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)

    filtered_data = data[mask]
    x_data = np.log10(filtered_data[:, x_col]) if log_x else filtered_data[:, x_col]
    y_data = filtered_data[:, y_col]
    z_data = filtered_data[:, z_col]

    x_vals = np.unique(x_data)
    y_vals = np.unique(y_data)
    width = 0.8 / len(y_vals)  # width of each bar within a group

    cmap = get_cmap('plasma')
    colors = [cmap(i / (len(y_vals) - 1)) for i in range(len(y_vals))]

    fig, ax = plt.subplots(figsize=(9, 6))

    for i, (yv, color) in enumerate(zip(y_vals, colors)):
        bar_heights = []
        for xv in x_vals:
            mask = (np.isclose(x_data, xv) & np.isclose(y_data, yv))
            if np.any(mask):
                bar_heights.append(np.mean(z_data[mask]))
            else:
                bar_heights.append(0)
        x_positions = np.arange(len(x_vals)) + i * width - 0.4 + width / 2
        ax.bar(x_positions, bar_heights, width=width, label=f'{yv:.1f}', color=color)

    ax.set_xticks(np.arange(len(x_vals)))
    #ax.set_xticklabels([f'{val:.2f}' for val in x_vals], fontsize=16)
    #ax.set_xticklabels([r'$100$', r'$200$', r'$300$'], fontsize=16)
    ax.set_xticklabels([r'$10^6$', r'$10^7$', r'$10^8$'], fontsize=16)
    #ax.set_xticklabels([r'$0.3$', r'$0.5$', r'$0.7$'], fontsize=16)
    #ax.set_xlabel(r'$n_{\rm gd}$ (cm$^{-3}$)', fontsize=16)
    #ax.set_xlabel(r'$f_{\rm g}$', fontsize=16)
    ax.set_xlabel(r'$M_{\rm 1} \ (M_{\odot})$', fontsize=16)
    
    #ax.set_ylabel(r'$t_{\rm p}$ (Gyr)', fontsize=16)

    # Set tick label font sizes
    ax.tick_params(axis='both', labelsize=16)

    # Move and format legend
    ax.legend(title=r'$\dot{\iota}_0$ (deg)', loc='center left', bbox_to_anchor=(1.02, 0.5),
              title_fontsize=16, fontsize=16)

    # Adjust layout to accommodate legend
    plt.tight_layout()
    plt.subplots_adjust(right=0.75)
    
    plt.show()

def plot_grouped_bar_vg(data, x_col, y_col, z_col, fixed_params, log_x=True):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)

    filtered_data = data[mask]
    x_data = np.log10(filtered_data[:, x_col]) if log_x else filtered_data[:, x_col]
    y_data = filtered_data[:, y_col]
    z_data = filtered_data[:, z_col]

    # Sort x_vals in decreasing order instead of increasing
    x_vals = np.unique(x_data)[::-1]  # <-- reverse the array here
    y_vals = np.unique(y_data)
    width = 0.8 / len(y_vals)

    cmap = get_cmap('plasma')
    colors = [cmap(i / (len(y_vals) - 1)) for i in range(len(y_vals))]

    fig, ax = plt.subplots(figsize=(8, 6))

    for i, (yv, color) in enumerate(zip(y_vals, colors)):
        bar_heights = []
        for xv in x_vals:
            mask = (np.isclose(x_data, xv) & np.isclose(y_data, yv))
            if np.any(mask):
                bar_heights.append(np.mean(z_data[mask]))
            else:
                bar_heights.append(0)
        x_positions = np.arange(len(x_vals)) + i * width - 0.4 + width / 2
        ax.bar(x_positions, bar_heights, width=width, label=f'{yv:.1f}', color=color)

    ax.set_xticks(np.arange(len(x_vals)))

    # Set custom tick labels (you manually input them here)
    ax.set_xticklabels([r'$0.7$', r'$0.5$', r'$0.3$'][::-1], fontsize=16)  # <-- reverse manual labels too

    ax.set_xlabel(r'$v_{\rm g} / v_{\rm c}(r,z)$', fontsize=16)
    ax.set_ylabel(r'$t_p$ (Gyr)', fontsize=16)
    ax.tick_params(axis='both', labelsize=16)

    #ax.legend(title=r'$\dot{\iota}_{\rm 0}$ (deg)', loc='center left', bbox_to_anchor=(1.02, 0.5),
    #          title_fontsize=20, fontsize=20)

    plt.tight_layout()
    plt.subplots_adjust(right=0.95)

    plt.show()

# Plot decay time as a tile plot with labeled values
def plot_tile_labeled(data, x_col, y_col, z_col, fixed_params, log_x=True):
    mask = np.ones(data.shape[0], dtype=bool)
    for col, val in fixed_params.items():
        mask &= np.isclose(data[:, col], val)

    filtered_data = data[mask]
    x_data = np.log10(filtered_data[:, x_col]) if log_x else filtered_data[:, x_col]
    y_data = filtered_data[:, y_col]
    z_data = filtered_data[:, z_col]

    x_vals = np.unique(x_data)
    y_vals = np.unique(y_data)
    z_grid = np.full((len(y_vals), len(x_vals)), np.nan)

    for i, xv in enumerate(x_vals):
        for j, yv in enumerate(y_vals):
            mask = (np.isclose(x_data, xv) & np.isclose(y_data, yv))
            if np.any(mask):
                z_grid[j, i] = np.mean(z_data[mask])

    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.viridis_r
    im = ax.imshow(z_grid, cmap=cmap, origin='lower')

    for i in range(len(y_vals)):
        for j in range(len(x_vals)):
            val = z_grid[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.1f}", ha='center', va='center', color='black', fontsize=8)

    ax.set_xticks(np.arange(len(x_vals)))
    ax.set_xticklabels([f"{x:.2f}" for x in x_vals])
    ax.set_yticks(np.arange(len(y_vals)))
    ax.set_yticklabels([f"{y:.0f}" for y in y_vals])

    #ax.set_xlabel(f'{"Log10(" if log_x else ""}Parameter {x_col}{")" if log_x else ""}')
    ax.set_xlabel(r'$f_g$')
    
    ax.set_ylabel(r'$\dot{\iota}$ ($\degree$)')
    #ax.set_title('Tile Plot of Decay Time with Labels')

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Decay Time (Gyr)')
    plt.tight_layout()
    plt.show()

def plot_decay_vs_angle(data, columns):
    angles = data[:, columns['angle']]
    decay_times = data[:, columns['decay_time']]

    plt.figure(figsize=(8, 5))
    plt.scatter(angles, decay_times, color='teal', alpha=0.7, edgecolor='k', s=30)

    plt.xlabel('Angle (degrees or radians)', fontsize=14)
    plt.ylabel('Decay Time (Gyr)', fontsize=14)
    plt.title('Decay Time vs. Angle', fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_decay_vs_angle_filtered(data, columns, fixed_param='fg', fixed_value=0.5, tolerance=1e-2):
    """
    Plots decay time vs. angle for a slice of data where the fixed_param is approximately fixed_value.
    
    Parameters:
    - data: ndarray with your full dataset
    - columns: dict mapping column names to indices
    - fixed_param: str, one of the keys in columns (e.g., 'velocity')
    - fixed_value: float, value of the parameter to fix
    - tolerance: float, allowable deviation when filtering the parameter
    """
    col_index = columns[fixed_param]
    mask = np.abs(data[:, col_index] - fixed_value) < tolerance
    filtered_data = data[mask]

    angles = filtered_data[:, columns['angle']]
    decay_times = filtered_data[:, columns['decay_time']]

    plt.figure(figsize=(8, 5))
    plt.scatter(angles, decay_times, color='teal', alpha=0.7, edgecolor='k', s=30)

    plt.xlabel('Angle (degrees or radians)', fontsize=14)
    plt.ylabel('Decay Time (Gyr)', fontsize=14)
    plt.title(f'Decay Time vs. Angle (Fixed {fixed_param} = {fixed_value})', fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
   
'''
def format_latex_coeff(c, power):
    """Format a polynomial coefficient into LaTeX-style with scientific notation."""
    if abs(c) < 1e-2 or abs(c) > 1e2:
        base = f"{c:.1e}"
        coeff, exp = base.split('e')
        exp = int(exp)
        coeff_str = f"{float(coeff):.1f}\\times10^{{{exp}}}"
    else:
        coeff_str = f"{c:.2f}"

    if power == 0:
        return f"{coeff_str}"
    elif power == 1:
        return f"{coeff_str}x"
    else:
        return f"{coeff_str}x^{{{power}}}"
'''
def format_latex_coeff(c, power):
    """Format a polynomial coefficient into LaTeX-style with scientific notation, using theta instead of x."""
    if abs(c) < 1e-2 or abs(c) > 1e2:
        base = f"{c:.1e}"
        coeff, exp = base.split('e')
        exp = int(exp)
        coeff_str = f"{float(coeff):.1f}\\times10^{{{exp}}}"
    else:
        coeff_str = f"{c:.2f}"

    if power == 0:
        return f"{coeff_str}"
    elif power == 1:
        return f"{coeff_str}i"
    else:
        return f"{coeff_str}i^{{{power}}}"
def plot_decay_vs_angle_filtered_bands(data, columns, fixed_param='fg', fixed_value=0.5, tolerance=1e-2):
    """
    Plots decay time vs. angle for a slice of data where the fixed_param is approximately fixed_value,
    fits a 3rd order polynomial, shows ±1σ and ±2σ bands, and displays the polynomial equation.
    """
    col_index = columns[fixed_param]
    mask = np.abs(data[:, col_index] - fixed_value) < tolerance
    filtered_data = data[mask]

    angles = filtered_data[:, columns['angle']]
    decay_times = filtered_data[:, columns['decay_time']]

    # Fit a 3rd-order polynomial
    coeffs = np.polyfit(angles, decay_times, deg=3)
    poly = np.poly1d(coeffs)

    # Smooth curve
    angle_range = np.linspace(min(angles), max(angles), 300)
    fit_vals = poly(angle_range)

    # Standard deviation of residuals
    residuals = decay_times - poly(angles)
    std_dev = np.std(residuals)

    # Format polynomial equation in LaTeX
    terms = [format_latex_coeff(c, p) for p, c in zip(range(3, -1, -1), coeffs)]
    equation_str = r"$\tau_{decay} = " + " + ".join(terms).replace("+-", "- ") + r"$"

    # Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(angles, decay_times, color='blue', alpha=0.7, edgecolor='k', s=30, label='Data')
    plt.plot(angle_range, fit_vals, color='darkred', linewidth=2.5, label='3rd Order Fit')

    # Bollinger-style bands
    for n, color, alpha in zip([1, 2], ['orange', 'gold'], [0.3, 0.2]):
        plt.fill_between(angle_range, fit_vals - n * std_dev, fit_vals + n * std_dev,
                         color=color, alpha=alpha, label=f'±{n}σ')

    # Labels and title
    plt.xlabel(r'$\dot{\iota}$ ($\degree$)', fontsize=16)
    plt.ylabel('Decay Time (Gyr)', fontsize=16)
    #plt.title(f'Decay Time vs. Angle (Fixed {fixed_param} = {fixed_value})', fontsize=18)
    plt.ylim([0, 14.5])

    # Display equation
    plt.text(0.95, 0.05, equation_str, transform=plt.gca().transAxes,
             fontsize=13, verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(facecolor='white', alpha=0.7))

    plt.grid(True)
    plt.tick_params(labelsize = 16)
    #plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()



def kernel_smooth_std(x, residuals, x_eval, bandwidth):
    """
    Estimate local standard deviation using Gaussian kernel smoothing.
    """
    stds = []
    for x0 in x_eval:
        weights = np.exp(-0.5 * ((x - x0) / bandwidth)**2)
        weights /= np.sum(weights)
        mean_res = np.sum(weights * residuals)
        variance = np.sum(weights * (residuals - mean_res)**2)
        stds.append(np.sqrt(variance))
    return np.array(stds)

def plot_decay_vs_angle_kernel_bands(data, columns, fixed_param='fg', fixed_value=0.5, tolerance=1e-2, bandwidth=5.0):
    """
    Plots decay time vs. angle with a polynomial fit and adaptive ±1σ, ±2σ bands using kernel-smoothed residual std.
    """
    col_index = columns[fixed_param]
    mask = np.abs(data[:, col_index] - fixed_value) < tolerance
    filtered_data = data[mask]

    angles = filtered_data[:, columns['angle']]
    decay_times = filtered_data[:, columns['decay_time']]

    # Fit a 3rd-order polynomial
    coeffs = np.polyfit(angles, decay_times, deg=3)
    poly = np.poly1d(coeffs)

    # Evaluate fit
    angle_range = np.linspace(min(angles), max(angles), 300)
    fit_vals = poly(angle_range)

    # Residuals for smoothing
    residuals = decay_times - poly(angles)

    # Kernel-smoothed standard deviation of residuals
    local_std = kernel_smooth_std(angles, residuals, angle_range, bandwidth)

    # Format polynomial equation
    def format_latex_coeff(c, p):
        if abs(c) < 1e-4:
            c_str = f"{c:.5e}"
        else:
            c_str = f"{c:.5e}"
        if p == 0:
            return c_str
        elif p == 1:
            return f"{c_str}i"
        else:
            return f"{c_str}i^{p}"
    terms = [format_latex_coeff(c, p) for p, c in zip(range(3, -1, -1), coeffs)]
    #equation_str = r"$\tau_{decay} = " + " + ".join(terms).replace("+-", "- ") + r"$"
    #print(equation_str)
    # Plot
    plt.figure(figsize=(9, 5))
    plt.scatter(angles, decay_times, color='blue', alpha=0.7, edgecolor='k', s=30)
    plt.plot(angle_range, fit_vals, color='darkred', linewidth=2.5)

    # ±1σ and ±2σ bands
    for n, color, alpha in zip([1, 2], ['orange', 'gold'], [0.3, 0.2]):
        upper = fit_vals + n * local_std
        lower = fit_vals - n * local_std
        plt.fill_between(angle_range, lower, upper, color=color, alpha=alpha)

    plt.xlabel(r'$i_{\rm 0}$ (deg)', fontsize=16)
    plt.ylabel(r'$t_{\rm p}$ (Gyr)', fontsize=16)
    plt.ylim([0, 14.7])
    
    #plt.text(0.95, 0.05, equation_str, transform=plt.gca().transAxes,
    #         fontsize=13, verticalalignment='bottom', horizontalalignment='right',
    #         bbox=dict(facecolor='white', alpha=0.7))

    plt.grid(True)
    plt.tick_params(labelsize=16)
    plt.tight_layout()
    plt.savefig('/Users/sghob/OneDrive/Desktop/Res_code/108.png', format='png', bbox_inches='tight')
    plt.show()
    

def sigmoid(x, A, B, C):
    """Sigmoid function: A + B / (1 + exp(-C x))."""
    return A + B / (1.0 + np.exp(-C * x))

def plot_density_kde(data, columns, fixed_param='velocity', fixed_value=300.0, tolerance=1e-2, levels=10):
    """
    Plots decay time vs. angle as a KDE contour plot.
    """
    col_index = columns[fixed_param]
    mask = np.abs(data[:, col_index] - fixed_value) < tolerance
    filtered_data = data[mask]

    angles = filtered_data[:, columns['angle']]
    decay_times = filtered_data[:, columns['decay_time']]

    plt.figure(figsize=(8, 5))
    sns.kdeplot(
        x=angles, y=decay_times,
        fill=True, cmap='magma', levels=levels, thresh=0
    )
    plt.xlim([0,85.0])
    plt.xlabel('Angle (degrees or radians)', fontsize=14)
    plt.ylabel('Decay Time (Gyr)', fontsize=14)
    plt.title(f'Density Contour (Fixed {fixed_param} = {fixed_value})', fontsize=16)
    plt.tight_layout()
    plt.show()

def plot_density_histogram(data, columns, fixed_param='velocity', fixed_value=300.0, tolerance=1e-2, bins=10):
    """
    Plots decay time vs. angle as a 2D histogram density plot.
    """
    col_index = columns[fixed_param]
    mask = np.abs(data[:, col_index] - fixed_value) < tolerance
    filtered_data = data[mask]

    angles = filtered_data[:, columns['angle']]
    decay_times = filtered_data[:, columns['decay_time']]

    plt.figure(figsize=(8, 5))
    plt.hist2d(angles, decay_times, bins=bins, cmap='plasma')
    plt.colorbar(label='Counts')
    plt.xlabel('Angle (degrees or radians)', fontsize=14)
    plt.ylabel('Decay Time (Gyr)', fontsize=14)
    plt.title(f'Density of Points (Fixed {fixed_param} = {fixed_value})', fontsize=16)
    plt.tight_layout()
    plt.show()

def plot_decay_vs_angle_given_fixed_params(data, columns, fixed_params):
    """
    Plots decay time vs. angle for rows in 'data' that match fixed parameter values.

    Parameters:
    - data: np.ndarray, the full dataset
    - columns: dict, mapping of parameter names to column indices
    - fixed_params: dict, can use either parameter names or column indices as keys
    """
    mask = np.ones(len(data), dtype=bool)

    for param_key, value in fixed_params.items():
        # If key is a string, use the columns dictionary to get index
        if isinstance(param_key, str):
            if param_key not in columns:
                raise KeyError(f"'{param_key}' is not a valid key in the columns dictionary.")
            col_index = columns[param_key]
        elif isinstance(param_key, int):
            col_index = param_key
        else:
            raise ValueError(f"Parameter key must be a string or int, got {type(param_key)}")

        mask &= (data[:, col_index] == value)

    filtered_data = data[mask]
    if len(filtered_data) == 0:
        print("No data matches the specified fixed parameters.")
        return

    angle_col = columns['angle']
    decay_col = columns['decay_time']

    angle = filtered_data[:, angle_col]
    decay_time = filtered_data[:, decay_col]

    plt.figure(figsize=(8, 6))
    plt.scatter(angle, decay_time, color='dodgerblue', edgecolor='k', s=60)
    plt.xlabel(r'$\dot{\iota}$ ($\degree$)', fontsize=18)
    plt.ylabel('Decay Time (Gyr)', fontsize=18)
    plt.tick_params(axis='both', labelsize = 18)
    #plt.title('Decay Time vs. Angle (Fixed Parameters)', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_three_column_data_dual_y(file_path):
    """
    Reads a tab-separated 3-column .txt file and creates a scatter plot with dual y-axes:
    - Column 1 is the x-axis
    - Column 2 is plotted on the left y-axis
    - Column 3 is plotted on the right y-axis
    """
    # Load the data
    data = np.loadtxt(file_path, delimiter='\t')

    x = data[:, 0]
    y1 = data[:, 1] / 3.15E16  # Column 2 (normalized)
    y2 = data[:, 2]            # Column 3

    fig, ax1 = plt.subplots(figsize=(8, 7))

    # Plot y1 on left y-axis
    ax1.set_xlabel(r'$\dot{\iota}_0$ (deg)', fontsize=18)
    ax1.set_ylabel(r'$t_{\rm p}$ (Gyr)', fontsize=18, color='black')
    scatter1 = ax1.scatter(x, y1, s = 150, color='purple', label=r'$t_{\rm p}$')
    ax1.tick_params(axis='both', labelsize=18, labelcolor='black')
    ax1.grid(True)

    # Plot y2 on right y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$\langle \mathcal{M}_{\rm z} \rangle$', fontsize=18, color='black')
    scatter2 = ax2.scatter(x, y2, s = 150, color='green', marker = 'd', label=r'$\langle \mathcal{M}_{\rm z} \rangle$')
    ax2.tick_params(axis='y', labelsize=18, labelcolor='black')

    # Combine legends from both axes
    handles = [scatter1, scatter2]
    labels = [h.get_label() for h in handles]
    ax1.legend(handles, labels, fontsize=18)

    # Layout
    fig.tight_layout()
    plt.savefig('/Users/sghob/OneDrive/Desktop/Res_code/twinscatter.png', format='png', bbox_inches='tight')
    plt.show()
    
here = os.path.dirname(os.path.abspath(__file__))
filename = os.path.join(here, 'decay_time_params_frfrfr.txt')
#filename = 'your_data.txt'  # Replace with your actual file path
data = load_data(filename)

columns = {
    'mass': 0,
    'fg': 1,
    'vg': 2,
    'ngd': 3,
    'angle': 4,
    'decay_time': 5,
    'flag': 6
}


x_param = columns['mass']
y_param = columns['angle']
z_param = columns['decay_time']

fixed_parameters = {
    columns['ngd']: 100,
    columns['fg']: 0.7,
    #columns['mass']: 1E6,
    columns['vg']: -0.7
}

fixed_parameters1 = {
    columns['ngd']: 100,
    columns['vg']: -0.7,
    columns['mass']: 1E6,
    columns['fg']: 0.7
}


mass_choice = 1E8
vg_omit = -0.7

mass_col = columns['mass']
vg_col = columns['vg']

mass_mask = np.isclose(data[:, mass_col], mass_choice)
vg_mask = np.isclose(data[:, vg_col], vg_omit)

subset_data = data[mass_mask & vg_mask]


#plot_3d_scatter(data, x_param, y_param, z_param, fixed_parameters)
#plot_2d_discrete_heatmap(data, x_param, y_param, z_param, fixed_parameters)
#plot_binned_scatter(data, x_param, y_param, z_param, fixed_parameters)
#plot_colormap(data, x_param, y_param, z_param, fixed_parameters)
#plot_colormap_space(data, x_param, y_param, z_param, fixed_parameters, log_x=False, spacing = 0.3)
#plot_grouped_bar_1(data, x_param, y_param, z_param, fixed_parameters, log_x=False)
#plot_grouped_bar_vg(data, x_param, y_param, z_param, fixed_parameters, log_x=False)
#plot_tile_labeled(data, x_param, y_param, z_param, fixed_parameters, log_x=False)
#plot_decay_vs_angle(data, columns)
#plot_decay_vs_angle_given_fixed_params(data, columns, fixed_parameters1)
#plot_decay_vs_angle_filtered_bands(data, columns, fixed_param='mass', fixed_value=1E8)
#plot_density_kde(data, columns, fixed_param='mass', fixed_value=1E8)
#plot_density_histogram(data, columns, fixed_param='fg', fixed_value=0.5)
#plot_decay_vs_angle_adaptive_bands(data, columns, fixed_param='mass', fixed_value=1E7)
#plot_decay_vs_angle_kernel_bands(data, columns, fixed_param='mass', fixed_value=1E8)
plot_decay_vs_angle_kernel_bands_sig(subset_data, columns, fixed_param='mass', fixed_value = 1E8)
#plot_three_column_data_dual_y('/Users/sghob/OneDrive/Desktop/Res_code/time_machavg.txt')
#feature_columns = [columns['mass'], columns['fg'], columns['vg'], columns['ngd'], columns['angle']]
#target_column = columns['decay_time']

#model, scaler = train_gp_model(data, feature_columns, target_column)
#plot_3d_prediction_surface(data, model, scaler, columns)