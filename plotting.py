import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import os
from scipy.signal import find_peaks

def average_over_periods(t, y):
    """
    Detect zero-crossings to define oscillation periods and 
    compute the average value of y within each period.
    
    Parameters:
        t, y: lists or arrays of equal length
    
    Returns:
        t_avg, y_avg: lists of averaged time and observable values
    """
    
    # Convert to numpy arrays
    t = np.asarray(t)
    y = np.asarray(y)
    
    # Find zero-crossings: where y changes sign
    sign_change_indices = np.where(np.diff(np.sign(y)) != 0)[0]
    
    t_avg = []
    y_avg = []
    
    # Loop over oscillation periods defined by zero-crossings
    for i in range(len(sign_change_indices) - 1):
        start = sign_change_indices[i]
        end = sign_change_indices[i + 1]
        
        # Slice data for this oscillation period
        period_t = t[start:end+1]
        period_y = y[start:end+1]
        
        # Compute averages
        t_avg.append(0.5 * (period_t[0] + period_t[-1]))  # midpoint time
        y_avg.append(np.mean(period_y))
    
    return t_avg, y_avg

here = os.path.dirname(os.path.abspath(__file__))
here1 = os.path.dirname(os.path.abspath(__file__))
here2 = os.path.dirname(os.path.abspath(__file__))
#filename = os.path.join(here, 'energy.txt')
filename = os.path.join(here, 'mass.txt')
#filename = os.path.join(here, "l_lat.txt")
#filename = os.path.join(here, 'decay_time.txt')
filename5 = os.path.join(here, 'pbz125.txt')
#filename5 = os.path.join(here, 'firstfew45.txt')

filename10 = os.path.join(here, 'pr.txt')
filename20 = os.path.join(here, 'theta20.txt')

filename1 = os.path.join(here, 'compress.txt')
filename2 = os.path.join(here, 'pot_fig.csv')
filename3 = os.path.join(here, 'phi_load200.csv')
#filename = os.path.join(here, 'mach.txt')

################### 3D PLOTTING
dats3d = np.genfromtxt(filename3, delimiter=',')
#log_norm = matplotlib.colors.LogNorm(vmin=-1, vmax=1)
x_start = -1.0
x_end = 1.0
y_start = -1.0
y_end = 1.0
#plt.xlabel(r'$r$ ($R_{g}$)')
#plt.ylabel(r'$z$ ($R_{g}$)')
#plt.title('Total Potential')
#plt.gca().invert_yaxis()

#'''
#im = plt.imshow(1.0 * dats3d, cmap='magma_r', interpolation='nearest', extent=[x_start, x_end, y_start, y_end], origin='lower', vmin = -0.4E-10, vmax = 0.4E-10)
#im = plt.imshow(1.0 * dats3d, cmap='magma_r', interpolation='nearest', extent=[x_start, x_end, y_start, y_end], origin='lower')

#plt.title('')

plt.xlabel(r'Radius ($r/R_{\rm g}$)', fontsize = 16)
plt.ylabel(r'Vertical height ($z/0.1R_{\rm g}$)', fontsize = 16)# ($z/R_g$)',)
plt.tick_params(axis='both', labelsize = 16)
# cbar = plt.colorbar(im)
# cbar.ax.tick_params(labelsize = 16)
# cbar.ax.yaxis.get_offset_text().set_fontsize(16)  # or whatever fontsize you want

#cbar.set_label(r'$\partial_z\Phi(r,z)$ (m$^2$ / s$^2$)', fontsize=16)
#cbar.set_label(r'$B(r,z)$ (G)', fontsize=16)

#plt.colorbar(label=r'$\Phi(r,z)$ (m$^2$ / s$^2$)')
plt.savefig("figure1.pdf")
plt.show()
#print(os.getcwd())
#'''

# read the data from the txt file
with open(filename, 'r') as f:
    data = f.readlines()
    
with open(filename5, 'r') as f:
    data5 = f.readlines()
with open(filename10, 'r') as f:
    data10 = f.readlines()
with open(filename20, 'r') as f:
    data20 = f.readlines()

# extract x and y values from data
x = []
y = []

xavg = []
yavg = []

x5 = []
y5 = []

x10 = []
y10 = []

x20 = []
y20 = []

yr_s = 3.15E7
pc_m = 3.086E16
#for time evolution
#'''
for line in data:
    values = line.split()
    datx = float(values[0]) / (3.15E7 * 1E9) 
   # datx = float(values[0]) / (3.15E7 * 1E9) 
    daty = (float(values[1]))# / pc_m
    #daty = (float(values[1]))# / (3.15E7 * 1E9) 
    x.append(datx)#1.59E-23
    y.append(daty)
#'''
def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

# Choose an appropriate window size (number of data points to average over)
window_size = 1000  # Adjust this based on how smooth you want the result

average_radius0 = moving_average(y, window_size)

# Adjust the time array for the moving average (since it is shorter)
average_time0 = x[(window_size - 1):]

for line in data5:
    values = line.split()
    datx1 = float(values[0])/ (3.15E7 * 1E9)
    daty1 = float(values[1])# / pc_m
    x5.append(datx1)
    y5.append(daty1)

average_radius5 = moving_average(y5, window_size)
average_time5 = x5[(window_size - 1):]

for line in data10:
    values = line.split()
    datx2 = float(values[0]) / (3.15E7 * 1E9)
    daty2 = float(values[1])# / pc_m
    x10.append(datx2)
    y10.append(daty2)

average_radius10 = moving_average(y10, window_size)
average_time10 = x10[(window_size - 1):]

for line in data20:
    values = line.split()
    datx3 = float(values[0]) / (3.15E7 * 1E9)
    daty3 = float(values[1]) / pc_m
    x20.append(datx3)
    y20.append(daty3)

average_radius20 = moving_average(y20, window_size)
average_time20 = x20[(window_size - 1):]

#for purely spatial data
'''
for line in data:
    values = line.split()
    datx = float(values[0])
    daty = float(values[1]) / pc_m
    x.append(daty * np.cos(datx))   
    y.append(daty * np.sin(datx)) 
'''
xsamp = np.linspace(0,300,10000)
y1l = 3 * np.exp(-(2/3)*0.01*xsamp)
y1s = 3 * np.exp(-(1/2)*0.01*xsamp)

xavg, yavg = average_over_periods(x, y)
'''
x = np.asarray(x)
y = np.asarray(y)

peaks,_ = find_peaks(y, distance=500)
time_peaks = x[peaks]
mach_peaks = y[peaks]
'''
'''

plt.plot(xsamp, y1s, color = 'orange', linestyle = 'dashed')
plt.plot(xsamp, y1l, color = 'red', linestyle = 'dashed')

plt.xlabel('t')
plt.ylabel(r'$z/z_g$')
plt.legend([r'$z$-position', r'Small $\varepsilon$ amplitude', r'Large $\varepsilon$ amplitude'])
plt.title(r'$z$-evolution')
'''
accline = np.full_like(x, 2.48E-5 / (2.2E-2 / 9))
plt.figure(figsize=(7,7))
#plt.plot(x, y, color = 'red', linestyle = 'solid')
plt.semilogy(x, y, color = 'orange', linestyle = 'solid')
#plt.plot(x, accline)
#plt.plot(time_peaks, mach_peaks, color = 'purple', linestyle = 'solid')

plt.grid()
#plt.xlabel(r'$\theta$  $(\degree)$')
#plt.ylabel(r'$z_{amp}$ (pc)')
plt.xlabel(r'$t$ (Gyr)', fontsize = 16)
#plt.xlabel(r'$r$ (pc)', fontsize = 16)
#plt.ylabel(r'$B_{\infty}$ (G)', fontsize = 16)
#plt.ylabel(r'$z$ (pc)', fontsize = 30)
#plt.ylabel(r'$\mathcal{M}$', fontsize = 14)
#plt.ylabel(r'$\dot{\iota} \ (^{\circ})$', fontsize = 16)
#plt.ylabel(r'$\dot{\iota} \ (^{\circ})$', fontsize = 16)
#plt.ylabel(r'$\dot{M} / \dot{M}_{E}$', fontsize = 16)
plt.ylabel(r'$P_{BZ}$ (erg/s)', fontsize = 16)
#plt.ylabel(r'$M$ ($M_{\odot}$)', fontsize = 16)
#plt.ylabel(r'$\varepsilon_{jet}$', fontsize = 16)
#plt.title(r'Angular Momentum Evolution at 75 $\degree$ Inclination')
#plt.title(r'Change in decay time as a function of inclination')
plt.tick_params(labelsize = 16)
#plt.ticklabel_format(axis='both', style='scientific', scilimits=(0, 0))
#plt.plot(x5, y5, color = 'red', linestyle = 'dashed')
#plt.plot(x10, y10, color = 'red', linestyle = 'dashed')
#plt.legend(['Secondary', 'Primary',], fontsize = 16)

'''
plt.plot(average_time0, average_radius0, color = 'red', linestyle = 'solid')

plt.plot(average_time10, average_radius10, color = 'orange', linestyle = 'solid')
plt.plot(average_time20, average_radius20, color = 'purple', linestyle = 'solid')
plt.plot(average_time5, average_radius5, color = 'blue', linestyle = 'solid')
plt.xlabel('t (Gyr)')
plt.ylabel('r (pc)')
plt.title(r'Radial evolution of secondary MBHs with varying orbital inclinations')
plt.legend([r'$\theta = 0\degree$', r'$\theta = 10\degree$', r'$\theta = 20\degree$', r'$\theta = 30\degree$'])
'''
'''
plt.xlabel('Time (Gyr)')
plt.ylabel('Relative Energy Error')

'''
plt.savefig('/Users/sghob/OneDrive/Desktop/Res_code/rz45.png', format='png', bbox_inches='tight')
plt.show()
