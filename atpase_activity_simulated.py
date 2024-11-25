import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

# Plot settings
plt.rcParams["font.family"] = "Garuda"
plt.rc('font', size=14)       # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=18)  # fontsize of the figure title
#params = {'mathtext.default': 'regular' }   # não deixa mais o texto no mathmode em itálico
#params = {'mathtext.default': 'it' }   # deixa o texto do mathmode em itálico
#plt.rcParams.update(params)   # Se fizer alguma das alterações acima, rodar este comando para atualizar
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors = [cycle[0], cycle[1]]



# - - - - METHODS - - - - #

# Method to generate simulated data
def simulate_data(x, avg, sigma):
    sim_data = np.zeros(len(x))
    for i in range(len(x)):
        if avg[i] < 0: sim_data[i] = 0
        else: sim_data[i] = np.random.normal(loc = avg[i], scale = sigma[i])
    return sim_data

# Select range at which limit the range (from start to index where y <= threshold)
def range_selector(x, y, start, threshold):
    start = 0
    end = 0
    for i in range(len(x)):
        if y[i] > 0: end += 1
        else: break
    if start == end: 
        print("Can't fit with just one point!")
        SystemExit()
    else: pass
    return start, end

# Method to perform linear fits by linear least squares
def linear_fit(x, y, yerr):
    d11 = d21 = m11 = m12 = m21 = m22 = 0
    #função da reta aqui: a + bx, ou seja, p[0][0] é coef linear e p[1][0] é coef angular
    for i in range (0, len(x)):
        d11 += y[i] / (yerr[i]**2)
        d21 += (x[i]*y[i]) / (yerr[i]**2)
        m11 += 1 / (yerr[i]**2)
        m12 += x[i] / (yerr[i]**2)
        m21 = m12
        m22 += (x[i]**2) / (yerr[i]**2)
    d = np.array([[d11], [d21]])
    m = np.array([[m11,m12],[m21,m22]])
    v = np.linalg.inv(m) #matriz de variância-covariância
    p = np.matmul(v,d) #vetor de estimativas para os coeficientes
    linear = p[0][0]
    angular = p[1][0]
    linear_std = v[0][0]
    angular_std = v[1][1]
    return angular, linear, angular_std, linear_std



# - - - - SIMULATIONS - - - - #

# Simulation parameters
c_1 = 1 #uM
c_2 = 2 #uM
volume = 150 #uL
epsilon_nadh_340nm = 1498.33 #M-1 cm-1
time = np.linspace(1e-1, 3600, num = 180) # 1h of experiment

# Initializing dataframes
data_1 = pd.DataFrame({'Time_s': time})
analysis_1 = {}
data_2 = pd.DataFrame({'Time_s': time})
analysis_2 = {}

# Protein features: "real" ATPase activities
ang_1 = -7e-4
lin_1 = 0.7
sigma_1 = 0.04*np.ones(len(time))
ang_2 = -4.0e-4
lin_2 = 0.7
sigma_2 = 0.02*np.ones(len(time))

# Generating three experiments for each of the two conditions (I will ignore background subtraction to make this simpler)
for j in [0,1,2]:
    exp_vals_1 = simulate_data(time, ang_1*time+lin_1, sigma_1)
    err_vals_1 = sigma_1 * np.random.normal(loc = 1, scale = 0.05, size=len(exp_vals_1))
    exp_vals_2 = simulate_data(time, ang_2*time+lin_2, sigma_2)
    err_vals_2 = sigma_2 * np.random.normal(loc = 1, scale = 0.05, size=len(exp_vals_2))
    key_uv = 'UV_' + str(j+1)
    key_err = 'Err_' + str(j+1)
    new_data_1 = {key_uv: exp_vals_1,
                   key_err: err_vals_1}
    data_1 = data_1.assign(**new_data_1)
    new_data_2 = {key_uv: exp_vals_2,
                   key_err: err_vals_2}
    data_2 = data_2.assign(**new_data_2)
    # There, I have my simulated data! Now we define the range where we will perform the fit (before NADH is depleted)
    index_fit_1 = range_selector(data_1['Time_s'], data_1[key_uv], 0, 0)
    data_to_fit_1 = [data_1['Time_s'][index_fit_1[0]:index_fit_1[1]], \
                     data_1[key_uv][index_fit_1[0]:index_fit_1[1]], \
                     data_1[key_err][index_fit_1[0]:index_fit_1[1]]]
    index_fit_2 = range_selector(data_2['Time_s'], data_2[key_uv], 0, 0)
    data_to_fit_2 = [data_2['Time_s'][index_fit_2[0]:index_fit_2[1]], \
                     data_2[key_uv][index_fit_2[0]:index_fit_2[1]], \
                     data_2[key_err][index_fit_2[0]:index_fit_2[1]]]
    # And now we perform the fit
    ang_1, lin_1, ang_std_1, lin_std_1 = linear_fit(data_to_fit_1[0], data_to_fit_1[1], data_to_fit_1[2])
    fit_1 = ang_1 * data_1['Time_s'] + lin_1
    fit_1_dict = {'Fit_' + str(j+1): fit_1}
    data_1 = data_1.assign(**fit_1_dict)
    ativ_1 = -ang_1 / epsilon_nadh_340nm * volume * 1e6 * 60 / (volume * c_1)
    ativ_std_1 = ang_std_1 / epsilon_nadh_340nm * volume * 1e6 * 60 / (volume * c_1)
    analysis_1.update({'Activity_' + str(j+1): ativ_1, 
                        'Activity_std_' + str(j+1): ativ_std_1,
                        'Fit_start_' + str(j+1): index_fit_1[0], 'Fit_end_' + str(j+1): index_fit_1[1]})
    ang_2, lin_2, ang_std_2, lin_std_2 = linear_fit(data_to_fit_2[0], data_to_fit_2[1], data_to_fit_2[2])
    fit_2 = ang_2 * data_2['Time_s'] + lin_2
    fit_2_dict = {'Fit_' + str(j+1): fit_2}
    data_2 = data_2.assign(**fit_2_dict)
    ativ_2 = -ang_2 / epsilon_nadh_340nm * volume * 1e6 * 60 / (volume * c_2)
    ativ_std_2 = ang_std_2 / epsilon_nadh_340nm * volume * 1e6 * 60 / (volume * c_2)    
    analysis_2.update({'Activity_' + str(j+1): ativ_2, 
                        'Activity_std_' + str(j+1): ativ_std_2,
                        'Fit_start_' + str(j+1): index_fit_2[0], 'Fit_end_' + str(j+1): index_fit_2[1]})

# Calculating averages and standard deviation of extracted activities
bar_labels = np.arange(2)
atpase_1_avg = np.average([analysis_1['Activity_1'],analysis_1['Activity_2'],analysis_1['Activity_3']])
atpase_1_std = np.std([analysis_1['Activity_1'],analysis_1['Activity_2'],analysis_1['Activity_3']], ddof=1)
atpase_2_avg = np.average([analysis_2['Activity_1'],analysis_2['Activity_2'],analysis_2['Activity_3']])
atpase_2_std = np.std([analysis_2['Activity_1'],analysis_2['Activity_2'],analysis_2['Activity_3']], ddof=1)
atpase_avg = [atpase_1_avg, atpase_2_avg]
atpase_std = [atpase_1_std, atpase_2_std]



# - - - - PLOT - - - - #

# Plot the results: one of the replicates (Time vs Abs), and ATPase activity on inset
fig,ax = plt.subplots()
ax.scatter(data_1['Time_s'],data_1['UV_1'], color = cycle[0], label = 'Sample 1')
ax.errorbar(data_1['Time_s'],data_1['UV_1'],data_1['Err_1'],fmt ='none',color = cycle[0])
ax.plot(data_1['Time_s'][analysis_1['Fit_start_1']: analysis_1['Fit_end_1']], \
        data_1['Fit_1'][analysis_1['Fit_start_1']: analysis_1['Fit_end_1']], linewidth = 2, color = 'k')
ax.scatter(data_2['Time_s'],data_2['UV_1'], color = cycle[1], label = 'Sample 2')
ax.errorbar(data_2['Time_s'],data_2['UV_1'],data_2['Err_1'],fmt ='none',color = cycle[1])
ax.plot(data_2['Time_s'][analysis_2['Fit_start_1']: analysis_2['Fit_end_1']], \
        data_2['Fit_1'][analysis_2['Fit_start_1']: analysis_2['Fit_end_1']], linewidth = 2, color = 'k')
sub_axes = plt.axes([.505, .5, .4, .4])
sub_axes.bar(bar_labels, atpase_avg, color = colors, alpha = 0.9)
sub_axes.errorbar(bar_labels, atpase_avg, yerr=atpase_std,fmt ='None',color='k')
sub_axes.tick_params(which = 'both', direction = 'in', \
                      labelbottom=True, labeltop=False, labelleft=True, labelright=False, \
                      bottom=True, top=True, left=True, right=True)
sub_axes.set_facecolor('whitesmoke')
sub_axes.tick_params(axis='x', labelsize=10)
sub_axes.tick_params(axis='y', labelsize=10)
sub_axes.set_xticks(bar_labels, ('Sample 1', 'Sample 2'))
sub_axes.set_ylabel('pmol P$_i$ min$^{-1}$ / \n pmol Protein', fontsize = 10)
sub_axes.set_ylim(0,32)
ax.tick_params(which = 'both', direction = 'in', \
                      labelbottom=True, labeltop=False, labelleft=True, labelright=False, \
                      bottom=True, top=True, left=True, right=True)
ax.set_xlim(-100,3000)
ax.set_ylim(0,0.8)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Abs$_{340nm}$')
ax.set_facecolor('whitesmoke')
ax.legend(loc = (0.65,0.05), frameon = False)
plt.tight_layout()
plt.savefig('atpase_simulation.png', dpi=300)
plt.show()

# Write dataframes to csv and ATPase activity in report
data_1.to_csv('atpase_data_1.csv')
data_2.to_csv('atpase_data_2.csv')

with open('atpase_analysis_1.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in analysis_1.items():
       writer.writerow([key, value])
with open('atpase_analysis_2.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in analysis_2.items():
       writer.writerow([key, value])