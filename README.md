# atpase_activity_simulations
Simulation of protein ATPase activity assay curves and analysis

Here we just simulate what happens when you add you ATPase protein to an ATP regenerating solution (experiment designed by Norby et al., Methods in Enzymology, 1988). Basically, you measure the depletion of NADH by absorbance at 340nm, and use it to infer how fast your protein hydrolyzes ATP. In this script, we set a few parameters, such as volume, protein concentration etc (I have set two different conditions for comparison), the the script simulates data and do the analysis for you, which includes linear least squares fitting of the linear decay of NADH absorbance and the conversion to ATPase activity. The simulations are done in triplicate and the final result contains the average and standard deviation of the set of three runs for each sample.

The plot is the simulated data for one of the replicates, and the inset contains a bar plot with errorbar for the average ATPase activity. The results are saved in a plot, a csv dataframe for the simulated results and fits, and a csv file por the activity parameters (average, stardard deviation and linear fitting range).
