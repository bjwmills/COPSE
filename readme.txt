COPSE V2.1 (Carbon Oxygen Phosphorus Sulfur Evolution)
As used in Tostevin and Mills (2020) Interface Focus
Coded by Benjamin JW Mills // b.mills@leeds.ac.uk

Contents:
COPSE_frontend.m // Model frontend file // run "COPSE_frontend(0)" to solve model and plot results // sets parameters and starting values, runs solver
COPSE_equations.m // Model equations file // do not run this code directly // contains flux and reservoir calculations
COPSE_sens.m // sensitivity analysis // individual sensitivity variables are contained in frontend and equations files // call this script to run sensitivity and plot
COPSE_plot.m // called by frontend to plot results // do not run this code directly
COPSE_plot_sens.m // called by sensitivity analysis to plot results // do not run this code directly

Reproduction of paper plots:
Paper shows a standard sensitivity analysis over 10,000 runs. 

Updates to code:

This version includes an update to the vegetation mass calculation to use low latitude rather than global temperautre.

This version also uses a lower 'pre-plant' weathering enhancement than assumed in the version that is published.

This version experiments with adding a direct P source from the oxidaiton of marine DOM during the Shuram excursion.


Notes:
Requires MATLAB. Tested in R2018a on win10 x64. Run time ~5 seconds. No installation required.
