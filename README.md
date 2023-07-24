# meso-scavenger

Supplemental materials for "Predicted short-term mesoscavenger release gives way to apex-scavenger dominance"

Authors: J. G. Donohue, P. T. Piiroinen and A. Kane 

Corresponding Author: A. Kane (adam.kane@ucd.ie), UCD School of Biology and Environmental Science, University College Dublin, Ireland


This repository provides a code to solve the scavenging community model for three different sets of parameter values and a code to plot the time-series output.

1. The Matlab code "solvescavengingmodel.m" was written in Matlab R2021b and solves the dynamical system given by Eq. (1). The scenarios 1,2 and 3 correspond to Figure 1 in the main text, Figure S1 in the supplement and Figure S2 in the supplement respectively.
2. The R project code "figures.R" plots the results for the main scenario and two sensitivity analyses. The data are contained within the .csv files.


Also included (in the "SupplementS5" folder) is a code to find relevant equilibrium solutions to the temporal-segregation-of-carrion extension to the original model and a code to plot the output.

1. The Matlab code "temporalsegregation_eqpts.m" was written in Matlab R2021b and follows the stable equilibrium point of the dynamical system given by Eq. (S20) as the n_f parameter is increased from 0 to 1.
2. The R project code "bifurcation_figure.R" plots the results, producing Figure S4 in the Supplement. The data is contained within the "increase_n_f.csv" file.
