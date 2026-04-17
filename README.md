# Evolution_Growth_Stress

Evolution of microbial population under successive rounds of growth, starvation, stress. The corresponding manuscript is titled "Adaptation to extreme stress under
the growth-survival fitness trade-off".

Trehalosesimulation_singleplot.m was used to generate Fig.2 and Fig.3 of the paper. The code follows the dynamics of the population as it undergoes repeated growth cycles with mutations (with the growth rate of each cell determined by its phenotype), followed by resource exhaustion and starvation (accompanied by probabilistic switch to quiescence determined by each cell's phenotype), and eventually exposure to a stress with a phenotype dependent survival probability.

Trehalosesimulation_comparison_LGR_g0.m and Trehalosesimulation_comparison_LGR_trl.m was used to generate Fig. 4 and Fig. 5.

LGR.nb is a mathematica file used to generate Fig. 6 as well as Fig. 5.

Coexistencetest.ipynb was used to generate Fig. 7.
