### Repository associated with the paper:

### Inspecting the interaction between HIV and the immune system through genetic
Andrea Mazzolini, Thierry Mora, and Aleksandra Walczak


The folders are organized as follows:
- **data** contains the pipeline and the instructions for downloading the data and processing them.
- **turnover_analysis** contains the notebooks for computing the "evolutionary properties" of the virus and the immune system, such as the measures of turnover.
It includes also all the procedures for the statistical analysis of the correlations which lead to the figure 3, 4, S3, S4, S5, S6, S7 of the paper.
- **wright-fisher** contains the population genetics simulations in c++ and the python notebooks for reading the trajectories and generating plots.
- **8state_mc** contains the Mathematica file for the computation of the analytical expression of the turnover correlation in the simple 8-states Markov Chain with equal selection coefficients and a python notebook for simulating the process.

This code can be used for reproducing most of the figure of the paper. Notice however that the code has not been written with the purpose of being easily readable from external users.
If you are interested in reproducing the results or a part of them and you're having hard time in deciphering the code, please send me an email at andrea.mazzolini.90@gmail.com.
