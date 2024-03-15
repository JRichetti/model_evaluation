# Model evaluation: the misuse of statistical techniques when evaluating observations versus predictions

Malcolm McPhee1, Jonathan Richetti2, Barry Croke3, and Brad Walmsley1,4
1 NSW Department of Primary Industries, Livestock Industries Centre, University of New England, Trevenna Road, Armidale, New South Wales, 2351, Australia
2 CSIRO, 147 Underwood Av, Floreat, Perth, 6014, WA, Australia
3 Mathematical Sciences Institute and Fenner School of Environmental and Society, The Australian National University, Australia
4 Animal Genetics and Breeding Unit (AGBU), University of New England, Armidale, NSW, 2351, Australia


This is a repository that have the data and functions used to calculate various metrics to evaluate model performance.

Please cite: 

Title: Model evaluation: the misuse of statistical techniques when evaluating observations versus predictions

Journal: SESMO

Year: 2024

Authors: Malcolm McPhee 1, Jonathan Richetti 2, Barry Croke 3, and Brad Walmsley 1,4


1- NSW Department of Primary Industries, Livestock Industries Centre, University of New England, Trevenna Road, Armidale, New South Wales, 2351, Australia
2- CSIRO, 147 Underwood Av, Floreat, Perth, 6014, WA, Australia
3- Mathematical Sciences Institute and Fenner School of Environmental and Society, The Australian National University, Australia
4- Animal Genetics and Breeding Unit (AGBU), University of New England, Armidale, NSW, 2351, Australia

DOI: xxx

### How to
The metrics are in the Metric_calculator.R - simply call download that and call:
 > source('Metric_calculator.R') # this loads the functions

then you can use the function on a dataframe `df` that has columns for observed and simulated or predicted values and a character or factor column for groups.

 > metrics_calc(df, Observations, Simulations, Group)

Note that the column names are not quoted, if you run with quotes you get an error! Don't do make_calc(df, 'Observations', 'Simulations', 'Datasets')

Check the other .R for three plots that I like to use.

- The standard Scatterplots
- These distribuition ridges plots that look nice.
- And the distribuition of the residuals (differences, Observations - Predictions)
