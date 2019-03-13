# iGraph Networks
Scripting to run iGraph networks.

## Notes
This is a pretty basic script to run network analyses using iGraph in R.
In the usage notes below, you will see the basic workflow outlines including:
1. Filtering predictor variables using options of "method" and "cutoff".
  These allow you to fitler based on beta estimate rank or p-value cutoff.
2. Z transform data from heterogenous sources to ensure everything is on the same scale prior to network build wuth the option "zed".
3. Easily exclude or force predictors into the network using the "perturb" or "force" options. For example you can list genes as predictors and simulate knock outs. This will be developed further in the future.
4. Cluster with correlation or euclidean distances via the "cluster" option.
5. Identify communities using edge betweenness, spinglass (aka Potter) or random walks via the "comms" option.

## Single command line
```
Rscript iGraphNetworkRunner.R [outcome] [predictors] [zed] [method] [cutoff] [perturb] [force] [cores] [output] [cluster] [comms]
```

## Example command line
```
Rscript iGraphNetworkRunner.R ppmiForNetwork_pheno.tab ppmiForNetwork_geneExp_symbolNames.tab no p 0.001 no no 2 testMonday euc spin
```

## Detailed command line options
```
outcome = 2 columns ... ID and PHENO as the header, ID is character sample ID and PHENO is a numeric phenotype
predictors = data used to predict the pheno, must be numeric / continuous for best results
zed = "yes" or "no", do you want your data Z transformed after filtering to make everything on the same scale
method = "effect" or "p", filter on either the absolute beta value from regression of each preditcor or its p-value
cutoff = decimal from 0-1, representing the minimum p-value for inclusion of a feature or quantile of effect estimate you want included in the filtering (ie in this case 0.25 will include the top 25% |beta|)
perturb = "no" or a file containing list of predictor names in a single column with no header to be used to make exclusions from data in analysis to simulate knock outs
force = "no" or a file containing list of predictor names in a single column with no header to be used to force inclusion into network estimates
cores = 1 - 24, how many cores you want this to run on?
output = text string to uniquely identify this run
cluster = "corr" or "euc" for correlation or euclidian based clustering
comms = "edge" or "spin" or "walk", denoting edge betweens, Potts spinglass method or random walk to identify communities
```

## Depends on
R > 3.5 with packages data.table, parallel and igraph

### Questions or comments
[Mike Nalls](mike@datatecnica.com)
