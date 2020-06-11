## Instructions

To reproduce the simulation section

```R
R CMD BATCH ./run_simulation_tl.R
cd ./code_simulation/
# run simulation
mkdir ./output/
mv ./df_metric.rda ./output/df_metric.rda
# create plots
R CMD BATCH ./plot_result.R
# the plots will be saved in `./code_simulation/output/` after the script
```
