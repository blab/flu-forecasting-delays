# Effects of delayed sequence submission and vaccine development on long-term forecast accuracy of seasonal influenza A/H3N2

Building on [previous work forecasting the evolution of seasonal influenza A/H3N2](https://github.com/blab/flu-forecasting), this workflow investigates the effects of delays in both sequence submission and vaccine development on forecasting accuracy.

## Analysis

Run analysis on a SLURM cluster.

``` bash
snakemake \
    -j 20 \
    --profile profiles/slurm
```

## Manuscript

Build manuscript (`manuscript/delays.pdf`).

``` bash
cd manuscript
./build.sh
```
