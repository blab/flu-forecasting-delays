# Timely vaccine strain selection and genomic surveillance improves evolutionary forecast accuracy of seasonal influenza A/H3N2

John Huddleston<sup>1,3</sup> and Trevor Bedford<sup>1,2</sup>

1. Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA
1. Howard Hughes Medical Institute, Seattle, WA, USA
1. Corresponding author (jhuddles@fredhutch.org)

DOI: https://doi.org/10.1101/2024.09.11.24313489 

## Abstract

For the last decade, evolutionary forecasting models have influenced seasonal influenza vaccine design.
These models attempt to predict which genetic variants circulating at the time of vaccine strain selection will be dominant 12 months later in the influenza season targeted by vaccination campaign.
Forecasting models depend on hemagglutinin (HA) sequences from the WHO’s Global Influenza Surveillance and Response System to identify currently circulating groups of related strains (clades) and estimate clade fitness for forecasts.
However, the average lag between collection of a clinical sample and the submission of its sequence to the Global Initiative on Sharing All Influenza Data (GISAID) EpiFlu database is ~3 months.
Submission lags complicate the already difficult 12-month forecasting problem by reducing understanding of current clade frequencies at the time of forecasting.
These constraints of a 12-month forecast horizon and 3-month average submission lags create an upper bound on the accuracy of any long-term forecasting model.
The global response to the SARS-CoV-2 pandemic revealed that modern vaccine technology like mRNA vaccines can reduce how far we need to forecast into the future to 6 months or less and that expanded support for sequencing can reduce submission lags to GISAID to 1 month on average.
To determine whether these recent advances could also improve long-term forecasts for seasonal influenza, we quantified the effects of reducing forecast horizons and submission lags on the accuracy of forecasts for A/H3N2 populations.
We found that reducing forecast horizons from 12 months to 6 or 3 months reduced average absolute forecasting errors to 25\% and 50\% of the 12-month average, respectively.
Reducing submission lags provided little improvement to forecasting accuracy but decreased the uncertainty in current clade frequencies by 50\%.
These results show the potential to substantially improve the accuracy of existing influenza forecasting models by modernizing influenza vaccine development and increasing global sequencing capacity.

Supplemental data are available on Zenodo at https://zenodo.org/records/13742375.

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
