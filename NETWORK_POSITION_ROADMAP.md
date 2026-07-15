# Network-Position Work Plan

This file is only for the HydroRIVERS network-position step.

## Goal

Use HydroRIVERS to describe where each sampling site sits in the downstream
river network.

This started as a distance-to-coast idea, but the current version is broader. It
looks at downstream river distance, ocean connection, drainage-area coverage,
estimated discharge, and stream order.

## What Exists Now

- The network-position script is `06_network-position.R`.
- Methods and output columns are described in `README_network-position.md`.
- Preliminary output tables are written to `outputs/network-position`.
- Preliminary tables were uploaded to the shared Drive folder.
- The current analysis compares a few possible ways to define near-mouth sites,
  using distance and drainage-area coverage cutoffs.

## Next Steps

1. Review `network_position_site_metrics.csv`.
2. Review `network_position_cutoff_comparison.csv`.
3. Decide which near-mouth definition is most useful.
4. Check the sites listed in `network_position_sites_needing_review.csv`.
5. Follow up with Abra about the distance-to-coast metric she had in mind.
6. Decide whether the network-position variables belong in the current
   long-term paper analysis or should stay as background for future flux work.
7. If the variables are used in the paper, add a short methods note describing
   HydroRIVERS matching and the main limitations.

## Useful References

- HydroRIVERS documentation
- Lehner and Grill (2013)
- Tank et al. (2020)
- McNicol et al. (2023)
