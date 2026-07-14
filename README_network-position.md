# Network Position Analysis

## Purpose

This analysis measures where each river sampling site sits within the downstream river network. It is intended to help distinguish sites that are far upstream from sites that may reasonably represent material carried toward the ocean.

This is not a straight-line distance to the coastline. Each site is matched to the nearest HydroRIVERS river segment, and distance is measured downstream through the river network.

## Script

`06_network-position.R` reads the two project input files, downloads any missing regional HydroRIVERS files, and writes three CSV files to `outputs/network-position`. A run takes roughly 10–15 minutes when the HydroRIVERS files are already downloaded.

## Current Analysis

The analysis uses:

- `all-data_si-extract_3_20260629.csv`, the current all-site spatial dataset
- `Site_Reference_Table - WRTDS_Reference_Table_LTER_V3.csv`, for site coordinates
- HydroRIVERS version 1.0, for the river network and connected drainage areas

Project inputs are stored in Box under `SiSyn/long-term-change/data` and mirrored in the [shared Google Drive data folder](https://drive.google.com/drive/folders/1F16i4-dMvIvd_jTHKKnbdJALGa4PSPOJ).

Current results:

- 543 sites across 32 LTER or project groups were included in the spatial input
- 446 sites were matched to HydroRIVERS
- 442 matched sites drain to the ocean
- 4 matched sites drain to inland basins
- Depending on the distance and basin-coverage cutoffs, 36 to 69 sites qualify as near-mouth, high-coverage sites

The cutoff comparison tests all combinations of:

- 10, 25, 50, and 100 river-km from the ocean
- 80%, 90%, and 95% of the connected outlet drainage area captured upstream of the sampling site

## Site Groups

- **All matched sites:** Use for network-position analyses and site-level WRTDS flux analyses
- **Ocean-connected sites:** Use to assess potential downstream silica transport
- **Near-mouth, high-coverage sites:** Use when a station needs to serve as a possible estimate of flux near the ocean

## Future Work

Network position could later be compared with WRTDS flux estimates. McNicol et al. (2023) provides background for that analysis.

## Files

### `network_position_site_metrics.csv`

One row per site successfully matched to HydroRIVERS.

| Column | Meaning |
|---|---|
| `LTER` | LTER |
| `Stream_Name` | Stream name |
| `Latitude` | Latitude used to match the site to HydroRIVERS |
| `Longitude` | Longitude used to match the site to HydroRIVERS |
| `HYRIV_ID` | Unique HydroRIVERS identifier for the matched river segment |
| `downstream_to_outlet_km` | River-network distance from the matched segment to its final downstream point, including inland basin endpoints |
| `downstream_to_ocean_km` | River-network distance to the ocean; blank for inland-draining sites |
| `upland_skm` | HydroRIVERS upstream drainage area at the matched segment, in km2 |
| `outlet_upland_skm` | HydroRIVERS upstream drainage area at the final basin segment, in km2 |
| `basin_coverage_fraction` | `upland_skm / outlet_upland_skm`; the fraction of the connected outlet basin captured at the matched segment |
| `reference_drainage_area_km2` | Drainage area recorded in the current spatial dataset |
| `reference_basin_coverage_fraction` | Reference drainage area divided by HydroRIVERS outlet drainage area |
| `discharge_est_cms` | HydroRIVERS long-term estimated discharge at the matched segment, in m3/s |
| `drains_to_ocean` | `TRUE` when the connected river network reaches the ocean |
| `ORD_STRA` | Strahler stream order; headwaters begin at 1 and order increases when streams of the same order meet |
| `snap_distance_km` | Distance between the supplied site coordinate and the matched HydroRIVERS segment |

### `network_position_cutoff_comparison.csv`

One row for each of the 12 distance and basin-coverage combinations.

| Column | Meaning |
|---|---|
| `distance_cutoff_km` | Maximum river-network distance from the ocean for that comparison |
| `basin_coverage_cutoff` | Minimum fraction of the outlet basin captured at the sampling site; `0.80` means 80% |
| `n_all_matched_sites` | Number of sites successfully matched to HydroRIVERS |
| `n_ocean_connected_sites` | Number of matched sites connected to the ocean |
| `n_sites_with_distance_and_coverage` | Number of ocean-connected sites with both required values available |
| `n_near_mouth_high_coverage_sites` | Number meeting both cutoffs |
| `near_mouth_percent_of_sites_with_metrics` | Percentage of eligible sites meeting both cutoffs |

### `network_position_sites_needing_review.csv`

Rows that could not be included cleanly or need a coordinate check. A site may appear more than once if it has more than one issue.

| Column | Meaning |
|---|---|
| `review_issue` | Reason the row needs review: missing coordinates, invalid coordinates, outside HydroRIVERS coverage, or a large snapping distance |
| `LTER` | LTER |
| `Stream_Name` | Stream name |
| `Latitude` | Latitude available for review |
| `Longitude` | Longitude available for review |
| `snap_distance_km` | Distance to the matched HydroRIVERS segment when a match was possible |

## Limitations

- HydroRIVERS is based on a global river network with an approximate 500 m grid at the equator
- Streams smaller than the HydroRIVERS extraction limits may not be represented
- HydroRIVERS quality is lower north of 60 degrees latitude because the underlying elevation data are coarser
- The downstream distance is measured from the outlet of the nearest HydroRIVERS river segment, rather than from the exact sampling location
- Basin coverage describes the fraction of connected drainage area represented by the site, not the amount of silica carried by the river

## References

- [HydroRIVERS technical documentation](https://data.hydrosheds.org/file/technical-documentation/HydroRIVERS_TechDoc_v10.pdf)
- [Lehner and Grill (2013)](https://doi.org/10.1002/hyp.9740), the HydroRIVERS data reference
- [McNicol et al. (2023)](https://doi.org/10.1029/2023GL103024), background for a possible future flux analysis
