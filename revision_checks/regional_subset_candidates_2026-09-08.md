# Regional subset selection audit

No model fits were run for this audit. Candidate selection used only country labels, study IDs, recorded coordinates, and geographic distances.

## Selection rule

Candidate regions were compared by effect-size count, study count, recorded-coordinate count, geographic coherence (longitude/latitude span and maximum site-to-site great-circle distance), and distance distortion under one defensible regional projected CRS. Model estimates, sampling variances, and any package-specific fit were not used.

## Candidate regions

| Candidate | Effect sizes | Studies | Recorded locations | Longitude span | Latitude span | Maximum great-circle distance |
|---|---:|---:|---:|---:|---:|---:|
| US | 1,047 | 152 | 144 | 90.06° | 38.53° | 6,905 km |
| Australia | 333 | 68 | 63 | 37.62° | 29.99° | 3,716 km |
| Spain | 186 | 30 | 32 | 11.96° | 7.50° | 999 km |
| Canada | 153 | 22 | 25 | 93.52° | 17.62° | 5,410 km |
| South Africa | 75 | 17 | 15 | 13.73° | 10.14° | 1,626 km |
| Brazil | 70 | 15 | 13 | 8.01° | 27.54° | 3,089 km |
| Argentina | 43 | 10 | 9 | 13.14° | 13.05° | 1,624 km |

Countries with fewer than eight recorded locations were not competitive under the same retention criterion.

## Projection check

For Spain, a WGS84 Lambert Conformal Conic CRS centred at longitude -3.5°, latitude 40.5°, with standard parallels 38° and 43°, gives a maximum absolute pairwise-distance distortion of 0.103% and a 95th-percentile absolute distortion of 0.088% relative to WGS84 ellipsoidal great-circle distances. The corresponding South Africa CRS (standard parallels -33° and -25°; centre -29°, 25°) gives 0.242% maximum and 0.240% at the 95th percentile. Spain therefore retains more data and has lower distortion.

## Selected candidate

Spain is the largest single-country candidate under the prespecified criteria. It retains 186 effect sizes, 30 studies, and 32 recorded locations. Twenty-nine studies have one recorded location and one study (`Fernandez-Garcia_2020`) has three recorded locations; no recorded coordinate location is shared by multiple studies in this subset. Thus the subset retains a within-study multi-location structure, but not the full dataset's cross-study shared-location structure.

An Iberian label-based candidate (`country %in% c("Spain", "Portugal")`) is also defensible: it adds one Portuguese study, three effect sizes, and one recorded location (189 effect sizes, 31 studies, 33 locations). Its maximum great-circle separation remains 998 km and the same Spain-centred LCC has maximum absolute pairwise distortion 0.103%. Because the gain is small and the Spanish subset has a simpler, single-country definition, Spain is retained as the primary candidate; Iberia is recorded as a sensitivity candidate, not fitted in this audit.

## Geographic-connectivity check within Spain

The country label was not used to silently delete isolated coordinates. Using a fixed 300-km great-circle link threshold and retaining the largest connected component gives 31 sites, 185 effect sizes, and 29 studies. The excluded singleton is the recorded coordinate `(37.06, -6.66)`, with one effect size and a 377-km nearest-neighbour gap. At the same threshold the coordinate `(44.55556, -0.111111)` remains connected to the largest component (its nearest-neighbour gap is 255 km). This confirms that the full Spain-labelled subset is largely coherent but contains one sparsely connected southern coordinate; the 186-effect-size full subset is retained unless a pre-specified clustering rule is adopted. No original papers were inspected and no model estimates were used for this check.

No model fits were run. The next gate is to choose between the transparent full Spain-labelled candidate (186/30/32), the objective 300-km largest-component version (185/29/31), or the slightly larger Iberian candidate (189/31/33) before fitting the three packages.
