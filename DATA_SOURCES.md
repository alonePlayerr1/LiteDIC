# Data sources

Every image and reference result used by LiteDIC, where it comes from, and on what terms it is
used. The LiteDIC source code itself is under the GNU AGPL v3 (`LICENSE`); the data below keep
their own terms.

## Bundled in this repository

| Folder | What | Source | Terms |
|---|---|---|---|
| `6_DIC/` (images) | Plate-with-hole image series | DICe, Sandia National Laboratories — <https://github.com/dicengine/dice> (`tests/examples/full_field`) | DICe license, copyright Sandia/NTESS; redistributed unchanged |
| `9_Validation_DICe/dice_reference/` | DICe gold solution and original settings | same | same |
| `9_Validation_DICe/stereo_d_spec/` | Six stereo image pairs, rig calibration, DICe stereo gold solutions | DICe (`tests/regression/stereo_d_spec`) | same |
| `9_Validation_DICe/calibration/` | 66 dot-target calibration pairs of the same two-camera rig | DICe | same |
| `9_Validation_DICe/stereo_sample/` | Two application pairs from the same experiment | DICe | same |

## Used but not redistributed

| Data | Source | Terms | How to get it |
|---|---|---|---|
| DIC Challenge 1.0 and 2.0 image sets and published results | Society for Experimental Mechanics (SEM) and iDICs — <https://sem.org/dicchallenge/> | Free for research use with attribution: P. L. Reu et al., *Exp. Mech.* 58 (2018) and 62 (2022) | Download and place as described in `10_Challenge2/README.md` |
| Four-point-bending images of a composite beam and of a sandwich panel | Mechanical Engineering Research Institute of the Russian Academy of Sciences (IMASH RAN) | Used with written permission of IMASH RAN for the LiteDIC article only; not published | Available from IMASH RAN on reasonable request |

The folders that hold data which is not redistributed (`DICe_examples/`, `2D-Challenge 1.0/`,
`2D-Challenge 2.0/`) are excluded by `.gitignore` and must never be committed.
