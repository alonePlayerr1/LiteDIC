# Changelog

All notable changes to LiteDIC are listed here. Versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- The repository holds LiteDIC only. The structural mechanics examples
  (plates, shells, vibration, rocket loads, sphere, FE beam) moved to a
  separate repository, `mechanics-lab`; `RUN_TESTS` keeps the 118 DIC and
  stereo checks.
- Licence: AGPL-3.0-or-later instead of MIT.
- `dic_rainbow`: red at the minimum, violet at the maximum, as in Vic-3D.
- `discover_pairs`: a folder with several series is ordered series by series.

### Added
- `dic_engine`: per-node displacement uncertainty `R.SU`, `R.SV`.
- `icgn_point`: convergence flag.
- `field_overlay`: contour bands with the Vic-3D colormap.

### Fixed
- `run_series_2d`: the ROI follows the specimen, and the rigid motion is fitted
  on converged points only.
- `icgn_point`: the per-point uncertainty could become complex when the
  noise-corrected Hessian was indefinite.
- `surf_strain`: the sign of the surface shear strain depended on the LAPACK
  eigenvector sign and could flip from node to node.
- `stereo_extrinsics`: the rig-rigidity check is measured on the independent
  per-view Zhang poses (after bundle adjustment it was identically zero);
  `rot2vec` is accurate near 180 degrees.
- `field_overlay`: colour limits applied after plane subtraction.
- `run_stereo_dataset`: fixes from the 122-state plate series.

## [1.0.0] - 2026-08-26

First tagged version.

### Added
- 2D subpixel DIC engine: ZNSSD, IC-GN with first-order shape function,
  spline interpolation, reliability-guided propagation with multi-seed restart.
- Pointwise-least-squares strains with a strain window defined in pixels.
- Validation harness against the DICe gold solution; DIC Challenge 2.0 benchmark.
- Interpolation-bias S-curve self-test.
- Stereo layers: dot-target detection, Zhang calibration, stereo matching,
  triangulation, surface strain.
- Regression suite `RUN_TESTS.m`.
