# Engineering Analysis Lab - MATLAB

A collection of structural mechanics and experimental mechanics analyses I am
building up step by step. The long-term goal is a validated 2D (and later
stereo) digital image correlation system running on a plain CPU, plus the
classical-mechanics models I compare it against. Everything runs on base
MATLAB - no toolboxes required (except the FE beam demo, which needs the PDE
Toolbox).

Quick start: open `START.m` in MATLAB and press Run. A launcher menu opens;
every analysis is one button.

`DICe_examples/` holds local working copies of the two canonical real-world
datasets (mono 2D large-motion series, stereo tension test) used to validate
`6_DIC`/`11_Stereo` against real data. They are large and third-party, so the
folder is git-ignored and they are not redistributed here. Nothing in the
default `RUN_TESTS` run depends on them - it uses only the images bundled in
this repository. The opt-in `RUN_TESTS('stereo-data')` mode does need an
external dataset, and prints how to point it at one.

## Contents

| Folder | Script | What it does |
|---|---|---|
| `1_Plates_CLT_FSDT` | `plate_clt.m` | [0/90/0] laminate bending, classical lamination theory, Navier series |
| `1_Plates_CLT_FSDT` | `plate_fsdt.m` | Same plate with transverse shear (Mindlin/FSDT), coupled 3x3 Navier system |
| `2_Damage_crack` | `matrix_cracking_ddm.m` | Matrix crack density evolution (discrete damage mechanics), saturating stiffness sigmoid with Weibull toughness statistics |
| `3_Cyl_shell` | `cyl_shell.m` | Pressurized sandwich cylinder with closed ends, membrane + edge effect |
| `4_Plate_vibration` | `plate_vibration.m` | Cantilever plate natural frequencies, Rayleigh-Ritz |
| `5_Rocket_strength` | `rocket_loads.m` | Rocket body sectional loads (N, Q, M) and tank stresses |
| `6_DIC` | `dic_coarse.m` | Fast integer-pixel DIC with parabolic subpixel peak fit |
| `6_DIC` | `dic_icgn.m` | Subpixel DIC engine: ZNSSD + inverse-compositional Gauss-Newton + reliability-guided propagation + pointwise-least-squares strains |
| `6_DIC` | `bias_scurve.m` | Interpolation-bias S-curve self-test: Fourier-phase-shift subpixel sweep (Schreier et al. 2000), bias and apparent strain per pattern/subset |
| `7_Sphere_geometry` | `sphere_demo.m` | Geodesic sphere triangulation from a subdivided icosahedron |
| `8_FEM_beam` | `beam_fem.m` | Transient FE vibration of a cantilever beam (PDE Toolbox) |
| `9_Validation_DICe` | `validate_dice.m` | One-button validation of the subpixel engine against DICe's official gold solution |
| `10_Challenge2` | `challenge2_benchmark.m` | Noise floor, spatial resolution and MEI on the DIC Challenge 2.0 Star images (images downloaded separately, see its README) |
| `12_Datasets` | `gen_dataset.m` | Factory of synthetic speckle pairs with EXACT dense ground truth (analytic pattern, interpolation-free warp) for training and testing DIC / DL-DIC |
| `13_Validator` | `validate_external.m` | Audit of any third-party DIC field (universal CSV import) against the bundled DICe gold: bias/RMS verdict, error maps, text report |

## The DIC engine (`6_DIC`)

This is the core of the project. The numerical kernels are plain function
files - `dic_engine.m` (correlation) and `pls_strain.m` (strains) - called
by the `dic_icgn.m` front end and reusable on arbitrary inputs (the DICe
validation and the Challenge 2.0 benchmark call the same kernels).
`dic_engine_batch.m` is a vectorized variant of the same IC-GN solve for
dense benchmark sweeps: all points are advanced simultaneously, so
MATLAB's built-in multithreading uses every CPU core - still base MATLAB,
no toolboxes. The pipeline is the standard high-accuracy local DIC method:

- ZNSSD correlation criterion (robust to illumination changes);
- subpixel registration with inverse-compositional Gauss-Newton - the
  Hessian is assembled once per point, the warp is updated by composition;
- spline interpolation of intensity;
- reliability-guided propagation: an automatic seed, then the solution
  grows to neighboring points in order of correlation quality (binary-heap
  front), so every point starts from a converged neighbor;
- strains by pointwise least squares with the window specified in pixels,
  computed only at valid nodes - this keeps the strain field independent
  of the grid step and free of edge artifacts at the specimen hole.

The script starts with a self-test (a known rigid subpixel shift must be
recovered with zero strain) and then processes a real loading series of a
plate with a hole. The frames are discovered by filename pattern
(`plate_ref.tiff` as the reference, `plate_def_*.tiff` as the series in
sorted order), so any image series can be analyzed by editing the two
pattern lines at the top of the script. The strain fields show the
expected stress concentration lobes at the hole.

The engine is validated against DICe (Sandia): the image series is
byte-identical to DICe's official `full_field` example, and against its
published gold solution the displacement agreement is |bias| < 0.01 px,
RMS < 0.01 px (see `9_Validation_DICe`).

Checkpoints I keep for myself while developing:

- self-test: shift recovery error < 0.02 px, spurious strains < 5e-4;
- loading series: displacement grows monotonically frame to frame,
  mean ZNCC stays ~0.997;
- the strain field must NOT change when the grid step changes (the strain
  window is defined in pixels for exactly this reason);
- strain concentration at the hole edge ~3x nominal, as theory predicts.

## Roadmap

- [x] Classical models (plates, shell, vibration, damage, FE beam)
- [x] Integer-pixel DIC
- [x] Subpixel DIC engine (IC-GN), strain computation, loading series
- [x] Validation against DICe (Sandia): bias < 0.01 px, RMS < 0.01 px vs
      the official gold result on the identical image series (see
      `9_Validation_DICe`); full-series cross-check tooling included
- [x] Synthetic benchmark: DIC Challenge 2.0 Star images (`10_Challenge2`).
      Displacement MEI 0.582 px^2 - on the published local-affine reference
      cluster; strain noise 2-4x below the reference at matched spatial
      resolution (PLS averages noise over the window area before
      differentiating; the scalar strain MEI is not parameter-invariant,
      so it is not quoted as a headline).
      Full sweep (subsets 9-59, strain windows 20/50/80) in 56 min on an
      8-core CPU, base MATLAB
- [x] Interpolation-bias S-curve sweep (`6_DIC/bias_scurve.m`):
      cubic-spline S-curve amplitude 2.2-2.4e-4 px (apparent strain 0.16 %)
      on the self-test speckle, 3.5e-3 px (2.3-2.4 %) on a harsh sigma-1 px
      pattern, subset-independent - the measured driver of the Star 6
      long-period amplitude excess; mini version in RUN_TESTS (51 checks)
- [ ] FE comparison of the hole-plate strain fields
- [x] mm/px calibration, then stereo calibration (`dic_scale.m`,
      `11_Stereo/calibrate_cameras.m` - toolbox-free Zhang on a dot target)
- [x] Stereo (3D) DIC on a two-camera sample (`11_Stereo/`, validated
      against DICe's stereo gold; demo on the GT4 application pair)
- [x] Surface strain from the stereo 3D points (`11_Stereo/surf_strain.m`:
      local tangent-plane least squares, Green-Lagrange; synthetic
      regression tests - rigid body = zero strain, known stretch recovered)
- [ ] CPU parallelization (per-frame and zone-based propagation)

Everything is CPU-only by design: the target machine is an ordinary laptop
with 16 GB RAM.

## License and citation

Released under the [GNU Affero General Public License v3.0 or later](LICENSE):
free to use, study, modify and redistribute, provided that derivative works -
including software offered to users over a network - are released under the
same licence with their complete source. See [`NOTICE`](NOTICE).

A separate commercial licence, without the copyleft obligations, is available
from the copyright holder: iliiawuest@gmail.com.

Citation metadata is in [`CITATION.cff`](CITATION.cff) - GitHub renders it as
a "Cite this repository" button.
