# LiteDIC vs. DICe - algorithmic comparison

Source-level comparison of the two correlation engines, for the
"related work / comparison" section of the LiteDIC paper.

Compared versions. LiteDIC as in this repository (`6_DIC/dic_engine.m`,
`dic_engine_batch.m`, `pls_strain.m`) against DICe (Sandia National Laboratories,
`dicengine/dice`), master snapshot reviewed 2026-07-29.
Every claim below is traceable to a named file - DICe paths are relative to its `src/`.

Companion documents. `validate_dice.m` in this folder performs the *quantitative*
comparison against DICe's official gold solution; this file explains *why* the two codes
differ. Literature backing is in the (non-committed) `literature/` knowledge base,
principally `notes/su2016_speckle_quality_rmse.md`.

---

## 1. Overview

DICe is a general-purpose DIC platform: local and global (FE-based) correlation,
occlusion handling, MPI parallelism, multiple image formats, a stereo camera model with
several distortion families, and a GUI - roughly 58,000 lines of C++ requiring
Trilinos, OpenCV, NetCDF and HDF5.

LiteDIC is a single-purpose local subset engine: ~430 lines of MATLAB for the
numerical core (`dic_engine.m` + `pls_strain.m`; ~10,000 lines for the whole repository
including validation, stereo, demos and tests), with no toolboxes and no build step.

The two are therefore not competitors in scope. Where they *are* directly comparable -
the local subset kernel - LiteDIC makes several choices that are, by the current
literature, the more conservative ones. Those choices are enumerated in §2 and their
consequences in §7.

---

## 2. The correlation kernel - where the mathematics actually differs

| | LiteDIC | DICe |
|---|---|---|
| Optimisation scheme | Inverse-compositional Gauss-Newton. Hessian `H = J'J` built from *reference* subset gradients once per point; update `W ← W · inv(W(Δp))` | Forward-additive Gauss-Newton. `H` re-accumulated every iteration and inverted by LAPACK `GETRF`+`GETRI` each time; update `p += Δp` (`Local_Shape_Function::update`) |
| Gradient source | reference frame, 5-point central stencil `[1 −8 0 8 −1]/12`, `O(h⁴)` | default `compute_def_gradients = false` → reference gradients rotated by θ (`Gx = cosθ·gx − sinθ·gy`); with the flag on, deformed-frame gradients. Stencil: 2-point central, optionally smoothed by a 5-point convolution |
| Residual minimised by the solver | full ZNSSD: `e = (fn/gn)·gz − fz` - invariant to intensity offset and gain | `GmF = (G − meanG) − (F − meanF)` - zero-mean SSD, offset-invariant only. The `denomF/denomG` normalisation exists in the source but is commented out (`core/DICe_Objective.cpp`). Full ZNSSD appears only as the reported metric `Subset::gamma()` |
| Convergence test | `‖Δp‖` weighted by the subset half-width over all six parameters, `P.conv_tol = 1e-3` | `max\|Δpᵢ\| < tol` over u, v and θ only - the strain parameters are not tested (`Affine_Shape_Function::test_for_convergence`) |
| Regularisation | none | optional Levenberg-Marquardt term `H(0,0) += α`, `H(1,1) += α`, and an optional momentum term on the update |
| Alternative solvers | none | `SIMPLEX` (Nelder-Mead, `core/DICe_Simplex.cpp`), `GRADIENT_THEN_SEARCH`, `SIMPLEX_THEN_GRADIENT_BASED`, `GRADIENT_BASED_THEN_SIMPLEX` |

### Why the scheme matters

Shao, Dai & He (2015), as summarised by Su et al. (2016) §5.3, show that the inverse
algorithm introduces no noise-induced bias, whereas the forward algorithm does. In the
error decomposition of Su et al. (2016) Eq. (12),

```
RMSE²(u₀) ≈ (σ_f² + σ_g²)/Σf′²  +  (C_ib + C_nb)² sin²(2π u₀)
```

the noise-induced-bias amplitude `C_nb = (σ_g²/f′²̄)·a₁` vanishes for LiteDIC (their
Eq. 28), leaving interpolation bias `C_ib` as the sole systematic term. DICe's default
gradient solver is forward-additive with a zero-mean SSD residual, i.e. the configuration
that retains `C_nb`; that term grows with the square of the image-noise level.

This is the main methodological difference between the two codes, and it favours
LiteDIC. Cite the sources in the paper.

---

## 3. Subset shape function

| | LiteDIC | DICe |
|---|---|---|
| Parameterisation | classical affine `p = [u uₓ u_y v vₓ v_y]`, 3×3 warp matrix | physically decomposed affine `(u, v, θ, εxx, εyy, γxy)`; individual DOFs can be switched off (`rotation_enabled`) |
| Orders available | first only | `AFFINE_SF` (6), `QUADRATIC_SF` (12), `PROJECTION_SF` (homography, 8), `RIGID_BODY_SF` |

DICe's quadratic and projective warps have no LiteDIC counterpart. Also a deliberate
choice: Su et al. (2016) §5.2, Eq. (27) - citing Schreier & Sutton
(2002), Xu et al. (2015) and Wang & Pan (2015) - report that a second-order shape function
doubles the standard-deviation error relative to zero- and first-order. It reduces
under-matching bias where deformation gradients are strong, at the cost of a 2× noise floor
everywhere. For a code whose stated aim is a low noise floor on smooth fields, first order
is the correct default; quadratic belongs in future work, with the penalty stated.

The projective warp remains a genuine gap for strongly out-of-plane / oblique 2D scenes.

---

## 4. Interpolation and image handling

| | LiteDIC | DICe |
|---|---|---|
| Interpolant | `griddedInterpolant(..., 'spline')` - cubic spline over the whole frame, built once per pair | `BILINEAR`, `BICUBIC`, `KEYS_FOURTH` (Keys biquartic), evaluated per pixel. DICe's own `full_field` example uses `KEYS_FOURTH` |
| Gradient method | 5-point central `O(h⁴)` | `FINITE_DIFFERENCE`, or `CONVOLUTION_5_POINT` (finite difference then 5-point smoothing) |
| Bit depth | any integer class scaled by its own `intmax`; floats normalised by their max | handled per format reader |
| Formats | whatever `imread` accepts (TIFF/PNG/JPEG…) | TIFF, PNG, JPEG plus `cine`, `netcdf`, `rawi`, and satellite imagery |

Interpolation is now the *only* systematic-error lever LiteDIC has left (see §2), so the
remaining headroom is cubic → quintic B-spline: Schreier et al. (2000) bound the relative
apparent strain at ~4-5 % for cubic and <2 % for quintic. Su et al. (2016) Eq. (7) gives
the closed-form transfer functions needed to compute the bias for either.

Note on pre-filtering: Su et al. (2016) §1, citing Baldi & Bertolino (2015), warns that
low-pass filtering removes gradient content and therefore raises the random term. It is
a bias/variance trade, not a free improvement.

---

## 5. Initialisation and traversal strategy

LiteDIC. One integer ZNCC search at the seed - optionally centred on a coarse prior
(`P.seed_shift`) and placed at a chosen location (`P.seed_xy`) - then reliability-guided
propagation: a binary max-heap ordered by ZNCC, with first-order (Taylor) transfer of
the converged warp to each 4-neighbour (Pan, Li & Tong 2013):

```
u ← u + uₓ·dx + u_y·dy      (gradients carried over unchanged)
```

A node is finalised only on success, so a failure from one neighbour does not block
convergence from another. When a correlation barrier (weld bead, tape, glare band) crosses
the entire grid, the engine re-seeds at the highest-texture unreached node
(`best_texture_node`, up to `P.max_seeds`); on barrier-free data no re-seed triggers and
the result is bit-identical to single-seed RG-DIC.

DICe. A pluggable `Initializer` hierarchy (`core/DICe_Initializer.h`):
`USE_FIELD_VALUES`, `USE_NEIGHBOR_VALUES`, `USE_NEIGHBOR_VALUES_FIRST_STEP_ONLY`,
`USE_PHASE_CORRELATION`, `USE_OPTICAL_FLOW`, `USE_ZEROS`, `USE_FEATURE_MATCHING` (OpenCV),
`USE_IMAGE_REGISTRATION`, `USE_SATELLITE_GEOMETRY`, `USE_SPACE_FILLING_ITERATIONS`,
`USE_PLANAR_PROJECTION`, `USE_RECTIFIED_CORRESPONDENCES`, plus `Path_Initializer` for
repeating trajectories and `Search_Initializer`.

Assessment. DICe's breadth is real and useful for hard acquisition cases (very large
motion, feature-poor scenes, satellite imagery). Against that, DICe's neighbour
initialisation transfers the neighbour's field value directly, whereas LiteDIC extrapolates
it one grid step with the measured gradients - a closer warm start under large strain or
rotation. LiteDIC's large-displacement path is the coarse estimator (`dic_coarse.m`) plus
`P.seed_shift`, validated on the four-point-bending series (3169 px total motion).

---

## 6. Strain computation and post-processing

| | LiteDIC | DICe |
|---|---|---|
| Primary method | PLS - pointwise least-squares plane fit in a window specified in pixels, differentiated analytically (`pls_strain.m`) | VSG - virtual strain gauge, same least-squares principle (`VSG_Strain_Post_Processor`) |
| Large rotations | optional Green-Lagrange | Green-Lagrange in `core/DICe_Strain.cpp` |
| Boundary behaviour | strains computed only at valid nodes; window must be covered along both axes and contain ≥ win²/3 points - no extrapolation into masked regions | - |
| Performance | separable fast path for fully valid interior windows (exact, not an approximation) | - |
| Additional processors | - | `NLVC` (nonlocal vector calculus) strains, `Distortion_Correction`, `Altitude`, `Crack_Locator`, `Plotly_Contour` |

The strain mathematics is essentially the same; LiteDIC's edge handling is stricter, DICe's
menu is wider.

---

## 7. Quality metrics and failure handling - LiteDIC's clearest deficit

| | LiteDIC | DICe |
|---|---|---|
| Per-point quality | `R.Q` (ZNCC), `R.NIT` (iterations), threshold `P.zncc_min` | `gamma` (ZNSSD residual), `sigma` (uncertainty estimate from image noise), `beta` (sensitivity), Hessian condition number, full `Status_Flag` taxonomy |
| Solver gates | `rcond(H) < 1e-12` rejection; subset-leaves-frame check | `skip_solve_gamma_threshold`, `initial_gamma_threshold`, `final_gamma_threshold`, `disp_jump_tol`, `theta_jump_tol`, `cond_2x2 > 1e12 → HESSIAN_SINGULAR` |

DICe reports a per-point uncertainty; LiteDIC does not. This is the cheapest gap to close:
the random-error formula of Su et al. (2016) Eq. (10),

```
Var(u) = (σ_f² + σ_g²) / Σf′²
```

uses `Σf′²`, which LiteDIC already computes - it is `H(1,1)` and `H(4,4)` of the
Gauss-Newton Hessian assembled in `icgn_point`. Only an image-noise estimate is missing.
Tracked as ROADMAP item 7a.

---

## 8. Capabilities present in DICe and absent from LiteDIC

Recorded for completeness; none of these are accuracy features.

- Global (FE-based) DIC - `src/global/`: Horn-Schunck, mixed Horn-Schunck,
  Levenberg-Marquardt, Lehoucq-Turner, unregularised, method of manufactured solutions;
  triangular mesh, boundary-condition manager, preconditioners.
- Conformal subsets - polygon, circle, rectangle and boolean combinations
  (`base/DICe_Shape.h`, `Conformal_Area_Def`). LiteDIC masks at grid-node granularity
  (`P.mask`, `P.roi`), not per pixel inside the subset.
- Occlusion handling - `obstructing_subset_ids`, subset evolution as pixels become
  visible, obstruction skin factor.
- Tracking routine - `TRACKING_ROUTINE` with motion windows, per-subset parameters,
  skip-solve flags.
- Parallelism - MPI decomposition (`core/DICe_Decomp.cpp`) over Trilinos
  Epetra/Tpetra. LiteDIC's answer is `dic_engine_batch.m`, which vectorises all points into
  paged matrix operations so MATLAB's built-in multithreading uses every core on one
  machine - comparable wall-clock benefit on a workstation, no cluster, no build.
- Stereo camera model - several lens-distortion families
  (`K1R1_K2R2_K3R3`, `K1R2_K2R4_K3R6`, `K1R3_K2R5_K3R7`, OpenCV 4/5/8/12/14-coefficient
  sets). LiteDIC's `11_Stereo/calibrate_zhang.m` fits `fx, fy, u0, v0` plus two radial
  terms only (`k1, k2`), no tangential distortion; triangulation by DLT.
- Output/format ecosystem - Exodus mesh output, NetCDF, cine, GUI.

## 9. Capabilities present in LiteDIC and absent from (or weaker in) DICe

- Full ZNSSD inside the solver, hence no noise-induced bias (§2).
- Cubic spline interpolation as the default, against Keys/bicubic.
- First-order Taylor transfer of the initial guess to neighbours (§5).
- Multi-seed restart across correlation barriers with texture-ranked re-seeding (§5).
- Batch-vectorised IC-GN for dense sweeps (`dic_engine_batch.m`), numerically identical
  to the point-by-point engine, giving multi-core speedup without MPI or a build system.
- Strain evaluation restricted to valid nodes, avoiding fabricated strain peaks at mask
  boundaries.
- Zero-dependency reproducibility: base MATLAB, no toolboxes, no compilation. The whole correlation
  kernel is short enough to read in one go.

---

## 10. Accuracy ceiling and physical units

Both codes are limited by the same physics.
Su et al. (2016) §1 puts the total error of subset-based DIC at "generally several
hundredth pixels", and their Eq. (12) shows it is set by image noise and intensity
gradients, not by the optimiser. Algorithmic work of the kind listed in the roadmap moves
the displacement floor by a factor of roughly 1.5-2, not by an order of magnitude.

Physical resolution is

```
resolution [nm] = RMSE [px] × scale [nm/px]
```

so nanometre-level figures follow from the field of view, not from the correlation
algorithm. At 0.01 px, a 100 mm field on a ~2500 px sensor gives ~400 nm; ~40 nm requires a
10 mm field; ~10 nm requires a field of a few millimetres, i.e. macro or microscope optics.
Any published nanometre figure - from any code - should be read together with the field of
view it was obtained at. LiteDIC's own stereo comparison against Vic-3D
(`11_Stereo/compare_vic3d.m`) agrees at the micrometre level on that rig, which is what
this arithmetic predicts for its field of view; quote that script's current output rather
than a number copied from here.

---

## 11. Quantitative agreement, as measured

`validate_dice.m` compares LiteDIC against DICe's official gold solution for the
`full_field` example, which is bundled here (`dice_reference/DICe_solution_0.txt`) and is
byte-identical to `6_DIC/plate_ref.tiff` / `plate_def_01.tiff`.

DICe reference settings: subset 27, step 10, ZNSSD, `KEYS_FOURTH` interpolation,
optimisation `GRADIENT_BASED_THEN_SIMPLEX`, VSG strain window 50 px. LiteDIC runs subset 41
with spline interpolation and gradient-based IC-GN only, so exact agreement is not expected;
the acceptance criterion is `|bias| < 0.02 px` and `RMS < 0.05 px` per frame. Refer to the
script's own output for the current numbers rather than quoting them here.

Note that the gold solution therefore comes from DICe's hybrid solver - gradient-based
first, Nelder-Mead simplex as fallback - not from the pure gradient solver analysed in §2.
The stereo counterpart (`11_Stereo/validate_dice_stereo.m`) compares against DICe's stereo
gold on the six `stereo_d_spec` pairs.

---

## 12. Positioning for the paper

Claim for the paper (not "better than DICe"):

> LiteDIC is a compact, dependency-free reference implementation of local subset DIC that
> reaches accuracy comparable to an established research code on a standard benchmark,
> while making two methodologically conservative choices - a fully normalised ZNSSD
> residual inside an inverse-compositional solver (which removes the noise-induced bias
> term) and spline intensity interpolation - and while remaining small enough to be read
> and audited in full.

Limitations to state with it: first-order shape function only
(with the 2× random-error penalty of second order given as the reason), no global DIC, no
occlusion handling, node-granularity masking, and a two-term radial distortion model in the
stereo path. Hild & Roux (2012) is the appropriate citation for positioning LiteDIC as a
local method rather than implementing a global one.

---

## Sources for the claims in this document

- LiteDIC: `6_DIC/dic_engine.m`, `6_DIC/dic_engine_batch.m`, `6_DIC/pls_strain.m`,
  `6_DIC/dic_coarse.m`, `11_Stereo/calibrate_zhang.m`, `9_Validation_DICe/validate_dice.m`.
- DICe: `src/core/DICe_Objective.cpp`, `src/core/DICe_Initializer.{h,cpp}`,
  `src/core/DICe_Schema.h`, `src/core/DICe_PostProcessor.h`, `src/core/DICe_Simplex.cpp`,
  `src/base/DICe.h` (enumerations), `src/base/DICe_LocalShapeFunction.{h,cpp}`,
  `src/base/DICe_Subset.cpp`, `src/base/DICe_Image.{h,cpp}`, `src/base/DICe_Shape.h`,
  `src/global/`, `CMakeLists.txt`.
- Literature: Su, Zhang, Xu & Gao, *Opt. Lasers Eng.* 86, 132-142 (2016),
  `10.1016/j.optlaseng.2016.05.019`; Shao, Dai & He, *Opt. Lasers Eng.* 71, 9-19 (2015);
  Schreier, Braasch & Sutton, *Opt. Eng.* 39(11), 2915-2921 (2000); Pan, Li & Tong,
  *Exp. Mech.* 53, 1277-1289 (2013); Wang & Pan, *Exp. Mech.* 55, 1717-1727 (2015);
  Xu et al., *Exp. Mech.* 55, 1575-1590 (2015); Hild & Roux, *Exp. Mech.* 52,
  1503-1519 (2012). DOIs staged in `literature/references_to_add.bib` and to be verified on
  Crossref before entering `references.bib`.
