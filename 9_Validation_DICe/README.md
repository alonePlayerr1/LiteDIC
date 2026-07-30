# Validation against DICe (Sandia)

`validate_dice.m` compares the LiteDIC subpixel engine (`6_DIC/dic_icgn.m`)
with [DICe](https://github.com/dicengine/dice) - an established, validated
open-source DIC code. The reference is DICe's own: the `6_DIC` image
series is byte-identical to DICe's official `full_field` example
(`tests/examples/full_field`), and DICe publishes a gold solution file
for it. That file is bundled here, so the validation needs no DICe install.

`COMPARISON_DICe.md` in this folder is the qualitative companion: a source-level
comparison of the two correlation kernels (solver scheme, criterion, shape functions,
interpolation, initialisation, quality metrics), what each code has that the other does
not, and how to position LiteDIC in the paper. This README and `validate_dice.m`
answer *how closely* the codes agree; that document answers *why they differ*.

## One-button run

Open `validate_dice.m` and press Run (or use the START menu). The script:

1. runs `6_DIC/dic_icgn.m` automatically if its results are missing/stale;
2. compares frame `plate_def_01` against the bundled official gold result
   (`dice_reference/DICe_solution_0.txt`);
3. additionally compares every frame for which a DICe solution file exists
   in `dice_results/`;
4. prints per-frame bias and RMS of the displacement difference with a
   PASS/FAIL verdict (|bias| < 0.02 px, RMS < 0.05 px) and plots agreement.

## Folder contents

- `dice_reference/` - verbatim copies from the official DICe repository
  (`tests/examples/full_field`): the gold solution file plus the original
  `input.xml`, `params.xml`, `subsets.txt`. Copyright Sandia/NTESS; kept as
  the authoritative reference for settings and output format.
- `dice_input/` - the same settings adapted to run DICe on the FULL series.
  With DICe installed: `cd dice_input && dice -i input.xml -v`
  Results land in `dice_results/` and are picked up automatically.
- `stereo_d_spec/` - verbatim copies from the official DICe repository
  (`tests/regression/stereo_d_spec`): six synchronized two-camera image
  pairs `GT4-0000..0005_{0,1}.tif`, the rig calibration `cal.xml`, the
  DICe settings (`input.xml`, `params.xml`), and DICe's stereo gold
  solutions (`gold/DICe_solution_0..5.txt`). Copyright Sandia/NTESS;
  kept as the authoritative reference for the stereo validation.
- `calibration/` - 66 synchronized dot-target image pairs
  (`AMCalB-%04d_{0,1}.jpeg`) of the same two-camera rig; used by the
  toolbox-free Zhang calibration in `11_Stereo/`.
- `stereo_sample/` - an application pair from the same experiment
  (`GT4-0110/0129`), ~19 load steps apart; no reference solution exists
  for it (see `11_Stereo/stereo_demo.m`).

## Stereo validation (`11_Stereo/`)

`11_Stereo/validate_dice_stereo.m` is the stereo counterpart of
`validate_dice.m`: one button, no DICe installation required. It
calibrates the rig from `calibration/` (cached in `cal_stereo.mat`),
runs the stereo pipeline on all six `stereo_d_spec` pairs, and compares
against DICe's stereo gold: cross-camera disparity at the gold subset
centres, temporal displacement, epipolar residuals under the calibrated
geometry, and the triangulated 3D surfaces. Note that DICe's gold
`CROSS_CORR` values are frame-constant (DICe cross-matches the cameras
once, on the reference pair), so the disparity comparison is
apples-to-apples at frame 0, while frames 1-5 exercise the temporal
channel. The same checks run in `RUN_TESTS.m` ('Stereo DIC vs DICe
gold'). `11_Stereo/stereo_demo.m` runs the pipeline on `stereo_sample/`
and reports internal-consistency indicators only.

## Reference settings (official example)

subset 27, step 10, ZNSSD, KEYS_FOURTH interpolation, optimization
GRADIENT_BASED_THEN_SIMPLEX, VSG strain window 50 px. LiteDIC runs with
subset 41 and spline interpolation, so small displacement differences are
expected - hence the 0.02/0.05 px criterion rather than exact equality.

DICe output format (verified against the gold file): `*` comment lines,
then a CSV header (`SUBSET_ID,COORDINATE_X,...`). The parser finds columns
by name and excludes points DICe marks failed (`SIGMA = -1`). Solution
files are numbered by deformed-frame index starting at 0 and are NOT
zero-padded, so frames are matched by the number in the file name, never
by sort order.
