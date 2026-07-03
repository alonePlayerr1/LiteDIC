# DIC Challenge 2.0 benchmark

Measures LiteDIC's metrological performance on the community-standard
DIC Challenge 2.0 "Star" images (Reu et al. 2022, *Experimental
Mechanics*) and scores it with the Metrological Efficiency Indicator
(MEI), so the engine can be compared directly against the ~20 published
commercial and open-source codes.

## What is computed

Along the center-row line-cut (1 px spacing, straight from the DIC grid):

- Measurement resolution `n` - 1σ of the quantity of interest on the
  undeformed noise-floor image pair (pure noise, no motion).
- Spatial resolution `ℓ10%` - the spatial period at which the measured
  sinusoid amplitude first drops 10 % below ground truth (±0.5 px
  displacement for Star 5, ±5 % Lagrangian ε_yy for Star 6). The Star
  field is a vertical sinusoid whose period λ(x) grows linearly across
  the image (verified directly on the images by optical flow); the
  measured amplitude at each 1-px column is extracted by a phase-free
  harmonic fit over a band of rows, and ℓ10% is read from a 12th-order
  polynomial fit of amplitude ratio vs period, as in the paper.
- MEI = `n^o · ℓ10%` with `o = 1` for displacement and `o = 2` for
  strain (lower is better). The headline number is the average of the
  three lowest MEI values, per the paper.

The sweep covers subset sizes 9-59 px and, for strain, PLS strain windows
20/50/80 px. The points are solved with the batched IC-GN kernel
(`6_DIC/dic_engine_batch.m`) - the same mathematics as the sequential
RG-DIC engine, vectorized over all points so MATLAB's built-in
multithreading uses every CPU core, still with no toolboxes. A preflight
asserts batch/sequential agreement on the actual images before anything
is measured. Full sweep on a modern 8-core desktop: about an hour - 56 min
measured (`quick_mode = true` gives a fast sanity pass; `use_batch = false`
forces the sequential engine, hours).

## Getting the images

Download the Star image sets from <https://sem.org/dicchallenge/> and
place them in a folder named `2D-Challenge 2.0` next to this folder, so
that e.g. `2D-Challenge 2.0/Star5LargeNoisy/Star5LargeNoisy/*.tif`
resolves. The images are free to use with attribution; they are not
committed to this repository (gitignored, ~50 MB).

If the published per-code results are also present
(`2D-Challenge 2.0/Results/...`), the benchmark overlays LiteDIC on the
`Local.Affine.B` reference curve - the "straightforward" local affine
implementation the paper uses as its reference. Note: in the official
download, the two Results folder names are swapped relative to their
content (`DisplacementResults_Star6` holds Star 5 displacement data,
`StrainResults_Star5` holds Star 6 strain data); the script accounts for
this.

## Running

```matlab
cd 10_Challenge2
challenge2_benchmark
```

Outputs (in this folder):

- printed per-subset table of `n`, `ℓ10%`, MEI, and the summary MEI;
- `challenge2_results.mat` - all results + raw line-cuts;
- `LiteDIC_linecuts_*.csv` - line-cuts in the Challenge submission format;
- figures: measured line-cuts per subset, and noise vs spatial resolution
  (the paper's Fig. 9/12 axes) with the published reference overlay.

## Attribution

Images and methodology: P. L. Reu et al., "DIC Challenge 2.0: Developing
Images and Guidelines for Evaluating Accuracy and Resolution of 2D
Analyses - Focus on the Metrological Efficiency Indicator," *Experimental
Mechanics* 62 (2022). Challenge 1.0: Reu et al., *Experimental Mechanics*
58 (2018). Images courtesy of the Society for Experimental Mechanics DIC
Challenge board.
