# vignettes/

Self-contained HTML write-ups of individual model versions. Each
vignette explains the simulation, the statistical model (with
equations), and the diagnostic + posterior-predictive output produced
by the corresponding pipeline run.

## What's here

| File | Source | Notes |
|---|---|---|
| `v3_vignette.html` | `outputs/whale_edna_output_v3/v3_vignette.qmd` | v3 simulation + HSGP joint qPCR / metabarcoding fit. |
| `v3.2_vignette.html` | `outputs/whale_edna_output_v3.2/v3.2_vignette.qmd` | v3.2 debugging case study. Eight PRs of iterative diagnosis on the v3 sampler pathology, walked through chronologically. Final fit + lessons learned. |

The HTML files are **self-contained** (Quarto's `embed-resources:
true`) — figures and styles are inlined, so you can drop them into a
browser, slack, or share by email without bringing additional assets
along.

## Regenerating

Each vignette's source `.qmd` lives next to its model artefacts in the
matching `outputs/whale_edna_output_v{N}/` folder. Re-render with:

```bash
# from the repo root
quarto render outputs/whale_edna_output_v3/v3_vignette.qmd --to html
cp outputs/whale_edna_output_v3/v3_vignette.html vignettes/v3_vignette.html
```

The qmd reads the per-version diagnostic PNGs and CSVs in place
(`diag_energy.png`, `param_recovery.png`, `sampler_summary.csv`,
etc.). If you've re-run `04_run_*` and `05_check_*` for that version,
just re-render and copy.

## Keeping figures up to date

The "simulated truth" page-by-page figures used in each vignette
(`outputs/whale_edna_output_v{N}/_vignette_figs/sim_page-*.png`) are
extracted from the per-version
`simulated_edna_fields_v{N}.pdf` plot via:

```bash
mkdir -p outputs/whale_edna_output_v3/_vignette_figs
pdftoppm -r 150 -png outputs/whale_edna_output_v3/simulated_edna_fields_v3.pdf \
  outputs/whale_edna_output_v3/_vignette_figs/sim_page
```

This requires `poppler` (`brew install poppler`).
