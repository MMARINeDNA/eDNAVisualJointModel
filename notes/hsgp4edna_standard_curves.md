# Request for HSGP4eDNA: estimate the qPCR standard curves in the model

Text to paste into HSGP4eDNA's roadmap (eDNAVisualJointModel ROADMAP_2026-10-21, section 5
item 3). Written 2026-10-21.

---

**Estimate the qPCR standard curves inside the eDNA model** (request from
eDNAVisualJointModel, for the real 2019 data). At present `alpha_ct`, `beta_ct`, `kappa`,
`gamma0_ct`, `gamma1_ct` and `sigma0_ct` are fixed in the data block of
`hsgp_2d_bathysp.stan`. For the real data they should be parameters, estimated jointly
from the qPCR standards and the field wells, as in the qPCR vignette's field-sample model.
Structure:
- **Per plate:** an intercept `alpha_ct[p]` and slope `beta_ct[p]`. The slope is bounded at
  `beta_ct >= 1/log(2)`, i.e. at most a doubling per cycle (the vignette's
  `beta1 <= -1/log(2)`). Priors as in the vignette: intercept N(40, 5), slope N(1/log 2, 0.2).
- **Shared across plates:** the Ct-noise terms `gamma0_ct`, `gamma1_ct`, `sigma0_ct` and the
  detection rate `kappa` (logit-normal(1, 3) prior). Each plate has only 14–21 standards, too
  few for its own noise terms.
- **New data:**
  - the standards: their plate index, detect flag, Ct and known copies per µL, with
    λ = 2 µL × copies/µL;
  - a plate index per field well.

  The standards' likelihood is the same hurdle as the field wells', with λ known.
- **Which standards:** all 33 plates, including the 8 with no field wells, because they
  still inform the shared terms. Standards with inhibition ≥ 0.75 are dropped (all the 10⁵
  standards fall in this group).
- **Reference implementation:** eDNAVisualJointModel's `stan/qpcr_standards.stan` (standards
  only) and `R/fit_qpcr_standards.r`. These give, for the 2019 hake data:
  - `alpha_ct` 36.1–39.9 and `beta_ct` 1.45–1.49 (efficiency 0.95–1.00);
  - γ₀ 0.26, γ₁ −0.37, σ₀ 0.14, κ 0.90;
  - posterior SD of about 0.1 Ct for the intercepts, and a within-plate intercept–slope
    correlation of about 0.8.

  The joint estimates should match these closely; a large shift would point to
  misspecification elsewhere. The real eDNA object carries `qpcr_plate_idx` per well and the
  stand-alone curves (`truth$qpcr_curves`).
- **Simulation:** the simulators keep one fixed curve unless per-plate truth is wanted. A
  simulator option to generate standards from the curve would let the joint estimation be
  tested on simulated data.
- **Related:** the planned bound on the Ct-noise SD at very low copies (§3) should be
  revisited with the γ₁ estimated from the real standards.
