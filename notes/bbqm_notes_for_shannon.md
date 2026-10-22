# Notes on bbqm for Shannon

Notes from reading [skachel/bbqm](https://github.com/skachel/bbqm)
(`Code/Share/bbqm_marginalized_stan`, as of 2026-10-21) while planning to port
the mechanistic φ into HSGP4eDNA (plan: `ROADMAP.md` in HSGP4eDNA). These are
things you might want to revisit. They come in three groups: small fixes in
the shared code, ideas from fitting it to our data, and questions.

## 1. Small fixes in the shared code

1. **β0 is described but not in the model.** The README (model equation,
   parameter table, "Expected output"), `simulate.R`'s header comment,
   `prepare_stan_data()` (`prior_beta0_mu`, `prior_beta0_sig`) and
   `init_stan()` (`beta0`) all refer to the global offset β0, but
   `single_marker_qm.stan` has no `beta0`: log φ = log φ_mech + ε. rstan
   ignores the extra data, so nothing breaks, but a reader can't tell which
   version is current. Either add β0 back (we plan to keep it, as the check on
   the mechanism you describe) or drop the references.
2. **`run_example.R` uses `sm` before defining it.** The diagnostics line calls
   `max(sm[, "Rhat"])` and `min(sm[, "n_eff"])`, but
   `sm <- rstan::summary(fit)$summary` is the script's last line. In a fresh
   session it stops with "object 'sm' not found". Moving that line up to just
   after `extract()` fixes it.
3. **`run_example.R` reports only σ_ε.** The README's "Expected output"
   describes β0 near 0. Once item 1 is settled, print whichever departure
   terms the model has.
4. **Repository housekeeping.** `.Rhistory`, `.Rproj.user/`, several
   `.DS_Store` files and a 630 KB `README.html` are committed. A `.gitignore`
   entry for each would keep them out.

## 2. Ideas from fitting it to our data

Our setting: 551 seawater bottles, 3 target species (hake, humpback, Pacific
white-sided dolphin) and everything else lumped as "junk", usually 3
metabarcoding replicates per bottle (1–6), and hake qPCR as an independent
absolute anchor.

1. **The taxon set doesn't have to be exhaustive.** If a taxon is dropped by
   removing its reads from its column *and* from the replicate total, the
   model stays exact. The reads are multinomial on the amplicon pool, so the
   reads of a subset of taxa, given their total, are multinomial on the
   subset's pool. The subset's shares w_s G_s / Σ_subset w_t G_t don't depend
   on the dropped taxa's copies, and `mech_phi_pi` for the subset uses only
   the subset's moments. So `prepare_stan_data()` could set
   `reads = rowSums(count)` (dropping aliquots left with no reads) instead of
   warning. This is useful for contamination: in our data *Homo sapiens* is
   27% of non-target reads, in 1,438 of 1,709 replicates.
2. **A lumped category can have its own free efficiency.** For a focal taxon,
   the other taxa enter `mech_phi_pi` only through Σ w_t F_t and
   Σ w_t² F_t. A single lumped category with efficiency w_J and expected
   copies F_J matches both exactly when
   `w_J = Σ w_t² F_t / Σ w_t F_t` and `F_J = (Σ w_t F_t)² / Σ w_t² F_t`. So
   the many rare or poorly characterised taxa can be one category with w_J
   free, instead of needing a known `amp_alpha` for every taxon. The
   approximations: w_J depends on the lumped composition, which can differ
   between samples, and F_J is an effective copy number no larger than the
   true one. We plan to build this first (292 non-target taxa in our data;
   17 make up 90% of their reads).
3. **Replicates at different dilutions.** Our metabarcoding replicates are run
   at 1–100× dilution, sometimes different within one bottle. Because
   `mech_phi_pi` already takes copies per observation, a dilution is just
   F_obs = λ_j / d_o, an offset in the data. Since φ scales with the total
   copies, a dilution series within a sample might also sharpen the absolute
   scale. Could be worth a simulation.
4. **Uncertain efficiencies.** `amp_alpha` is fixed. If the mock-community
   fit leaves real uncertainty, a prior on `amp_alpha` (or a joint fit with
   the mock-community data) would carry it into λ. That matters most where
   the efficiencies differ most.
5. **Structure for the departure terms.** The full ε[j,i] gives every
   sample × taxon its own φ offset. In our model that's about 2,200 extra
   latent variables. Many are informed only by the prior (the cetaceans are
   in 3–5% of replicates), and a large σ_ε would take away the information
   about absolute scale that the model is built on. Cheaper options with
   more structure, which we will try first if a global β0 isn't enough:
   - β0 per taxon (primer mismatch, taxon-specific PCR drift);
   - covariates on the offset: copies (drift is worst at low copy numbers),
     inhibition, dilution;
   - **clumping inside the mechanism** rather than an offset: aliquot copies
     negative binomial instead of Poisson, variance F(1 + F/k), so v_s and the
     zero probability (1 + F/k)^−k change together. It keeps the zeros and
     the replicate spread consistent, which an offset on φ doesn't.
6. **What real replicates look like.** Across 511 bottle × dilution groups
   with at least two replicates, the share of the bottle's most abundant
   non-target taxon (within the non-target, non-human reads) varies between
   replicates a median 570× more than multinomial sampling (middle half
   130–2,550×). If that were all Poisson, it would imply a few tens of
   effective copies per aliquot against about 15,000 reads, well inside the
   range where your README says scale is identified. Human and cattle shares
   vary even more (median 600×), which fits contamination added
   per replicate rather than DNA in the extract.
7. **With few taxa, the others can all be absent.** The zero inflation
   conditions the focal taxon on presence, and `mech_phi_pi` marginalises
   the others. With 15–30 taxa they are never all absent at once. With our
   four categories they often are: an aliquot with no hake or cetacean
   copies gives junk every read. A beta-binomial can't put mass at
   reads = total, so the fit pulls that category's φ down. On simulated
   reads, junk's best β0 at the true copies was −1.0; it was +0.1 with those
   replicates left out. Our fix mirrors the zero inflation on the other side.
   Each taxon's reads are:
   - 0 if it is absent;
   - all the reads if it is present and every other taxon is absent
     (probability e^−F_o, with F_o the others' total copies);
   - otherwise beta-binomial, with the others' mean and variance taken given
     at least one of them is present: mean Σ w F / p_o, and variance
     Σ w² F / p_o − mean² (1 − p_o), where p_o = 1 − e^−F_o.

   When the others always have copies, this reduces exactly to your φ. It
   only matters with few taxa, or in low-concentration samples where every
   taxon but one can drop out. It could matter in your lowest dilution
   rungs.
8. **A departure in the simulator.** `simulate_marker()` generates exactly
   from the mechanism, so it shows that β0 and ε aren't needed when the
   mechanism is true. An option to generate a departure (clumped copies or
   random per-replicate PCR efficiency) would show the other half: how
   biased λ is without the offset terms, and whether they recover it.

## 3. Questions

1. **Composite likelihood.** The likelihood multiplies separate
   zero-inflated beta-binomials over the taxa in an aliquot, although the
   counts must add up to `reads`. Have you checked interval coverage for λ
   in simulation (the README reports coverage of the totals)? Our current
   model does the same, so we'll be checking coverage too.
2. **The mock communities where β0 and ε were needed.** Were the departures
   similar across taxa (pointing to a global or per-taxon offset) or
   specific to sample × taxon? That decides how much structure we give the
   offset.
3. **Cycle numbers.** In your data, do all libraries get the same number of
   PCR cycles, and does the dilution step change it? Ours is unknown so far.
