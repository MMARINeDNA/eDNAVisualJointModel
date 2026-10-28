#!/usr/bin/env bash
# =============================================================================
# stats3_reruns.sh
#
# Runs the 35N re-runs (ROADMAP_2026-10-06, 6.0; ROADMAP_2026-10-21, section 3)
# as jobs in one tmux session, one window
# per job. Every job runs under nice from the repo root (HSGP4eDNA jobs from
# the submodule) and logs to logs/<job>.log. See notes/stats3_guide.md.
#
#   tools/stats3_reruns.sh list            stages, jobs and cores
#   tools/stats3_reruns.sh stage <1-5>     launch every job in a stage
#   tools/stats3_reruns.sh launch <job>    launch one job
#   tools/stats3_reruns.sh run <job>       run one job in this shell
#   tools/stats3_reruns.sh status          last line of every log
#
# Env: SESSION (tmux session name, default "reruns"), NICE (default 10),
#      CMP_SEEDS (joint-comparison seeds for stage 3, default "2 3"; e.g.
#      CMP_SEEDS="4 5" tools/stats3_reruns.sh stage 3 for a second wave),
#      CMP_SEA (1: stage 3 runs the sea-state comparison instead: eDNA-only,
#      visual-only, joint and joint_conv on simulate_joint(sea_state = TRUE),
#      fitted with the Beaufort covariate and, for contrast with CMP_NOBFT=1
#      (default), visual-only, joint and joint_conv without it; 16-28 cores
#      per seed), CMP_NOBFT (0: skip the without-covariate contrast, e.g. for
#      seeds 3-5). R/compare_joint.r's defaults apply: NB counts for both
#      cetaceans, non-centred junk noise, adapt_delta 0.9.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
SESSION=${SESSION:-reruns}
NICE=${NICE:-10}
CMP_SEEDS=${CMP_SEEDS:-2 3}
CMP_SEA=${CMP_SEA:-0}
CMP_NOBFT=${CMP_NOBFT:-1}
CMP_MODES="visual edna joint joint_conv"
SEA_MODES="edna visual joint joint_conv"   # the eDNA-only fit ignores sea state, but is on the same data
NOBFT_MODES="visual joint joint_conv"
mkdir -p logs

# Stage 1 must finish first: it rebuilds the real data and compiles the models.
# Stages 2 to 4 only need stage 1; stage 5 needs 3 and 4. Stage 3 is the joint
# comparison: one job per seed in CMP_SEEDS and mode, 4 chains each (16 cores
# per seed; the visual-only job ends after about 30 min).
# expl_edna_real is in no stage: its data are in, but its script
# (docs/explainer/04_fit_edna_real.R) is still the earlier metabarcoding-only
# fit, to be rewritten first (ROADMAP_2026-10-21, section 5 item 7).
STAGE1="prep"
STAGE2="hsgp_nobathy hsgp_bathysp_surface hsgp_bathysp_depthspref hsgp_bathygp_surface
        visual_val visual_val_sea joint_sim"
if [ "$CMP_SEA" = 1 ]; then
  STAGE3=$(for s in $CMP_SEEDS; do
             for m in $SEA_MODES; do printf 'cmpsea_s%s_%s ' "$s" "$m"; done
             if [ "$CMP_NOBFT" = 1 ]; then for m in $NOBFT_MODES; do printf 'cmpseanobft_s%s_%s ' "$s" "$m"; done; fi
           done)
else
  STAGE3=$(for s in $CMP_SEEDS; do for m in $CMP_MODES; do printf 'cmp_s%s_%s ' "$s" "$m"; done; done)
fi
STAGE4="expl_edna_sim expl_visual_sim expl_real_humpback expl_real_pwsd
        expl_real1600_humpback expl_real1600_pwsd"
STAGE5="finish"
# Stage 6: the joint model on the real data (R/fit_joint_real.r), one job per
# mode, 4 cores each; needs only 'compile' (the real objects are tracked)
STAGE6="jreal_edna jreal_visual jreal_joint_conv"

hsgp() { (cd external/HSGP4eDNA && Rscript "$@"); }

# The processed objects were regenerated on the 35N domain on main (#85), so the
# loader should reproduce them; report any difference rather than overwrite.
check_processed() {
  for f in muri_edna muri_visual; do git show "HEAD:data/processed/$f.rds" > "logs/$f.committed.rds"; done
  Rscript -e 'for (f in c("muri_edna", "muri_visual")) {
    r <- all.equal(readRDS(sprintf("logs/%s.committed.rds", f)), readRDS(sprintf("data/processed/%s.rds", f)))
    cat(sprintf("data/processed/%s.rds: %s\n", f,
        if (isTRUE(r)) "matches the committed copy" else paste0("DIFFERS from the committed copy: ", paste(head(r, 5), collapse = "; "))))
  }'
}

job() {
  case "$1" in
    # ---- stage 1: real data, tests, compile every model (about 1 h, 2 cores)
    prep)
      for f in data/qPCR_2019A.rds data/M3_compiled_taxon_table_wide.csv data/Hake_2019_metadata.csv \
               data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv; do
        [ -f "$f" ] || { echo "missing $f: copy it from the Mac (notes/stats3_guide.md, Part 4)"; exit 1; }
      done
      Rscript R/01_load_muri.r
      check_processed
      Rscript R/02_check_muri.r
      Rscript tests/smoke.R
      Rscript tests/visual_detection.R
      hsgp -e 'cmdstanr::cmdstan_model("stan/hsgp_nd.stan", include_paths = "stan/include")' ;;
    # ---- stage 6: the joint model on the real data (4 cores each)
    #      jreal_<mode>_ctcond: the Ct SD bounded at low copies (CT_SD=conditional),
    #      into outputs/joint/real_ctcond, keeping the chain CSVs there
    jreal_*)
      [[ $1 =~ ^jreal_(edna|visual|joint|joint_conv)(_ctcond)?$ ]] ||
        { echo "bad job name $1 (jreal_<edna|visual|joint|joint_conv>[_ctcond])"; exit 2; }
      if [ -n "${BASH_REMATCH[2]}" ]; then
        MODES=${BASH_REMATCH[1]} CT_SD=conditional OUT=outputs/joint/real_ctcond \
          CSV_DIR=outputs/joint/real_ctcond/csv_${BASH_REMATCH[1]} Rscript R/fit_joint_real.r
      else
        MODES=${BASH_REMATCH[1]} Rscript R/fit_joint_real.r
      fi ;;
    # compile the visual and joint models once (no data needed), before
    # launching jobs that would otherwise all compile them at the same time
    compile)
      Rscript -e 'source("R/fit_joint.r"); source("R/fit_visual.r"); invisible(compile_joint()); invisible(compile_visual()); cat("compiled\n")' ;;
    # ---- stage 2: HSGP4eDNA 2-D validation (3 cores each)
    hsgp_nobathy)
      hsgp R/01_sim_nobathy_surface.r
      hsgp R/03_format_stan_data_gp2d.r
      hsgp R/04_fit_gp2d.r
      hsgp R/05_check_gp2d.r ;;
    hsgp_bathysp_surface|hsgp_bathysp_depthspref|hsgp_bathygp_surface)
      local sc=${1#hsgp_}
      hsgp "R/01_sim_${sc}.r"
      SCENARIO=$sc hsgp R/09_fit_gp2d_bathysp.r ;;
    # ---- stage 2: visual validation (4 cores each), joint simulator (3 cores)
    visual_val)     SEEDS=1:5 Rscript R/validate_visual.r ;;
    visual_val_sea) SEA_STATE=1 SEEDS=1:3 Rscript R/validate_visual.r ;;
    joint_sim)      Rscript R/validate_joint_sim.r ;;
    # ---- stage 3: joint comparison, one job per seed and mode (4 cores each)
    cmp_s*)
      [[ $1 =~ ^cmp_s([0-9]+)_(edna|visual|joint|joint_conv)$ ]] ||
        { echo "bad job name $1 (cmp_s<seed>_<edna|visual|joint|joint_conv>)"; exit 2; }
      SEED=${BASH_REMATCH[1]} MODES=${BASH_REMATCH[2]} Rscript R/compare_joint.r ;;
    # ---- stage 3 with CMP_SEA=1: sea-state comparison, with and without the covariate
    cmpsea_s*|cmpseanobft_s*)
      [[ $1 =~ ^cmpsea(nobft)?_s([0-9]+)_(edna|visual|joint|joint_conv)$ ]] ||
        { echo "bad job name $1 (cmpsea[nobft]_s<seed>_<edna|visual|joint|joint_conv>)"; exit 2; }
      local nobft=${BASH_REMATCH[1]} seed=${BASH_REMATCH[2]} mode=${BASH_REMATCH[3]}
      SEA_STATE=1 SEED=$seed MODES=$mode BFT_SPECIES=${nobft:+none} Rscript R/compare_joint.r ;;
    # ---- stage 4: explainer fits (2-3 cores each; expl_edna_real waits for its script rewrite)
    expl_edna_sim)
      Rscript docs/explainer/01_fit_edna_sim.R
      Rscript docs/explainer/05_fig_simulation.R
      Rscript docs/explainer/06_fig_edna_sim.R ;;
    expl_visual_sim) Rscript docs/explainer/02_fit_visual_sim.R ;;
    expl_real_humpback|expl_real_pwsd)
      SP=${1#expl_real_} Rscript docs/explainer/03_fit_visual_real.R ;;
    expl_real1600_humpback|expl_real1600_pwsd)
      SP=${1#expl_real1600_} M_MAX=1600 TAG=real1600 Rscript docs/explainer/03_fit_visual_real.R ;;
    expl_edna_real)
      Rscript docs/explainer/04_fit_edna_real.R
      Rscript docs/explainer/09_fig_edna_real.R ;;
    # ---- stage 5: comparison tables and the remaining figures (minutes)
    finish)
      # tables and figure for every seed with comparison fits (with sea state:
      # compare_sea_s<seed>, compare_sea_nobft_s<seed>)
      for d in outputs/joint/compare_s[0-9]* outputs/joint/compare_sea_s[0-9]* \
               outputs/joint/compare_sea_nobft_s[0-9]*; do
        [ -d "$d" ] && ls "$d"/fit_*.rds > /dev/null 2>&1 || continue
        b=$(basename "$d"); seed=${b##*_s}
        case "$b" in
          compare_sea_nobft_s*) SEA_STATE=1 BFT_SPECIES=none SEED=$seed MODES= Rscript R/compare_joint.r ;;
          compare_sea_s*)       SEA_STATE=1 SEED=$seed MODES= Rscript R/compare_joint.r ;;
          *)                    SEED=$seed MODES= Rscript R/compare_joint.r ;;
        esac
        cp "$d/compare_joint.png" "docs/figures/joint_${b}.png"
      done
      Rscript docs/explainer/07_fig_visual.R
      Rscript docs/explainer/08_fig_key_choice.R
      Rscript docs/explainer/10_fig_sea_state.R ;;
    *) echo "unknown job: $1 (try: $0 list)" >&2; return 2 ;;
  esac
}

stage_jobs() {
  case "$1" in
    1) echo $STAGE1 ;; 2) echo $STAGE2 ;; 3) echo $STAGE3 ;;
    4) echo $STAGE4 ;; 5) echo $STAGE5 ;; 6) echo $STAGE6 ;;
    *) echo "unknown stage: $1" >&2; return 2 ;;
  esac
}

launch() {
  tmux has-session -t "$SESSION" 2>/dev/null || tmux new-session -d -s "$SESSION" -n home
  tmux new-window -d -t "$SESSION" -n "$1" "$ROOT/tools/stats3_reruns.sh run $1; exec bash"
  echo "launched $1 (tmux window $SESSION:$1, log logs/$1.log)"
}

case "${1:-list}" in
  run)
    renice -n "$NICE" $$ >/dev/null
    log="logs/$2.log"
    echo "[$(date '+%F %T')] start $2 on $(hostname)" | tee "$log"
    # a fresh process, so set -e stops the job at its first failing step
    if "$ROOT/tools/stats3_reruns.sh" _job "$2" 2>&1 | tee -a "$log"; [ "${PIPESTATUS[0]}" -eq 0 ]; then
      echo "[$(date '+%F %T')] DONE $2" | tee -a "$log"
    else
      echo "[$(date '+%F %T')] FAILED $2" | tee -a "$log"; exit 1
    fi ;;
  _job)   job "$2" ;;
  launch) launch "$2" ;;
  stage)  for j in $(stage_jobs "$2"); do launch "$j"; done ;;
  status)
    for f in logs/*.log; do
      [ -e "$f" ] || { echo "no logs yet"; break; }
      printf '%-26s %s\n' "$(basename "$f" .log)" "$(tail -n 1 "$f" | cut -c1-90)"
    done ;;
  list)
    cat <<EOF
stage 1  (run alone first; ~1 h)                 $STAGE1
         (or, for simulation-only jobs such as stage 3: 'run compile', a few minutes)
stage 2  (~23 cores; hours to ~2 days)           $(echo $STAGE2)
stage 3  (16 cores per seed; ~4-6 h per job)     $(echo $STAGE3)
         seeds from CMP_SEEDS (now "$CMP_SEEDS"); CMP_SEA=1 for the sea-state
         comparison (now CMP_SEA=$CMP_SEA; 24 cores per seed)
stage 4  (~14 cores; hours)                      $(echo $STAGE4)
stage 5  (after 3 and 4; minutes)                $STAGE5
stage 6  (real data; 12 cores; visual ~20 h, eDNA and joint ~37 h on a busy stats3)
                                                 $(echo $STAGE6)
         launch jreal_<mode>_ctcond for the Ct SD bounded at low copies
not in a stage (needs the ROADMAP_2026-10-21 section 5 data): expl_edna_real
EOF
    ;;
  *) echo "usage: $0 list | stage <1-5> | launch <job> | run <job> | status" >&2; exit 2 ;;
esac
