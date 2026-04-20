#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15:00
#SBATCH --job-name=miicsem_pilot1
#SBATCH --output=logs/miicsem_pilot1_%j.log
#SBATCH --error=logs/miicsem_pilot1_%j.err
#SBATCH --exclude=csphbiostats.ucdenver.pvt

# =============================================================================
# MIICSEM TIER 1 PILOT
# =============================================================================
#
# Purpose: 60-second smoke test of the miicsem SEM simulation package.
#          Runs 10 replications at N=250, 25% MCAR, writes per-condition
#          and combined .rds files to:
#            /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results/
#
# Expected:
#   - Wall time ~100 seconds on 32 cores
#   - 10/10 replications succeed
#   - Files: results_n=250_mr=0.25.rds, results_combined.rds
#
# =============================================================================

echo "=============================================="
echo "MIICSEM TIER 1 PILOT"
echo "=============================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo ""

cd /biostats_share/waldmanm/git-repositories/MI-IC/hpc

/home/hillandr/share/R-4.4.1/bin/Rscript pilot_tier1.R

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "EXECUTION COMPLETED"
echo "=============================================="
echo "Job completed at: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
  echo "SUCCESS: Tier 1 pilot complete"
  echo ""
  echo "Results written to:"
  echo "  /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results/"
  echo ""
  echo "Next steps:"
  echo "  1. Inspect results_combined.rds for structure"
  echo "  2. Check chi-square summary for M1 in the .log file"
  echo "  3. If clean, scale to tier 2 (100 reps) or full grid"
else
  echo "ERROR: Script failed with exit code $EXIT_CODE"
  echo "  Check logs/miicsem_pilot1_${SLURM_JOB_ID}.err"
fi

echo ""
