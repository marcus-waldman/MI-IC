#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=100G
#SBATCH --time=03:00:00
#SBATCH --job-name=miicsem_decomp
#SBATCH --output=miicsem_decomp_%j.log
#SBATCH --error=miicsem_decomp_%j.err
#SBATCH --exclude=csphbiostats.ucdenver.pvt

# =============================================================================
# MIICSEM THREE-TERM DECOMPOSITION RUN
# =============================================================================
#
# Purpose: Targeted run to identify which of the three bias-decomposition
#          terms has the wrong sign.
#            N         in {100, 250, 1000}
#            miss_rate = 0.40
#            M         = 50
#            n_reps    = 500
#
#          Writes per-condition and combined .rds files to:
#            /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-decomp/
#
# Usage: submit from the directory that contains both this .sh and
#        decomp_run.R. Log/err files land in the same directory.
#
#        cd /biostats_share/waldmanm/simulation-studies/MI-IC/SeM
#        sbatch run_decomp.sh
#
# Expected: wall time ~45–75 minutes (3 conditions x 500 reps at M=50 on
# 100 cores, with the added overhead of 2 extra complete-data log-lik
# evaluations per model per rep).  3-hour limit gives headroom.
# =============================================================================

echo "=============================================="
echo "MIICSEM THREE-TERM DECOMPOSITION RUN"
echo "=============================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo "Cores allocated (SLURM_CPUS_PER_TASK): $SLURM_CPUS_PER_TASK"
echo ""

cd "$SLURM_SUBMIT_DIR" || { echo "FATAL: cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }

if [ ! -f decomp_run.R ]; then
  echo "FATAL: decomp_run.R not found in $SLURM_SUBMIT_DIR"
  exit 1
fi

/home/hillandr/share/R-4.4.1/bin/Rscript decomp_run.R

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "EXECUTION COMPLETED"
echo "=============================================="
echo "Job completed at: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
  echo "SUCCESS: Decomposition simulation complete"
  echo ""
  echo "Results at:"
  echo "  /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-decomp/"
  echo ""
  echo "Next steps:"
  echo "  1. Transfer results-decomp/ to local machine"
  echo "  2. Run hpc/analyze_decomp.R to tabulate Term 1 / 2 / 3 biases"
  echo "     vs theoretical predictions across all conditions"
else
  echo "ERROR: Script failed with exit code $EXIT_CODE"
  echo "  Check miicsem_decomp_${SLURM_JOB_ID}.err"
  echo ""
  echo "Note: run_simulation() checkpoints each condition separately."
  echo "      Resubmit to resume from the last completed condition."
fi

echo ""
