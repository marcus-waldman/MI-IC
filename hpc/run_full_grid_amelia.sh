#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=100G
#SBATCH --time=16:00:00
#SBATCH --job-name=miicsem_amelia
#SBATCH --output=miicsem_amelia_%j.log
#SBATCH --error=miicsem_amelia_%j.err
#SBATCH --exclude=csphbiostats.ucdenver.pvt

# =============================================================================
# MIICSEM FULL-GRID STUDY 2 SIMULATION (AMELIA BACKEND)
# =============================================================================
#
# Purpose: amelia x empri sweep — congeniality dose-response.
#          Same N x mr grid run twice with two empri levels:
#            empri = 0.1 * N   (mostly congenial — modest stabilizer)
#            empri = 100 * N   (severely uncongenial — diagonal prior
#                               drowns out CFA structure)
#
#          2000 replications x 6 conditions x 2 empri = 24 cells
#          (N = 100, 500, 1000 x missing = 20%, 40% MCAR x empri sweep)
#          M = 100 imputations
#          Imputer: Amelia (joint MVN, EMB)
#
#          Writes per-condition and combined .rds files to:
#            /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-full-M100-amelia/
#
# Usage: submit from the directory that contains both this .sh and
#        full_grid_amelia.R.  Log/err files land in the same directory.
#
#        cd /biostats_share/waldmanm/simulation-studies/MI-IC/SeM
#        sbatch run_full_grid_amelia.sh
#
# Expected wall time: ~10-13 hours on 100 cores (24 cells; back to the
# original scale of the PMM grid since we doubled by adding empri).
# =============================================================================

echo "=============================================="
echo "MIICSEM FULL-GRID AMELIA SIMULATION"
echo "=============================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo "Cores allocated (SLURM_CPUS_PER_TASK): $SLURM_CPUS_PER_TASK"
echo ""

cd "$SLURM_SUBMIT_DIR" || { echo "FATAL: cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }

if [ ! -f full_grid_amelia.R ]; then
  echo "FATAL: full_grid_amelia.R not found in $SLURM_SUBMIT_DIR"
  exit 1
fi

/home/hillandr/share/R-4.4.1/bin/Rscript full_grid_amelia.R

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "EXECUTION COMPLETED"
echo "=============================================="
echo "Job completed at: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
  echo "SUCCESS: Full-grid amelia simulation complete"
  echo ""
  echo "Results at:"
  echo "  /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-full-M100-amelia/"
  echo ""
  echo "Next steps:"
  echo "  1. Transfer results-full-M100-amelia/ to local"
  echo "  2. Run analyze comparison vs results-full-M100/ (PMM)"
  echo "  3. Confirm r_term1 ~ +1 across all 12 conditions"
else
  echo "ERROR: Script failed with exit code $EXIT_CODE"
  echo "  Check miicsem_amelia_${SLURM_JOB_ID}.err"
  echo ""
  echo "Note: run_simulation() checkpoints each condition separately."
  echo "      Resubmit to resume from the last completed condition."
fi

echo ""
