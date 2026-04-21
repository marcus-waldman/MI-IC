#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=100G
#SBATCH --time=30:00:00
#SBATCH --job-name=miicsem_fullgrid
#SBATCH --output=miicsem_fullgrid_%j.log
#SBATCH --error=miicsem_fullgrid_%j.err
#SBATCH --exclude=csphbiostats.ucdenver.pvt

# =============================================================================
# MIICSEM FULL-GRID STUDY 2 SIMULATION
# =============================================================================
#
# Purpose: Run the full Study 2 SEM simulation on HPC.
#          2000 replications x 12 conditions
#          (N = 100, 250, 500, 1000 x missing = 10%, 25%, 40% MCAR)
#          Writes per-condition and combined .rds files to:
#            /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-full/
#
# Usage: submit from the directory that contains both this .sh and
#        full_grid.R. Log/err files land in the same directory.
#
#        cd /biostats_share/waldmanm/simulation-studies/MI-IC/SeM
#        sbatch run_full_grid.sh
#
# Expected: wall time ~18-20 hours at M=100 (5x the M=20 run).
# Wall-time limit here is 30 hours — generous headroom for slow reps.
# =============================================================================

echo "=============================================="
echo "MIICSEM FULL-GRID STUDY 2 SIMULATION"
echo "=============================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo "Cores allocated (SLURM_CPUS_PER_TASK): $SLURM_CPUS_PER_TASK"
echo ""

cd "$SLURM_SUBMIT_DIR" || { echo "FATAL: cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }

if [ ! -f full_grid.R ]; then
  echo "FATAL: full_grid.R not found in $SLURM_SUBMIT_DIR"
  exit 1
fi

/home/hillandr/share/R-4.4.1/bin/Rscript full_grid.R

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "EXECUTION COMPLETED"
echo "=============================================="
echo "Job completed at: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
  echo "SUCCESS: Full-grid simulation complete"
  echo ""
  echo "Results at:"
  echo "  /biostats_share/waldmanm/simulation-studies/MI-IC/SeM/results-full/"
  echo ""
  echo "Next steps:"
  echo "  1. Inspect results_combined.rds"
  echo "  2. Run load_chi_squares() and selection_accuracy_table() locally"
  echo "  3. Begin Study 2 writeup in manuscript Section 4"
else
  echo "ERROR: Script failed with exit code $EXIT_CODE"
  echo "  Check miicsem_fullgrid_${SLURM_JOB_ID}.err"
  echo ""
  echo "Note: run_simulation() checkpoints each condition separately."
  echo "      Resubmit to resume from the last completed condition."
fi

echo ""
