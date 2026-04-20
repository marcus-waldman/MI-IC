#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=15:00
#SBATCH --job-name=miicsem_pilot1
#SBATCH --output=miicsem_pilot1_%j.log
#SBATCH --error=miicsem_pilot1_%j.err
#SBATCH --exclude=csphbiostats.ucdenver.pvt

# =============================================================================
# MIICSEM TIER 1 PILOT
# =============================================================================
#
# Purpose: 60-second smoke test of the miicsem SEM simulation package.
#          Runs 10 replications at N=250, 25% MCAR, writes per-condition
#          and combined .rds files to the results directory specified
#          inside pilot_tier1.R.
#
# Usage: submit from the directory that contains both this .sh and
#        pilot_tier1.R. Log files land in the same directory.
#
#        cd /biostats_share/waldmanm/simulation-studies/MI-IC/SeM
#        sbatch run_pilot_tier1.sh
#
# Expected:
#   - Wall time ~100 seconds on 32 cores
#   - 10/10 replications succeed
# =============================================================================

echo "=============================================="
echo "MIICSEM TIER 1 PILOT"
echo "=============================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo ""

cd "$SLURM_SUBMIT_DIR" || { echo "FATAL: cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }

if [ ! -f pilot_tier1.R ]; then
  echo "FATAL: pilot_tier1.R not found in $SLURM_SUBMIT_DIR"
  exit 1
fi

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
else
  echo "ERROR: Script failed with exit code $EXIT_CODE"
  echo "  Check miicsem_pilot1_${SLURM_JOB_ID}.err"
fi

echo ""
