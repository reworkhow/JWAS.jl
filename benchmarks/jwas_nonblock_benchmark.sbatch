#!/bin/bash
#SBATCH -J jwas_nonblock
#SBATCH -p low
#SBATCH -A publicgrp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=220G
#SBATCH --time=06:00:00
#SBATCH -o /home/qtlcheng/jwas_nonblock_%j.out
#SBATCH -e /home/qtlcheng/jwas_nonblock_%j.err

set -euo pipefail
source /etc/profile >/dev/null 2>&1 || true
module load julia/1.12.1

export BENCH_N=50000
export BENCH_P_LIST=100000,200000
export BENCH_L_LIST=2,4
export BENCH_REPS=2
export WARMUP_N=2000
export WARMUP_P=2000
export WARMUP_L=5
export TARGET_P=2000000
export TARGET_L=2000
export BENCH_TXT=/home/qtlcheng/jwas_nonblock_results_${SLURM_JOB_ID}.txt

julia --project=/home/qtlcheng/Github/JWAS.jl /home/qtlcheng/jwas_nonblock_benchmark.jl
