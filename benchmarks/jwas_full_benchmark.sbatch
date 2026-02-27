#!/bin/bash
#SBATCH -J jwas_full_bench
#SBATCH -p low
#SBATCH -A publicgrp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=220G
#SBATCH --time=04:00:00
#SBATCH -o /home/qtlcheng/jwas_full_bench_%j.out
#SBATCH -e /home/qtlcheng/jwas_full_bench_%j.err

set -euo pipefail
source /etc/profile >/dev/null 2>&1 || true
module load julia/1.12.1

export BENCH_N=50000
export BENCH_P_LIST=100000,200000
export BENCH_L_LIST=223,446
export TARGET_P=2000000
export TARGET_L=2000
export BENCH_JSON=/home/qtlcheng/jwas_full_bench_${SLURM_JOB_ID}.json

julia --project=/home/qtlcheng/Github/JWAS.jl /home/qtlcheng/jwas_full_benchmark.jl
