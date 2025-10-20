#!/bin/bash --login
#SBATCH --job-name=af2_pmpnn_array
#SBATCH --partition=ll
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=5          # 每个任务分配 5 个 CPU
#SBATCH --time=12-00:00:00
#SBATCH --array=1-2%10              # 1-10 个任务，同时最多 10 个并行,根据输入文件的多少修改numnber%10
#SBATCH --mail-user=hengbin.gao@uwa.edu.au
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/af2_pmpnn_%A_%a.out
#SBATCH --error=logs/af2_pmpnn_%A_%a.err

# ========================================
# 初始化环境
# ========================================
source activate
conda activate /group/ll010/hgao/conda_env/scFv

cd /group/ll010/hgao/chromatin_remodel/scFv/histone_train

# ========================================
# 输入目录列表（根据你的 fasta 文件夹修改）
# ========================================
FASTA_DIRS=(
"/group/ll010/hgao/chromatin_remodel/scFv/histone_train/Flag"
"/group/ll010/hgao/chromatin_remodel/scFv/histone_train/H3K27Ac"
)

# ========================================
# 根据数组索引选择目标目录
# ========================================
TARGET_DIR=${FASTA_DIRS[$SLURM_ARRAY_TASK_ID-1]}

if [ -z "$TARGET_DIR" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID is out of range!"
    exit 1
fi

echo "[$(date)] Starting AF2+PMPNN for: $TARGET_DIR"
echo "Job ID: $SLURM_JOB_ID  |  Task ID: $SLURM_ARRAY_TASK_ID"
echo "Using CPUs: $SLURM_CPUS_PER_TASK"

# ========================================
# 运行你的主程序
# ========================================
bash af2_pmpnn_parallel.sh "$TARGET_DIR"

echo "[$(date)] Finished processing $TARGET_DIR"
