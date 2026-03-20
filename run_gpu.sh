#!/bin/bash
# ================================================================
# run_gpu.sh — 自動尋找空閒 GPU 並啟動 FluidX3D
# 用法: bash run_gpu.sh [GPU數量]    (預設=4)
# ================================================================
cd ~/Reference7_FluidX3D

NEED=${1:-4}  # 預設需要 4 顆 GPU

# 找出記憶體使用量 < 100 MiB 的 GPU（只有 Xorg ~4MiB 的算空閒）
FREE_GPUS=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits \
    | awk -F', ' '$2 < 100 {print $1}' \
    | head -n "$NEED" \
    | tr '\n' ' ')

NUM_FREE=$(echo $FREE_GPUS | wc -w)

if [ "$NUM_FREE" -lt "$NEED" ]; then
    echo "=========================================="
    echo " ERROR: 只找到 ${NUM_FREE} 顆空閒 GPU，需要 ${NEED} 顆"
    echo " 空閒 GPU: ${FREE_GPUS:-無}"
    echo "=========================================="
    echo ""
    echo "目前 GPU 使用狀況:"
    nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv,noheader
    exit 1
fi

echo "=========================================="
echo " 找到空閒 GPU: ${FREE_GPUS}"
echo " 啟動 FluidX3D (${NEED} GPU)..."
echo "=========================================="

nohup bin/FluidX3D $FREE_GPUS > "log$(date +%Y%m%d)" 2>&1 &
echo " PID: $!"
echo " Log: log$(date +%Y%m%d)"
