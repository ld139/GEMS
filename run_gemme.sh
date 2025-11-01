#!/bin/bash

# 检查是否提供了任务名称参数
if [ $# -eq 0 ]; then
    echo "错误：请提供任务名称作为参数"
    echo "用法: $0 <task_name>"
    exit 1
fi

# 获取任务名称参数
task_name=$1

# 设置基础目录
BASE_DIR=$PWD

# 设置输出目录
OUTPUT_DIR="${BASE_DIR}/DATASET/${task_name}/GEMME"

# 创建输出目录（如果不存在）
mkdir -p "${OUTPUT_DIR}"

# 转换 a2m 到 a3m
perl reformat.pl "${BASE_DIR}/MSA/${task_name}.a2m" "${BASE_DIR}/MSA/${task_name}.a3m"

# 在 Docker 容器中执行 GEMME 命令
docker run -it --rm \
    --volume "${BASE_DIR}:${BASE_DIR}" \
    99b37549eec8 \
    sh -c "GEMME_PATH=/opt/GEMME && \
           cd ${OUTPUT_DIR} && \
           python2.7 \$GEMME_PATH/gemme.py \
           ${BASE_DIR}/MSA/${task_name}.a3m \
           -r input \
           -f ${BASE_DIR}/MSA/${task_name}.a3m \
           -m ${BASE_DIR}/DATASET/${task_name}/${task_name}.txt \
           -N 20000"

# 清理输出目录，只保留关键文件
cd "${OUTPUT_DIR}"
# 删除除关键文件外的所有文件
find . -type f ! -name "${task_name}_normPred_evolCombi.txt" -delete
# 删除空目录
find . -type d -empty -delete

# echo "GEMME 处理完成！结果保存在: ${OUTPUT_DIR}/${task_name}_normPred_evolCombi.txt"