#!/bin/bash

if [ $# -eq 0 ]; then
    echo "错误：请提供任务名称作为参数"
    echo "用法: $0 <task_name>"
    exit 1
fi


task_name=$1


BASE_DIR=$PWD

OUTPUT_DIR="${BASE_DIR}/DATASET/${task_name}/GEMME"


mkdir -p "${OUTPUT_DIR}"


perl reformat.pl "${BASE_DIR}/MSA/${task_name}.a2m" "${BASE_DIR}/MSA/${task_name}.a3m"


udocker run --rm \
    --volume "${BASE_DIR}:${BASE_DIR}" \
    elodielaine/gemme:gemme \
    sh -c "GEMME_PATH=/opt/GEMME && \
           cd ${OUTPUT_DIR} && \
           python2.7 \$GEMME_PATH/gemme.py \
           ${BASE_DIR}/MSA/${task_name}.a3m \
           -r input \
           -f ${BASE_DIR}/MSA/${task_name}.a3m \
           -m ${BASE_DIR}/DATASET/${task_name}/${task_name}.txt \
           -N 20000"


cd "${OUTPUT_DIR}"

find . -type f ! -name "*normPred_evolCombi.txt" -delete

find . -type d -empty -delete

echo "GEMME processing complete!"