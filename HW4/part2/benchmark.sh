DATA_PATH="/data/pp/hw4/test"
DATASETS=(
    "${DATA_PATH}/data1_1"
    "${DATA_PATH}/data1_2"
    "${DATA_PATH}/data1_3"
    "${DATA_PATH}/data1_4"
    "${DATA_PATH}/data1_5"
    "${DATA_PATH}/data1_6"
    "${DATA_PATH}/data1_7"
    "${DATA_PATH}/data1_8"
    "${DATA_PATH}/data1_9"
    "${DATA_PATH}/data1_10"
    "${DATA_PATH}/data1_11"
    "${DATA_PATH}/data1_12"
    "${DATA_PATH}/data1_13"
    "${DATA_PATH}/data1_14"
    "${DATA_PATH}/data1_15"
    "${DATA_PATH}/data1_16"
    "${DATA_PATH}/data1_17"
    "${DATA_PATH}/data1_18"
    "${DATA_PATH}/data1_19"
    "${DATA_PATH}/data1_20"
    "${DATA_PATH}/data2_1"
    "${DATA_PATH}/data2_2"
    "${DATA_PATH}/data2_3"
    "${DATA_PATH}/data2_4"
    "${DATA_PATH}/data2_5"
    "${DATA_PATH}/data2_6"
    "${DATA_PATH}/data2_7"
    "${DATA_PATH}/data2_8"
    "${DATA_PATH}/data2_9"
    "${DATA_PATH}/data2_10"
    "${DATA_PATH}/data2_11"
    "${DATA_PATH}/data2_12"
    "${DATA_PATH}/data2_13"
    "${DATA_PATH}/data2_14"
    "${DATA_PATH}/data2_15"
    "${DATA_PATH}/data2_16"
    "${DATA_PATH}/data2_17"
    "${DATA_PATH}/data2_18"
    "${DATA_PATH}/data2_19"
    "${DATA_PATH}/data2_20"
)

echo "DATASET,TIME_SECONDS"

make clean >/dev/null 2>&1
make TILE_SIZE=${tile} >/dev/null 2>&1

for data in "${DATASETS[@]}"
do
    dataname=$(basename ${data})

    time_sec=$(run --mpi=pmix -N 4 -n 4 -- ./matmul ${data} | tail -n 1 | awk '{print $4}')

    echo "${dataname},${time_sec}"
done