#!/usr/bin/env bash

set -o pipefail

# ==============================================
# User settings
# Edit the default values below, then run ./main.sh.
# Environment variables with the same names may
# be used for a temporary override if desired.
# ==============================================

# Required: input root containing YYYY/MM directories
#DATA_DIR=${DATA_DIR:-"/path/to/2B.AC_CLP.vCb"}
DATA_DIR=${DATA_DIR:-"/raid_zeus/EarthCARE/Data/CPR+ATLID/2B.AC_CLP.vCb"}

# Required: inclusive analysis period (YYYYMM)
#START_YM=${START_YM:-"202409"}
#END_YM=${END_YM:-"202508"}
START_YM=${START_YM:-"202506"}
END_YM=${END_YM:-"202508"}

# Optional: JAXA input product version
VERSION=${VERSION:-"vCb"}

# Optional: 0 refuses existing NetCDF; 1 permits replacement
#OVERWRITE_OUTPUT=${OVERWRITE_OUTPUT:-0}
OVERWRITE_OUTPUT=${OVERWRITE_OUTPUT:-1}

# Optional: automatic visualization (0: OFF, 1: ON)
#RUN_VISUALIZATION=${RUN_VISUALIZATION:-0}
RUN_VISUALIZATION=${RUN_VISUALIZATION:-1}

# Optional: Python command used only when visualization is ON
PYTHON=${PYTHON:-"python3.11"}

# Internal settings
data=${DATA_DIR}
version=${VERSION}
exe=${ECARE_EXE:-ecare_driver}

# Always run relative-path-dependent commands from this script directory.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
cd "${SCRIPT_DIR}" || exit 1

MAKE=${MAKE:-make}
LEVEL3_DIR="${SCRIPT_DIR}/../output/Level3"
RUN_START_EPOCH=$(date +%s)
RUN_START_TIME=$(date '+%Y-%m-%d %H:%M:%S %z')
RUN_ID="$(date '+%Y%m%dT%H%M%S')_$$"
RUN_DIR="${SCRIPT_DIR}/log/runs/${RUN_ID}"
filelist="${RUN_DIR}/filelist.txt"
summary="${RUN_DIR}/run_summary.txt"

if ! mkdir -p "${RUN_DIR}"; then
    echo "ERROR: Level 3 processing failed." >&2
    exit 1
fi

exec > >(tee -a "${RUN_DIR}/main.stdout") \
     2> >(tee -a "${RUN_DIR}/main.stderr" >&2)

run_logged() {
    local stdout_file=$1
    local stderr_file=$2
    shift 2
    "$@" > >(tee "${stdout_file}") 2> >(tee "${stderr_file}" >&2)
}

valid_ym() {
    local value=$1 year month
    [[ "${value}" =~ ^[0-9]{6}$ ]] || return 1
    year=${value:0:4}
    month=${value:4:2}
    (( 10#${year} >= 1 && 10#${month} >= 1 && 10#${month} <= 12 ))
}

finish_run() {
    local exit_code=$1
    local processing_status=$2
    local visualization_status=$3
    local end_epoch end_time
    end_epoch=$(date +%s)
    end_time=$(date '+%Y-%m-%d %H:%M:%S %z')
    {
        echo "end_time=${end_time}"
        echo "elapsed_seconds=$((end_epoch - RUN_START_EPOCH))"
        echo "processing_status=${processing_status}"
        echo "visualization_status=${visualization_status}"
        echo "exit_code=${exit_code}"
    } >> "${summary}"
    exit "${exit_code}"
}

{
    echo "start_time=${RUN_START_TIME}"
    echo "main_script=${SCRIPT_DIR}/main.sh"
    echo "main_script_sha256=$(sha256sum "${SCRIPT_DIR}/main.sh" | awk '{print $1}')"
    echo "input_data=${data}"
    echo "product_version=${version}"
    echo "start_ym=${START_YM}"
    echo "end_ym=${END_YM}"
    echo "overwrite_output=${OVERWRITE_OUTPUT}"
    echo "run_visualization=${RUN_VISUALIZATION}"
    echo "python=${PYTHON}"
} > "${summary}"

if [ -z "${DATA_DIR}" ] || [ -z "${START_YM}" ] || [ -z "${END_YM}" ]; then
    echo "ERROR: DATA_DIR, START_YM, and END_YM are required." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_configuration" "not_run"
fi
if [ -z "${VERSION}" ]; then
    echo "ERROR: VERSION must not be empty." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_configuration" "not_run"
fi
if ! valid_ym "${START_YM}" || ! valid_ym "${END_YM}"; then
    echo "ERROR: START_YM and END_YM must be valid months in YYYYMM format." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_configuration" "not_run"
fi
if (( 10#${START_YM} > 10#${END_YM} )); then
    echo "ERROR: START_YM must not be later than END_YM." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_configuration" "not_run"
fi
if [ "${RUN_VISUALIZATION}" != "0" ] && [ "${RUN_VISUALIZATION}" != "1" ]; then
    echo "ERROR: RUN_VISUALIZATION must be 0 or 1." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_configuration" "not_run"
fi
if [ "${OVERWRITE_OUTPUT}" != "0" ] && [ "${OVERWRITE_OUTPUT}" != "1" ]; then
    echo "ERROR: OVERWRITE_OUTPUT must be 0 or 1." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_configuration" "not_run"
fi

if [ ! -d "${data}" ]; then
    echo "ERROR: Input data directory does not exist: ${data}" >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_input_directory" "not_run"
fi

for output_dir in map2d map3d isccp cfad; do
    if ! mkdir -p "${LEVEL3_DIR}/${output_dir}" || [ ! -w "${LEVEL3_DIR}/${output_dir}" ]; then
        echo "ERROR: Level 3 output directory is unavailable: ${LEVEL3_DIR}/${output_dir}" >&2
        echo "ERROR: Level 3 processing failed." >&2
        finish_run 1 "failed_output_directory" "not_run"
    fi
done

year_months=()
month_periods=()
annual_periods=()
current_year=$((10#${START_YM:0:4}))
current_month=$((10#${START_YM:4:2}))
last_year=""
while :; do
    year_token=$(printf '%04d' "${current_year}")
    month_token=$(printf '%02d' "${current_month}")
    current_ym="${year_token}${month_token}"
    year_months+=("${year_token}/${month_token}")
    month_periods+=("${current_ym}")
    if [ "${year_token}" != "${last_year}" ]; then
        annual_periods+=("${year_token}")
        last_year=${year_token}
    fi
    [ "${current_ym}" = "${END_YM}" ] && break
    current_month=$((current_month + 1))
    if [ "${current_month}" -eq 13 ]; then
        current_month=1
        current_year=$((current_year + 1))
    fi
done

output_periods=("${month_periods[@]}" "${annual_periods[@]}" "all")
expected_output_files=()
for output_spec in \
    "map2d:ecare_map2d" \
    "map3d:ecare_map3d" \
    "map3d:ecare_cld3d" \
    "isccp:ecare_isccp" \
    "cfad:ecare_cfad" \
    "map3d:ecare_c4_3d"; do
    output_dir=${output_spec%%:*}
    output_prefix=${output_spec#*:}
    for output_period in "${output_periods[@]}"; do
        expected_output_files+=("${LEVEL3_DIR}/${output_dir}/${output_prefix}_${output_period}.nc")
    done
done
echo "month_count=${#month_periods[@]}" >> "${summary}"

if [ "${OVERWRITE_OUTPUT}" -eq 0 ]; then
    existing_outputs=()
    for output_file in "${expected_output_files[@]}"; do
        [ -e "${output_file}" ] && existing_outputs+=("${output_file}")
    done
    if [ "${#existing_outputs[@]}" -gt 0 ]; then
        echo "ERROR: Refusing to overwrite existing Level 3 NetCDF files:" >&2
        printf '  %s\n' "${existing_outputs[@]}" >&2
        echo "Set OVERWRITE_OUTPUT=1 only if replacement is intentional." >&2
        echo "ERROR: Level 3 processing failed." >&2
        finish_run 1 "failed_existing_output" "not_run"
    fi
fi

#-------------------------------
if ! run_logged "${RUN_DIR}/make_clean_start.stdout" \
        "${RUN_DIR}/make_clean_start.stderr" "${MAKE}" clean; then
    echo "ERROR: Initial make clean failed." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_build" "not_run"
fi
if ! run_logged "${RUN_DIR}/make.stdout" "${RUN_DIR}/make.stderr" "${MAKE}"; then
    echo "ERROR: make failed." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_build" "not_run"
fi
if [ ! -x "${SCRIPT_DIR}/${exe}" ]; then
    echo "ERROR: Executable was not created: ${SCRIPT_DIR}/${exe}" >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_build" "not_run"
fi
echo "executable_sha256=$(sha256sum "${SCRIPT_DIR}/${exe}" | awk '{print $1}')" >> "${summary}"

#--- Preprocess for Input Data Preparation:
: > "${filelist}" || {
    echo "ERROR: Cannot create run-specific filelist: ${filelist}" >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_filelist" "not_run"
}

for year_month in "${year_months[@]}"; do
        year=${year_month%/*}
        month=${year_month#*/}
        input_dir="${data}/${year}/${month}"
        echo "================================================"
        echo "             Input Data:  ${year}/${month}"
        echo "================================================"
        if [ ! -d "${input_dir}" ]; then
            echo "ERROR: Input directory does not exist: ${input_dir}" >&2
            echo "ERROR: Level 3 processing failed." >&2
            finish_run 1 "failed_input_directory" "not_run"
        fi
        file_count_before=$(wc -l < "${filelist}")
        if ! find "${input_dir}" -type f -name "*${version}.h5" | sort >> "${filelist}"; then
            echo "ERROR: Filelist generation failed for ${input_dir}" >&2
            echo "ERROR: Level 3 processing failed." >&2
            finish_run 1 "failed_filelist" "not_run"
        fi
        file_count_after=$(wc -l < "${filelist}")
        if [ "${file_count_after}" -eq "${file_count_before}" ]; then
            echo "ERROR: No *${version}.h5 input files found in ${input_dir}" >&2
            echo "ERROR: Level 3 processing failed." >&2
            finish_run 1 "failed_empty_filelist" "not_run"
        fi
done

if [ ! -s "${filelist}" ]; then
    echo "ERROR: Input filelist is empty: ${filelist}" >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_empty_filelist" "not_run"
fi
echo "filelist=${filelist}" >> "${summary}"
echo "filelist_sha256=$(sha256sum "${filelist}" | awk '{print $1}')" >> "${summary}"
echo "file_count=$(wc -l < "${filelist}")" >> "${summary}"

#--- Main Procedures:
echo '*** EXECUTION MAIN'
run_logged "${RUN_DIR}/ecare_driver.stdout" "${RUN_DIR}/ecare_driver.stderr" \
    "${SCRIPT_DIR}/${exe}" "${filelist}" "${version}"
driver_status=$?
echo "ecare_driver_exit_code=${driver_status}" >> "${summary}"
echo "stderr_line_count=$(wc -l < "${RUN_DIR}/ecare_driver.stderr")" >> "${summary}"

#--- Finalize without hiding the driver status:
run_logged "${RUN_DIR}/make_clean_end.stdout" \
    "${RUN_DIR}/make_clean_end.stderr" "${MAKE}" clean
cleanup_status=$?
echo "cleanup_exit_code=${cleanup_status}" >> "${summary}"

if [ "${driver_status}" -ne 0 ]; then
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_ecare_driver" "not_run"
fi
if [ "${cleanup_status}" -ne 0 ]; then
    echo "ERROR: Final make clean failed." >&2
    echo "ERROR: Level 3 processing failed." >&2
    finish_run 1 "failed_cleanup" "not_run"
fi

for output_file in "${expected_output_files[@]}"; do
    if [ ! -s "${output_file}" ] || [ "$(stat -c %Y "${output_file}")" -lt "${RUN_START_EPOCH}" ]; then
        echo "ERROR: Expected Level 3 output was not created by this run: ${output_file}" >&2
        echo "ERROR: Level 3 processing failed." >&2
        finish_run 1 "failed_output_validation" "not_run"
    fi
done

if [ "${RUN_VISUALIZATION}" -eq 1 ]; then
    LEVEL3_RUN_START_EPOCH=${RUN_START_EPOCH} PYTHON=${PYTHON} \
        run_logged "${RUN_DIR}/visualization.stdout" \
        "${RUN_DIR}/visualization.stderr" bash "${SCRIPT_DIR}/plot_all.sh"
    visualization_status=$?
    echo "visualization_exit_code=${visualization_status}" >> "${summary}"
    if [ "${visualization_status}" -ne 0 ]; then
        echo "Level 3 NetCDF generation succeeded, but visualization failed." >&2
        finish_run 2 "success" "failed"
    fi
    echo "Level 3 processing completed successfully."
    finish_run 0 "success" "success"
fi

echo "visualization_exit_code=not_run" >> "${summary}"
echo "Level 3 processing completed successfully."
finish_run 0 "success" "not_run"
