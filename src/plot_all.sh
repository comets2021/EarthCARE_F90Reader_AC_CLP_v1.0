#!/usr/bin/env bash

# Generate the standard All-period Level 3 PDF figures.

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
LEVEL3_DIR=${LEVEL3_DIR:-"${SCRIPT_DIR}/../output/Level3"}
FIGURE_ROOT=${FIGURE_ROOT:-"${LEVEL3_DIR}/figures/all"}
PYTHON=${PYTHON:-python3.11}

if ! command -v "${PYTHON}" >/dev/null 2>&1; then
    echo "ERROR: Python command not found: ${PYTHON}" >&2
    exit 1
fi

MAP2D_NC="${LEVEL3_DIR}/map2d/ecare_map2d_all.nc"
MAP3D_NC="${LEVEL3_DIR}/map3d/ecare_map3d_all.nc"
CLD3D_NC="${LEVEL3_DIR}/map3d/ecare_cld3d_all.nc"
ISCCP_NC="${LEVEL3_DIR}/isccp/ecare_isccp_all.nc"
CFAD_NC="${LEVEL3_DIR}/cfad/ecare_cfad_all.nc"
C4_NC="${LEVEL3_DIR}/map3d/ecare_c4_3d_all.nc"

for input_file in "${MAP2D_NC}" "${MAP3D_NC}" "${CLD3D_NC}" \
                  "${ISCCP_NC}" "${CFAD_NC}" "${C4_NC}"; do
    if [ ! -f "${input_file}" ] || [ ! -s "${input_file}" ] || [ ! -r "${input_file}" ]; then
        echo "ERROR: Required All-period NetCDF is unavailable: ${input_file}" >&2
        exit 1
    fi
    if [ -n "${LEVEL3_RUN_START_EPOCH:-}" ] && \
       [ "$(stat -c %Y "${input_file}")" -lt "${LEVEL3_RUN_START_EPOCH}" ]; then
        echo "ERROR: All-period NetCDF predates the current Level 3 run: ${input_file}" >&2
        exit 1
    fi
done

run_plot() {
    local diagnostic=$1
    local script=$2
    shift 2
    if [ ! -r "${script}" ]; then
        echo "ERROR: Visualization script is unavailable for ${diagnostic}: ${script}" >&2
        exit 1
    fi
    echo "*** VISUALIZATION: ${diagnostic}"
    if ! "${PYTHON}" "${script}" "$@"; then
        echo "ERROR: Visualization failed for ${diagnostic}: ${script}" >&2
        exit 1
    fi
}

mkdir -p \
    "${FIGURE_ROOT}/map2d" \
    "${FIGURE_ROOT}/map3d" \
    "${FIGURE_ROOT}/cld3d" \
    "${FIGURE_ROOT}/isccp" \
    "${FIGURE_ROOT}/cfad" \
    "${FIGURE_ROOT}/c4_3d" || {
        echo "ERROR: Cannot create figure output directories under ${FIGURE_ROOT}" >&2
        exit 1
    }

run_plot "map2d" "${LEVEL3_DIR}/map2d/plot_map2d.py" \
    "${MAP2D_NC}" -o "${FIGURE_ROOT}/map2d"
run_plot "map3d" "${LEVEL3_DIR}/map3d/plot_zonalmean.py" \
    "${MAP3D_NC}" -o "${FIGURE_ROOT}/map3d"
run_plot "cld3d" "${LEVEL3_DIR}/map3d/plot_zonalmean_cldfrc.py" \
    "${CLD3D_NC}" -o "${FIGURE_ROOT}/cld3d"
run_plot "ISCCP legacy" "${LEVEL3_DIR}/isccp/plot_isccp.py" \
    "${ISCCP_NC}" -o "${FIGURE_ROOT}/isccp" --region glob --diagnostic legacy
run_plot "ISCCP C4" "${LEVEL3_DIR}/isccp/plot_isccp.py" \
    "${ISCCP_NC}" -o "${FIGURE_ROOT}/isccp" --region glob --diagnostic c4
run_plot "CFAD" "${LEVEL3_DIR}/cfad/plot_cfad.py" \
    "${CFAD_NC}" -o "${FIGURE_ROOT}/cfad"
run_plot "C4 3D 500 hPa" "${LEVEL3_DIR}/map3d/plot_c4_3d.py" \
    "${C4_NC}" -o "${FIGURE_ROOT}/c4_3d" \
    --variable c4_cloud_fraction --plev 500
run_plot "C4 3D zonal mean" "${LEVEL3_DIR}/map3d/plot_c4_zonalmean.py" \
    "${C4_NC}" -o "${FIGURE_ROOT}/c4_3d" \
    --variable c4_cloud_fraction

echo "All-period visualization completed successfully: ${FIGURE_ROOT}"
