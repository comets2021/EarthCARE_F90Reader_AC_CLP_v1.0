#!/usr/bin/env python3
"""Plot count-pooled latitude-pressure sections from C4 3D diagnostics."""

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


CLASS_OPTIONS = {
    "particle_type": "particle_type",
    "particle_category": "particle_category",
    "type_qc_class": "type_qc",
}

FREQUENCY_DEFINITIONS = {
    "c4_cloud_fraction": ("n_c4_cloud", "n_c4_valid",
                          "C4 Cloud Fraction"),
    "c4_type_occurrence_frequency": ("n_c4_type", "n_c4_valid",
                                     "C4 Type Occurrence Frequency"),
    "c4_type_relative_frequency": ("n_c4_type", "n_c4_type_assigned",
                                   "C4 Type Relative Frequency"),
    "c4_category_occurrence_frequency": ("n_c4_category", "n_c4_valid",
                                         "C4 Category Occurrence Frequency"),
    "c4_category_relative_frequency": ("n_c4_category", "n_c4_category_assigned",
                                       "C4 Category Relative Frequency"),
}

COUNT_VARIABLES = {
    "n_c4_valid", "n_c4_cloud", "n_c4_type", "n_c4_category",
    "n_c4_type_qc", "n_c4_type_unassigned", "n_c4_category_unassigned",
    "n_c4_mask_invalid", "n_c4_category_invalid",
}

COUNT_TITLES = {
    "n_c4_valid": "Zonal-summed C4 Valid Count",
    "n_c4_cloud": "Zonal-summed C4 Cloud Count",
    "n_c4_type": "Zonal-summed C4 Type Count",
    "n_c4_category": "Zonal-summed C4 Category Count",
    "n_c4_type_qc": "Zonal-summed Type QC Count",
    "n_c4_type_unassigned": "Zonal-summed C4 Type-unassigned Count",
    "n_c4_category_unassigned": "Zonal-summed C4 Category-unassigned Count",
    "n_c4_mask_invalid": "Zonal-summed C4 Mask-invalid Count",
    "n_c4_category_invalid": "Zonal-summed Category-invalid Count",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
})


def period_token(nc_path):
    stem = Path(nc_path).stem
    prefix = "ecare_c4_3d_"
    if not stem.startswith(prefix):
        raise ValueError(f"cannot identify period token from input filename: {stem}")
    token = stem[len(prefix):]
    if token != "all" and re.fullmatch(r"\d{4}(?:\d{2})?", token) is None:
        raise ValueError(f"unsupported period token in input filename: {token}")
    return token


def coordinate_index(nc, coordinate, code):
    values = np.asarray(nc.variables[coordinate][:], dtype=int)
    match = np.flatnonzero(values == int(code))
    if match.size != 1:
        available = ", ".join(str(value) for value in values)
        raise ValueError(
            f"code {code} is not present in coordinate '{coordinate}'; "
            f"available codes: {available}"
        )
    return int(match[0])


def class_meaning(nc, coordinate, index, code):
    meanings = getattr(nc.variables[coordinate], "flag_meanings", "").split()
    if index < len(meanings):
        return meanings[index].replace("_", " ")
    return f"code {code}"


def read_count(nc, variable, selections):
    """Read a count field and return it in plev,lat,lon order."""
    if variable not in nc.variables:
        raise KeyError(f"count variable '{variable}' not found")
    var = nc.variables[variable]
    dims = list(var.dimensions)
    required = {"plev", "lat", "lon"}
    if not required.issubset(dims):
        raise ValueError(
            f"count variable '{variable}' must contain plev/lat/lon; got {tuple(dims)}"
        )

    indexer = []
    selected = []
    remaining = []
    for dim in dims:
        if dim in CLASS_OPTIONS:
            code = selections[dim]
            if code is None:
                option = CLASS_OPTIONS[dim].replace("_", "-")
                raise ValueError(f"variable '{variable}' requires --{option}")
            index = coordinate_index(nc, dim, code)
            indexer.append(index)
            selected.append((dim, int(code), class_meaning(nc, dim, index, code)))
        elif dim in required:
            indexer.append(slice(None))
            remaining.append(dim)
        else:
            raise ValueError(f"unsupported dimension '{dim}' in variable '{variable}'")

    data = np.ma.asarray(var[tuple(indexer)], dtype=np.float64)
    fill = getattr(var, "_FillValue", None)
    if fill is not None:
        data = np.ma.masked_where(data == float(fill), data)
    data = np.ma.masked_invalid(data)
    axes = [remaining.index(name) for name in ("plev", "lat", "lon")]
    data = np.ma.transpose(data, axes)
    return data, selected, getattr(var, "long_name", variable)


def validate_class_options(variable, selections, selected):
    used = {dim for dim, _, _ in selected}
    for dim, code in selections.items():
        if code is not None and dim not in used:
            option = CLASS_OPTIONS[dim].replace("_", "-")
            raise ValueError(f"--{option} does not apply to variable '{variable}'")


def zonal_sum(count):
    """Sum native-bin counts over longitude using int64 precision."""
    return np.ma.sum(count, axis=2, dtype=np.int64)


def ratio_of_sums(numerator, denominator):
    """Return 100*sum_lon(numerator)/sum_lon(denominator)."""
    numerator_sum = zonal_sum(numerator)
    denominator_sum = zonal_sum(denominator)
    invalid = np.ma.getmaskarray(denominator_sum) | (denominator_sum <= 0)
    result = np.ma.masked_all(numerator_sum.shape, dtype=np.float64)
    valid = ~invalid
    result[valid] = 100.0 * numerator_sum[valid] / denominator_sum[valid]
    return result, numerator_sum, denominator_sum


def compute_zonal_field(nc, variable, selections):
    """Compute a pooled frequency or a zonal-summed count."""
    if variable in FREQUENCY_DEFINITIONS:
        numerator_name, denominator_name, title = FREQUENCY_DEFINITIONS[variable]
        numerator, selected, _ = read_count(nc, numerator_name, selections)
        denominator, denominator_selected, _ = read_count(nc, denominator_name, selections)
        if denominator_selected:
            raise ValueError(f"unexpected classified denominator '{denominator_name}'")
        validate_class_options(variable, selections, selected)
        field, numerator_sum, denominator_sum = ratio_of_sums(numerator, denominator)
        return field, "%", title, selected, numerator_sum, denominator_sum

    if variable not in COUNT_VARIABLES:
        supported = sorted(FREQUENCY_DEFINITIONS) + sorted(COUNT_VARIABLES)
        raise ValueError(
            f"unsupported variable '{variable}'; supported variables: {', '.join(supported)}"
        )
    count, selected, _ = read_count(nc, variable, selections)
    validate_class_options(variable, selections, selected)
    field = zonal_sum(count)
    title = COUNT_TITLES[variable]
    return field, "1", title, selected, field, None


def plot_limits(units, data, user_vmin, user_vmax):
    valid = data.compressed()
    if not valid.size:
        raise ValueError("the requested field has no valid values")
    vmin = user_vmin if user_vmin is not None else 0.0
    vmax = user_vmax if user_vmax is not None else (
        100.0 if units == "%" else float(valid.max())
    )
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin > vmax:
        raise ValueError(f"invalid color range: vmin={vmin}, vmax={vmax}")
    if vmin == vmax:
        vmax = vmin + 1.0
    return vmin, vmax


def output_stem(args, selected):
    parts = ["c4_3d", period_token(args.nc), args.variable]
    for dim, code, _ in selected:
        label = {"particle_type": "type", "particle_category": "category",
                 "type_qc_class": "typeqc"}[dim]
        parts.append(f"{label}{code}")
    parts.append("zonalmean")
    return "_".join(parts)


def draw(args):
    selections = {
        "particle_type": args.particle_type,
        "particle_category": args.particle_category,
        "type_qc_class": args.type_qc,
    }
    with Dataset(args.nc, "r") as nc:
        for coordinate in ("lat", "lon", "plev"):
            if coordinate not in nc.variables:
                raise KeyError(f"coordinate '{coordinate}' not found")
        lat = np.asarray(nc.variables["lat"][:], dtype=float)
        plev = np.asarray(nc.variables["plev"][:], dtype=float)
        field, units, title, selected, _, _ = compute_zonal_field(
            nc, args.variable, selections
        )

    vmin, vmax = plot_limits(units, field, args.vmin, args.vmax)
    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad((0.0, 0.0, 0.0, 0.0))
    valid = field.compressed()
    extend = "max" if valid.size and float(valid.max()) > vmax else "neither"

    fig = plt.figure(figsize=(8.0, 4.8))
    ax = fig.add_subplot(111)
    mesh = ax.pcolormesh(lat, plev, field, shading="auto", cmap=cmap,
                         vmin=vmin, vmax=vmax)
    colorbar = fig.colorbar(mesh, ax=ax, extend=extend)
    class_text = ""
    if selected:
        dim, code, meaning = selected[0]
        label = {"particle_type": "Type", "particle_category": "Category",
                 "type_qc_class": "QC"}[dim]
        class_text = f" ({label} {code}: {meaning})"
    colorbar_label = "Pooled zonal frequency [%]" if units == "%" else "Zonal-summed count [1]"
    colorbar.set_label(colorbar_label, fontsize=12)
    colorbar.ax.tick_params(labelsize=11)

    ax.set_xlim(-90, 90)
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xticklabels(["90S", "60S", "30S", "Eq", "30N", "60N", "90N"])
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=12)
    ax.set_ylabel("Pressure (hPa)", fontsize=13)
    ax.set_yscale("linear")
    ax.set_ylim(1000, 0)
    ax.set_yticks(np.arange(0, 1001, 200))
    ax.tick_params(axis="y", labelsize=12)
    ax.set_title(f"{title}{class_text}", fontsize=13)
    fig.tight_layout()

    os.makedirs(args.outdir, exist_ok=True)
    stem = Path(args.outdir) / output_stem(args, selected)
    outputs = []
    if args.write_pdf:
        output = Path(f"{stem}.pdf")
        fig.savefig(output, dpi=args.dpi)
        outputs.append(output)
    if args.write_png:
        output = Path(f"{stem}.png")
        fig.savefig(output, dpi=args.dpi)
        outputs.append(output)
    plt.close(fig)
    print(f"[OK] {args.variable} -> {', '.join(str(path) for path in outputs)}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot count-pooled latitude-pressure sections from an EarthCARE C4 3D NetCDF. "
            "Frequencies are ratios of longitude-summed counts, not means of percentages."
        )
    )
    parser.add_argument("nc", help="Input NetCDF file (ecare_c4_3d_<period>.nc)")
    parser.add_argument("-o", "--outdir", default="zonal_mean_figs",
                        help="Output directory (default: zonal_mean_figs)")
    parser.add_argument("--variable", default="c4_cloud_fraction",
                        help="C4 frequency or count variable")
    parser.add_argument("--particle-type", type=int,
                        help="Official particle_type coordinate code (1-14)")
    parser.add_argument("--particle-category", type=int,
                        help="Official particle_category coordinate code (1-12)")
    parser.add_argument("--type-qc", type=int,
                        help="Official type_qc_class coordinate code (15,16,17,-999)")
    parser.add_argument("--cmap", default="turbo",
                        help="Matplotlib colormap (default: turbo)")
    parser.add_argument("--vmin", type=float, help="Colorbar minimum")
    parser.add_argument("--vmax", type=float, help="Colorbar maximum")
    parser.add_argument("--format", choices=("png", "pdf"), default=None,
                        help="Deprecated compatibility format: png-only or pdf-only")
    parser.add_argument("--png", action="store_true",
                        help="Also save PNG files (PDF is the standard output)")
    parser.add_argument("--dpi", type=int, default=200,
                        help="Output resolution in dpi (default: 200)")
    args = parser.parse_args()

    if not os.path.isfile(args.nc):
        parser.error(f"input NetCDF not found: {args.nc}")
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    if args.format == "png" and args.png:
        parser.error("--format png and --png cannot be used together")
    period_token(args.nc)
    args.write_pdf = args.format != "png"
    args.write_png = args.png or args.format == "png"
    draw(args)


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
