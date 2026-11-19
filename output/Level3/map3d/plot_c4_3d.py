#!/usr/bin/env python3
"""Plot a pressure-level map from the Work 38/39 C4 3D NetCDF."""

import argparse
import os
import re
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from netCDF4 import Dataset


CLASS_OPTIONS = {
    "particle_type": "particle_type",
    "particle_category": "particle_category",
    "type_qc_class": "type_qc",
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


def read_level_map(nc, variable, pressure, selections):
    if variable not in nc.variables:
        raise KeyError(f"variable '{variable}' not found")
    var = nc.variables[variable]
    dims = list(var.dimensions)
    required = {"plev", "lat", "lon"}
    if not required.issubset(dims):
        raise ValueError(
            f"variable '{variable}' must contain plev/lat/lon dimensions; got {tuple(dims)}"
        )

    plev = np.asarray(nc.variables["plev"][:], dtype=float)
    pmatch = np.flatnonzero(np.isclose(plev, float(pressure), rtol=0.0, atol=1.0e-6))
    if pmatch.size != 1:
        available = ", ".join(f"{value:g}" for value in plev)
        raise ValueError(
            f"pressure level {pressure:g} hPa is not present; available levels: {available}"
        )

    indexer = []
    selected = []
    for dim in dims:
        if dim == "plev":
            indexer.append(int(pmatch[0]))
        elif dim in CLASS_OPTIONS:
            code = selections[dim]
            if code is None:
                option = CLASS_OPTIONS[dim].replace("_", "-")
                raise ValueError(f"variable '{variable}' requires --{option}")
            cindex = coordinate_index(nc, dim, code)
            indexer.append(cindex)
            selected.append((dim, int(code), class_meaning(nc, dim, cindex, code)))
        elif dim in ("lat", "lon"):
            indexer.append(slice(None))
        else:
            raise ValueError(f"unsupported dimension '{dim}' in variable '{variable}'")

    for dim, code in selections.items():
        if code is not None and dim not in dims:
            option = CLASS_OPTIONS[dim].replace("_", "-")
            raise ValueError(f"--{option} does not apply to variable '{variable}'")

    data = np.ma.asarray(var[tuple(indexer)], dtype=float)
    remaining = [dim for dim, item in zip(dims, indexer) if isinstance(item, slice)]
    if remaining == ["lon", "lat"]:
        data = data.T
    elif remaining != ["lat", "lon"]:
        raise ValueError(f"unexpected remaining dimension order: {remaining}")
    fill = getattr(var, "_FillValue", None)
    if fill is not None:
        data = np.ma.masked_where(data == float(fill), data)
    data = np.ma.masked_invalid(data)
    return data, var, selected, float(plev[pmatch[0]])


def plot_limits(variable, units, data, user_vmin, user_vmax):
    valid = data.compressed()
    if not valid.size:
        raise ValueError(f"variable '{variable}' has no valid values at the requested selection")
    if user_vmin is not None:
        vmin = user_vmin
    elif units == "%":
        vmin = 0.0
    else:
        vmin = float(valid.min())
    if user_vmax is not None:
        vmax = user_vmax
    elif units == "%":
        vmax = 100.0
    else:
        vmax = float(valid.max())
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin > vmax:
        raise ValueError(f"invalid color range: vmin={vmin}, vmax={vmax}")
    if vmin == vmax:
        vmax = vmin + 1.0
    return vmin, vmax


def lat_formatter(value, position=None):
    rounded = int(round(value))
    if rounded == 0:
        return "Eq"
    return f"{abs(rounded)}N" if rounded > 0 else f"{abs(rounded)}S"


def lon_formatter(value, position=None):
    rounded = int(round(((value + 180) % 360) - 180))
    if rounded == 0:
        return "0"
    if abs(rounded) == 180:
        return "180"
    if rounded == 90:
        return "90E"
    if rounded == -90:
        return "90W"
    return ""


def output_stem(args, selected, pressure):
    parts = ["c4_3d", period_token(args.nc), args.variable]
    for dim, code, meaning in selected:
        label = {"particle_type": "type", "particle_category": "category",
                 "type_qc_class": "typeqc"}[dim]
        parts.append(f"{label}{code}")
    parts.append(f"{pressure:g}hPa")
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
        lon = np.asarray(nc.variables["lon"][:], dtype=float)
        data, var, selected, pressure = read_level_map(
            nc, args.variable, args.plev, selections
        )
        long_name = getattr(var, "long_name", args.variable)
        units = getattr(var, "units", "")

    vmin, vmax = plot_limits(args.variable, units, data, args.vmin, args.vmax)
    cmap = plt.get_cmap(args.cmap).copy()
    # Match the existing map2d behavior: missing cells expose the map background,
    # while valid zero values retain the lowest colormap color.
    cmap.set_bad((0.0, 0.0, 0.0, 0.0))
    lon2d, lat2d = np.meshgrid(lon, lat)

    fig = plt.figure(figsize=(9.2, 4.9))
    projection = ccrs.Robinson(central_longitude=180.0)
    data_crs = ccrs.PlateCarree()
    ax = plt.axes(projection=projection)
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.6)
    grid = ax.gridlines(crs=data_crs, draw_labels=True, linewidth=0.4,
                        color="gray", alpha=0.55, linestyle="--")
    grid.top_labels = False
    grid.bottom_labels = True
    grid.left_labels = True
    grid.right_labels = True
    grid.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    grid.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60])
    grid.xformatter = mticker.FuncFormatter(lon_formatter)
    grid.yformatter = mticker.FuncFormatter(lat_formatter)
    grid.xlabel_style = {"size": 10}
    grid.ylabel_style = {"size": 10}
    grid.xpadding = 3
    grid.ypadding = 3

    mesh = ax.pcolormesh(lon2d, lat2d, data, transform=data_crs,
                         shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    class_text = ""
    if selected:
        _, code, meaning = selected[0]
        class_text = f", code {code}: {meaning}"
    ax.set_title(f"{long_name} ({pressure:g} hPa{class_text})", fontsize=13)
    colorbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.06,
                            fraction=0.035, aspect=45)
    colorbar.ax.tick_params(labelsize=11)
    if units:
        colorbar.set_label(units, fontsize=11)

    os.makedirs(args.outdir, exist_ok=True)
    stem = Path(args.outdir) / output_stem(args, selected, pressure)
    outputs = []
    if args.write_pdf:
        output = Path(f"{stem}.pdf")
        fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
        outputs.append(output)
    if args.write_png:
        output = Path(f"{stem}.png")
        fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
        outputs.append(output)
    plt.close(fig)
    print(f"[OK] {args.variable} -> {', '.join(str(path) for path in outputs)}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot one exact pressure level from an EarthCARE C4 3D NetCDF file."
    )
    parser.add_argument("nc", help="Input NetCDF file (ecare_c4_3d_<period>.nc)")
    parser.add_argument("-o", "--outdir", default="figs",
                        help="Output directory (default: figs)")
    parser.add_argument("--variable", default="c4_cloud_fraction",
                        help="C4 variable to plot (default: c4_cloud_fraction)")
    parser.add_argument("--plev", type=float, default=500.0,
                        help="Exact pressure coordinate in hPa (default: 500)")
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
