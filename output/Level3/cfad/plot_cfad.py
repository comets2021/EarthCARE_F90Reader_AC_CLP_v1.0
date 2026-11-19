#!/usr/bin/env python3
"""Plot EarthCARE CFAD panels using the legacy GrADS layout and values."""

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
from netCDF4 import Dataset


LEVELS = np.array([0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20])
GRADS_COLORS = [
    "#1e3cff", "#00a0ff", "#00dc00", "#a0e632", "#e6dc32",
    "#e6af2d", "#f08228", "#fa3c3c", "#f00082",
]
REGION_LONG = {
    "trwp": "Tropical Warm Pool", "itcz": "ITCZ", "spcz": "SPCZ",
    "nepa": "N.E. Pacific", "cast": "California StCu Deck", "peru": "Peruvian",
    "natl": "North Atlantic", "nami": "Namibian", "aust": "Australian",
    "japn": "Japan", "eqct": "Eqt. Cold Tongue", "easi": "Eastern Asia",
    "soce": "Southern Ocean", "arct": "Arctic", "anta": "Antarctic",
    "glob": "Global",
}
PAGES = [
    ["trwp", "itcz", "spcz", "nepa", "cast", "peru", "natl", "nami", "aust"],
    ["japn", "eqct", "easi", "soce", "arct", "anta", "glob"],
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
})


def period_token(nc_path):
    stem = Path(nc_path).stem
    prefix = "ecare_cfad_"
    if not stem.startswith(prefix):
        raise ValueError(f"cannot identify period token from input filename: {stem}")
    token = stem[len(prefix):]
    if token != "all" and re.fullmatch(r"\d{4}(?:\d{2})?", token) is None:
        raise ValueError(f"unsupported period token in input filename: {token}")
    return token


def read_field(nc, region):
    name = f"cfad_{region}"
    if name not in nc.variables:
        raise KeyError(f"variable '{name}' not found")
    var = nc.variables[name]
    data = np.ma.asarray(var[:], dtype=float)
    fill = getattr(var, "_FillValue", None)
    if fill is not None:
        data = np.ma.masked_where(data == float(fill), data)
    data = np.ma.masked_invalid(data)
    if data.ndim != 2:
        raise ValueError(f"variable '{name}' must be 2D, got shape {data.shape}")
    return data


def draw_page(nc, dbz, lev, regions, page_number, output_paths, dpi, cmap_name,
              combined_pdf=None):
    fig, axes = plt.subplots(3, 3, figsize=(11.0, 8.5), squeeze=False)
    fig.subplots_adjust(left=0.065, right=0.91, bottom=0.08, top=0.95,
                        wspace=0.22, hspace=0.42)

    if cmap_name == "grads":
        cmap = ListedColormap(GRADS_COLORS, name="legacy_grads_cfad")
    else:
        cmap = plt.get_cmap(cmap_name, len(GRADS_COLORS))
    cmap = cmap.copy()
    cmap.set_bad("0.85")
    norm = BoundaryNorm(LEVELS, cmap.N, extend="both")
    last_by_row = {}
    mesh_by_row = {}

    for index, region in enumerate(regions):
        row, col = divmod(index, 3)
        ax = axes[row, col]
        field = read_field(nc, region)
        if field.shape == (len(dbz), len(lev)):
            field = field.T
        if field.shape != (len(lev), len(dbz)):
            raise ValueError(
                f"cfad_{region} shape {field.shape} does not match lev/dbz "
                f"({len(lev)}, {len(dbz)})"
            )

        # GrADS uses "gxout shade2": discrete filled contours at the stated
        # thresholds, without any scientific re-normalization.
        mesh = ax.contourf(dbz, lev, field, levels=LEVELS,
                           cmap=cmap, norm=norm, extend="both")
        ax.set_xlim(-50, 20)
        ax.set_ylim(0, 20)
        ax.set_xticks(np.arange(-50, 21, 10))
        ax.set_yticks(np.arange(0, 21, 5))
        ax.tick_params(labelsize=11)
        ax.set_xlabel("Reflectivity [dBZe]", fontsize=12)
        if col == 0:
            ax.set_ylabel("Altitude [km]", fontsize=12)
        else:
            ax.set_ylabel("")
        ax.text(0.00, 1.01, "EarthCARE", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=10)
        ax.text(1.00, 1.01, REGION_LONG[region], transform=ax.transAxes,
                ha="right", va="bottom", fontsize=10)
        letter = chr(ord("a") + (page_number - 1) * 9 + index)
        ax.text(0.045, 0.94, f"({letter})", transform=ax.transAxes,
                ha="left", va="top", fontsize=12)
        last_by_row[row] = ax
        mesh_by_row[row] = mesh

    for index in range(len(regions), 9):
        row, col = divmod(index, 3)
        axes[row, col].axis("off")

    # The GrADS script places one vertical key to the right of the last panel
    # in every occupied row (beside the global panel on page 2, row 3).
    for row, ax in last_by_row.items():
        box = ax.get_position()
        cax = fig.add_axes([box.x1 + 0.010, box.y0, 0.012, box.height])
        cb = fig.colorbar(mesh_by_row[row], cax=cax, ticks=LEVELS,
                          spacing="uniform")
        cb.ax.tick_params(labelsize=8)

    for output in output_paths:
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
    if combined_pdf is not None:
        combined_pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] page {page_number} -> {', '.join(str(path) for path in output_paths)}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot two-page EarthCARE CFAD panels matching plot_cfad.gs."
    )
    parser.add_argument("nc", help="Input NetCDF file (e.g., ecare_cfad_all.nc)")
    parser.add_argument("-o", "--outdir", default="figs",
                        help="Output directory (default: figs)")
    parser.add_argument("--format", choices=("png", "pdf"), default=None,
                        help="Deprecated compatibility format: png-only or pdf-only")
    parser.add_argument("--png", action="store_true",
                        help="Also save PNG files (PDF is the standard output)")
    parser.add_argument("--dpi", type=int, default=200,
                        help="Output resolution in dpi (default: 200)")
    parser.add_argument("--cmap", default="grads",
                        help="Colormap name, or 'grads' for the legacy colors")
    args = parser.parse_args()

    if not os.path.isfile(args.nc):
        parser.error(f"input NetCDF not found: {args.nc}")
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    if args.format == "png" and args.png:
        parser.error("--format png and --png cannot be used together")
    write_pdf = args.format != "png"
    write_png = args.png or args.format == "png"
    token = period_token(args.nc)
    os.makedirs(args.outdir, exist_ok=True)

    with Dataset(args.nc, "r") as nc:
        for coordinate in ("dbz", "lev"):
            if coordinate not in nc.variables:
                raise KeyError(f"coordinate '{coordinate}' not found")
        dbz = np.asarray(nc.variables["dbz"][:], dtype=float)
        lev = np.asarray(nc.variables["lev"][:], dtype=float)
        merged_path = Path(args.outdir) / f"all_cfad_{token}.pdf"
        combined_pdf = PdfPages(merged_path) if write_pdf and len(PAGES) > 1 else None
        for page_number, regions in enumerate(PAGES, start=1):
            stem = Path(args.outdir) / f"cfad_{token}_page{page_number}"
            outputs = []
            if write_pdf:
                outputs.append(Path(f"{stem}.pdf"))
            if write_png:
                outputs.append(Path(f"{stem}.png"))
            draw_page(nc, dbz, lev, regions, page_number, outputs,
                      args.dpi, args.cmap, combined_pdf)
        if combined_pdf is not None:
            combined_pdf.close()
            print(f"[OK] merged PDF -> {merged_path}")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
