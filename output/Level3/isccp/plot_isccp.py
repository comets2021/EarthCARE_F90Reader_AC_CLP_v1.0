#!/usr/bin/env python3
import argparse
import os
import re
import sys
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


# -----------------------------
# Fixed ranges / ticks
# -----------------------------
FQ_VMIN, FQ_VMAX = 0.0, 16.0
FQ_TICKS = [0, 4, 8, 12, 16]

DV_VMIN, DV_VMAX = -2.0, 2.0
DV_TICKS = [-2.0, -1.0, 0.0, 1.0, 2.0]

# Colormaps
FQ_CMAP_DEFAULT = "viridis"   # keep as-is
DV_CMAP_NAME = "RdBu_r"       # revert DV scheme


# -----------------------------
# Text sizes
# -----------------------------
FS_TITLE = 13
FS_LABEL = 13
FS_TICK  = 12
FS_CBAR_LABEL = 12
FS_CBAR_TICK = 11

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
})


# -----------------------------
# Region long names
# -----------------------------
REG_LONG = {
    "trwp": "Tropical Warm Pool",
    "itcz": "ITCZ",
    "spcz": "SPCZ",
    "nepa": "North East Pacific",
    "cast": "California Stratocumulus",
    "peru": "Peruvian",
    "natl": "North Atlantic",
    "nami": "Namibian",
    "aust": "Australian",
    "japn": "Japan",
    "eqct": "Equatorial Cold Tongue",
    "easi": "Eastern Asia",
    "soce": "Southern Ocean",
    "arct": "Arctic",
    "anta": "Antarctic",
    "glob": "Global",
}


def read_var(nc: Dataset, name: str):
    if name not in nc.variables:
        raise KeyError(f"Variable '{name}' not found. Available: {list(nc.variables.keys())}")
    v = nc.variables[name]
    a = np.array(v[:], dtype=float)

    fv = getattr(v, "_FillValue", None)
    mv = getattr(v, "missing_value", None)
    if fv is not None:
        a = np.where(a == float(fv), np.nan, a)
    if mv is not None:
        a = np.where(a == float(mv), np.nan, a)

    # common fallback for your files
    a = np.where(a <= -900.0, np.nan, a)
    return a, v


def bounds_to_edges(bnds: np.ndarray) -> np.ndarray:
    b = np.array(bnds, dtype=float)
    if b.ndim != 2:
        raise ValueError(f"bounds must be 2D, got {b.shape}")

    if b.shape[0] == 2:
        lower = b[0, :]
        upper = b[1, :]
    elif b.shape[1] == 2:
        lower = b[:, 0]
        upper = b[:, 1]
    else:
        raise ValueError(f"bounds must have a dimension of size 2, got {b.shape}")

    return np.concatenate([lower, [upper[-1]]])


def ensure_ctp_cot_order(data2d: np.ndarray, cot_n: int, ctp_n: int) -> np.ndarray:
    if data2d.ndim != 2:
        raise ValueError(f"Expected 2D array, got {data2d.ndim}D")
    ny, nx = data2d.shape
    if (ny == ctp_n) and (nx == cot_n):
        return data2d
    if (ny == cot_n) and (nx == ctp_n):
        return data2d.T
    raise ValueError(f"Shape mismatch: data {data2d.shape}, cot {cot_n}, ctp {ctp_n}")


def add_grid_lines(ax, cot_edges, ctp_edges,
                   cot_solid=(3.6, 23.0), ctp_solid=(440.0, 680.0)):
    nxc = len(cot_edges) - 1
    nyc = len(ctp_edges) - 1

    dotted_style = (0, (1.0, 3.0))  # black thin dotted
    col = "0.0"

    # dotted interior boundaries ONLY
    for k in range(1, nxc):
        ax.axvline(k, color=col, linestyle=dotted_style, linewidth=0.6, zorder=2)
    for k in range(1, nyc):
        ax.axhline(k, color=col, linestyle=dotted_style, linewidth=0.6, zorder=2)

    # solid classification lines (overwrite dotted)
    for val in cot_solid:
        idx = np.where(np.isclose(cot_edges, val))[0]
        if idx.size > 0:
            k = int(idx[0])
            ax.axvline(k, color="0.0", linestyle="-", linewidth=1.3, zorder=4)

    for val in ctp_solid:
        idx = np.where(np.isclose(ctp_edges, val))[0]
        if idx.size > 0:
            k = int(idx[0])
            ax.axhline(k, color="0.0", linestyle="-", linewidth=1.3, zorder=4)


def plot_one(data, cot_edges, ctp_edges, region_long, output_paths,
             cmap, norm=None, vmin=None, vmax=None,
             cbar_label="", cbar_ticks=None, dpi=200, combined_pdf=None):
    nxc = len(cot_edges) - 1
    nyc = len(ctp_edges) - 1

    x = np.arange(nxc + 1)
    y = np.arange(nyc + 1)

    fig, ax = plt.subplots(figsize=(6, 5.2), constrained_layout=False)

    m = ax.pcolormesh(
        x, y, data,
        shading="flat",
        cmap=cmap,
        norm=norm,
        vmin=vmin, vmax=vmax,
        zorder=1
    )

    # invert y so small pressure is at top
    ax.set_ylim(nyc, 0)
    ax.set_xlim(0, nxc)

    # axis labels
    ax.set_xlabel(r"Cloud Optical Thickness $\tau_c$", fontsize=FS_LABEL)
    ax.set_ylabel(r"Cloud-Top Pressure $p_{ctop}$ (hPa)", fontsize=FS_LABEL)

    # boundary ticks
    ax.set_xticks(np.arange(nxc + 1))
    ax.set_xticklabels([f"{v:g}" for v in cot_edges], fontsize=FS_TICK)
    ax.set_yticks(np.arange(nyc + 1))
    ax.set_yticklabels([f"{v:g}" for v in ctp_edges], fontsize=FS_TICK)

    # grid + classification lines
    add_grid_lines(ax, cot_edges, ctp_edges,
                   cot_solid=(3.6, 23.0),
                   ctp_solid=(440.0, 680.0))

    # frame/spines drawn last
    for spine in ax.spines.values():
        spine.set_linewidth(1.4)
        spine.set_zorder(6)

    # colorbar with same height as axes
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.10)
    cb = fig.colorbar(m, cax=cax, ticks=cbar_ticks)
    cb.set_label(cbar_label, fontsize=FS_CBAR_LABEL)
    cb.ax.tick_params(labelsize=FS_CBAR_TICK)

    # region label aligned to diagram right edge (not over colorbar)
    ax.text(
        1.00, 1.01, region_long,
        transform=ax.transAxes,
        ha="right", va="bottom",
        fontsize=FS_TITLE
    )

    fig.subplots_adjust(left=0.14, right=0.86, bottom=0.14, top=0.90)
    for output in output_paths:
        fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    if combined_pdf is not None:
        combined_pdf.savefig(fig, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    for path in output_paths:
        output_type = "PDF" if path.suffix.lower() == ".pdf" else "PNG"
        print(f"[OK] {output_type} saved -> {path}")


def period_token(nc_path):
    stem = Path(nc_path).stem
    prefix = "ecare_isccp_"
    if not stem.startswith(prefix):
        raise ValueError(f"cannot identify period token from input filename: {stem}")
    token = stem[len(prefix):]
    if token != "all" and re.fullmatch(r"\d{4}(?:\d{2})?", token) is None:
        raise ValueError(f"unsupported period token in input filename: {token}")
    return token


def output_stem(args, metric, region):
    if args.legacy_prefix:
        return Path(f"{args.legacy_prefix}_{metric}_{region}")
    name = f"isccp_{args.diagnostic}_{period_token(args.nc)}_{metric}_{region}"
    return Path(args.outdir) / name


def output_paths(args, metric, region, write_pdf, write_png):
    stem = output_stem(args, metric, region)
    paths = []
    if write_pdf:
        paths.append(Path(f"{stem}.pdf"))
    if write_png:
        paths.append(Path(f"{stem}.png"))
    return paths


def main():
    ap = argparse.ArgumentParser(
        description="Plot legacy or C4 EarthCARE ISCCP-like histograms without re-normalization."
    )
    ap.add_argument("nc", help="Input NetCDF file (e.g., ecare_isccp_all.nc)")
    ap.add_argument("-o", "--outdir", default="figs",
                    help="Output directory (default: figs)")
    ap.add_argument(
        "--region",
        default="all",
        help="Region code (e.g., arct) or 'all' to process all regions"
    )
    ap.add_argument("--diagnostic", choices=("legacy", "c4"), default="legacy",
                    help="ISCCP population to plot (default: legacy)")
    ap.add_argument("--fq-cmap", default=FQ_CMAP_DEFAULT,
                    help="Frequency colormap (e.g., viridis, cividis, YlGnBu, magma)")
    ap.add_argument("--format", choices=("png", "pdf"), default=None,
                    help="Deprecated compatibility format: png-only or pdf-only")
    ap.add_argument("--png", action="store_true",
                    help="Also save PNG files (PDF is the standard output)")
    ap.add_argument("--dpi", type=int, default=200,
                    help="Output resolution in dpi (default: 200)")
    ap.add_argument("--out", dest="legacy_prefix", default=None,
                    help=argparse.SUPPRESS)
    args = ap.parse_args()

    if not os.path.isfile(args.nc):
        ap.error(f"input NetCDF not found: {args.nc}")
    if args.dpi <= 0:
        ap.error("--dpi must be positive")
    if args.format == "png" and args.png:
        ap.error("--format png and --png cannot be used together")
    token = period_token(args.nc)
    if args.legacy_prefix and args.format is None and not args.png:
        write_pdf = False
        write_png = True
    else:
        write_pdf = args.format != "png"
        write_png = args.png or args.format == "png"
    if args.legacy_prefix is None:
        os.makedirs(args.outdir, exist_ok=True)

    # decide target regions
    r = args.region.strip().lower()
    if r in ("all", "*"):
        regions = list(REG_LONG.keys())  # preserves dict insertion order (Py3.7+)
    else:
        regions = [args.region.strip()]

    written = 0
    with Dataset(args.nc, "r") as nc:
        # read common axes only once
        cot, _ = read_var(nc, "cot")
        ctp, _ = read_var(nc, "ctp")
        cot_bnds, _ = read_var(nc, "cot_bnds")
        ctp_bnds, _ = read_var(nc, "ctp_bnds")

        cot_edges = bounds_to_edges(cot_bnds)
        ctp_edges = bounds_to_edges(ctp_bnds)

        available_regions = []
        for region in regions:
            suffix = f"c4_{region}" if args.diagnostic == "c4" else region
            if (f"isccp_fq_{suffix}" in nc.variables and
                    f"isccp_dv_{suffix}" in nc.variables):
                available_regions.append(region)
            else:
                print(f"[SKIP] region='{region}': required variables are missing",
                      file=sys.stderr)

        if args.legacy_prefix:
            merged_path = Path(f"{args.legacy_prefix}_all.pdf")
        else:
            merged_path = Path(args.outdir) / f"all_isccp_{args.diagnostic}_{token}.pdf"
        combined_pdf = PdfPages(merged_path) if write_pdf and len(available_regions) * 2 > 1 else None

        for region in available_regions:
            region_long = REG_LONG.get(region, region)

            suffix = f"c4_{region}" if args.diagnostic == "c4" else region
            vfq = f"isccp_fq_{suffix}"
            vdv = f"isccp_dv_{suffix}"

            fq, _ = read_var(nc, vfq)
            dv, _ = read_var(nc, vdv)

            fq = ensure_ctp_cot_order(fq, cot_n=len(cot), ctp_n=len(ctp))
            dv = ensure_ctp_cot_order(dv, cot_n=len(cot), ctp_n=len(ctp))

            # Colormaps with explicit "bad" color so NaN is distinguishable from 0
            fq_cmap = plt.get_cmap(args.fq_cmap).copy()
            fq_cmap.set_bad("0.85")  # light gray for missing

            dv_cmap = plt.get_cmap(DV_CMAP_NAME).copy()  # reverted DV scheme
            dv_cmap.set_bad("0.85")  # light gray for missing

            dv_norm = TwoSlopeNorm(vmin=DV_VMIN, vcenter=0.0, vmax=DV_VMAX)

            plot_one(
                data=fq,
                cot_edges=cot_edges,
                ctp_edges=ctp_edges,
                region_long=region_long,
                output_paths=output_paths(args, "fq", region, write_pdf, write_png),
                cmap=fq_cmap,
                vmin=FQ_VMIN, vmax=FQ_VMAX,
                cbar_label="Frequency (%)",
                cbar_ticks=FQ_TICKS,
                dpi=args.dpi,
                combined_pdf=combined_pdf,
            )
            written += 1

            plot_one(
                data=dv,
                cot_edges=cot_edges,
                ctp_edges=ctp_edges,
                region_long=region_long,
                output_paths=output_paths(args, "dv", region, write_pdf, write_png),
                cmap=dv_cmap,
                norm=dv_norm,
                cbar_label=r"Doppler Velocity (m s$^{-1}$)",
                cbar_ticks=DV_TICKS,
                dpi=args.dpi,
                combined_pdf=combined_pdf,
            )
            written += 1

        if combined_pdf is not None:
            combined_pdf.close()
            print(f"[OK] PDF saved -> {merged_path}")

    if written == 0:
        raise ValueError(
            f"no ISCCP variables found for diagnostic='{args.diagnostic}' and region='{args.region}'"
        )


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)

