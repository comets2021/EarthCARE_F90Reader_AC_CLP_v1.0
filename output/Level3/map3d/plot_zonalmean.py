#!/usr/bin/env python3
#
# USAGE: python3.11 plot_zonalmean.py [filename.nc] -o figs --pdf

import os
import argparse
import re
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patheffects as pe


def _get_fill_value(da):
    """Safely get FillValue from attrs/encoding."""
    fv = da.attrs.get("_FillValue", None)
    if fv is None:
        fv = da.encoding.get("_FillValue", None)
    return fv


def format_units(u: str) -> str:
    """
    Convert common unit strings to mathtext superscripts.
    """
    if not u:
        return ""

    s = u.strip()

    # normalize spacing variants
    s = s.replace("kg kg -1", "kg kg-1")
    s = s.replace("m s -1", "m s-1")
    s = s.replace("# cm -3", "# cm-3")
    s = s.replace("g m -3", "g m-3")

    if s.lower() == "micron":
        return r"$\mu$m"

    s = s.replace("m s-1", r"m s$^{-1}$")
    s = s.replace("kg kg-1", r"kg kg$^{-1}$")
    s = s.replace("# cm-3", r"# cm$^{-3}$")
    s = s.replace("g m-3", r"g m$^{-3}$")

    return s


def nice_ticks(vmin, vmax, n=7):
    """Create ticks including endpoints."""
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        return None
    return np.linspace(vmin, vmax, n)


def get_plot_spec(vname: str, z: np.ndarray):
    """
    Decide color limits, colormap, and colorbar extension.
    Returns dict: vmin, vmax, cmap, extend, ticks
    """
    spec = {"vmin": None, "vmax": None, "cmap": "turbo", "extend": "neither", "ticks": None}

    finite = np.isfinite(z)
    if np.any(finite):
        zmin = float(np.nanmin(z))
        zmax = float(np.nanmax(z))
    else:
        zmin, zmax = np.nan, np.nan

    if vname == "temp":
        spec.update(vmin=200.0, vmax=300.0, cmap="turbo", extend="both")
        spec["ticks"] = nice_ticks(200.0, 300.0, n=6)

    elif vname == "shum":
        spec.update(vmin=0.0, vmax=0.015, cmap="turbo", extend="max")
        spec["ticks"] = nice_ticks(0.0, 0.015, n=6)

    elif vname == "rhum":
        spec.update(vmin=0.0, vmax=100.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 100.0, n=6)

    elif vname in ("cwc1", "cwc2", "ciwc"):
        spec.update(vmin=0.0, vmax=0.1, cmap="turbo", extend="max")
        spec["ticks"] = nice_ticks(0.0, 0.1, n=6)

    elif vname == "cer1":
        spec.update(vmin=0.0, vmax=100.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 100.0, n=6)

    elif vname == "cer2":
        spec.update(vmin=0.0, vmax=150.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 150.0, n=6)

    elif vname == "cerliq":
        spec.update(vmin=0.0, vmax=150.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 150.0, n=6)

    elif vname == "cerice":
        spec.update(vmin=0.0, vmax=300.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 300.0, n=6)

    elif vname in ("cldavel", "clddvel"):
        spec.update(vmin=-3.0, vmax=3.0, cmap="RdBu_r", extend="both")
        spec["ticks"] = nice_ticks(-3.0, 3.0, n=7)

    elif vname in ("cldtvel1", "cldtvel2"):
        spec.update(vmin=-0.3, vmax=0.0, cmap="Blues_r", extend="min")
        spec["ticks"] = nice_ticks(-0.3, 0.0, n=6)

    elif vname == "cldtvelt":
        spec.update(vmin=-0.5, vmax=0.0, cmap="Blues_r", extend="min")
        spec["ticks"] = nice_ticks(-0.5, 0.0, n=6)

    elif vname == "numc1":
        spec.update(vmin=0.0, vmax=100.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 100.0, n=6)

    elif vname == "numc2":
        spec.update(vmin=0.0, vmax=1.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 1.0, n=6)

    elif vname == "numliq":
        spec.update(vmin=0.0, vmax=300.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 300.0, n=6)

    elif vname == "numice":
        spec.update(vmin=0.0, vmax=5.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 5.0, n=6)

    else:
        if np.isfinite(zmin) and np.isfinite(zmax) and (zmin < 0.0) and (zmax > 0.0):
            vmax = max(abs(zmin), abs(zmax))
            spec.update(vmin=-vmax, vmax=vmax, cmap="RdBu_r", extend="neither")
            spec["ticks"] = nice_ticks(spec["vmin"], spec["vmax"], n=7)
        else:
            if np.isfinite(zmin) and np.isfinite(zmax):
                spec.update(vmin=zmin, vmax=zmax, cmap="turbo", extend="neither")
                spec["ticks"] = nice_ticks(spec["vmin"], spec["vmax"], n=7)
            else:
                spec.update(vmin=None, vmax=None, cmap="turbo", extend="neither")
                spec["ticks"] = None

    return spec


def draw_isotherm_subtle_halo(ax, lat, plev, field2d, level, color,
                             lw_color=1.1, lw_halo=2.0, alpha_halo=0.75,
                             z_halo=9, z_color=10):
    """
    Draw isotherm as white thin underlay + colored line.
    Subtle halo (not too strong).
    Returns cs_color for labeling.
    """
    ax.contour(
        lat, plev, field2d,
        levels=[level],
        colors=["white"],
        linestyles=["solid"],
        linewidths=lw_halo,
        alpha=alpha_halo,
        zorder=z_halo
    )
    cs_color = ax.contour(
        lat, plev, field2d,
        levels=[level],
        colors=[color],
        linestyles=["solid"],
        linewidths=lw_color,
        zorder=z_color
    )
    return cs_color


def period_token(nc_path):
    stem = os.path.splitext(os.path.basename(nc_path))[0]
    prefix = "ecare_map3d_"
    if not stem.startswith(prefix):
        raise ValueError(f"cannot identify period token from input filename: {stem}")
    token = stem[len(prefix):]
    if not re.fullmatch(r"(?:[0-9]{6}|[0-9]{4}|all)", token):
        raise ValueError(f"unsupported period token '{token}' in input filename")
    return token


def plot_zonal_means(nc_path: str, out_dir: str, write_png: bool = False,
                     dpi: int = 200):
    os.makedirs(out_dir, exist_ok=True)
    token = period_token(nc_path)

    ds = xr.open_dataset(nc_path, decode_times=False)

    for k in ("lat", "lon", "plev"):
        if k not in ds:
            raise ValueError(f"Dataset must contain coordinate '{k}'.")

    lat = ds["lat"].values
    plev = ds["plev"].values  # hPa

    lat_ticks = [-90, -60, -30, 0, 30, 60, 90]
    lat_ticklabels = ["90S", "60S", "30S", "Eq", "30N", "60N", "90N"]

    three_d_vars = [name for name, da in ds.data_vars.items() if set(da.dims) == {"plev", "lat", "lon"}]

    overlay_vars = set()
    if "cldavel" in three_d_vars:
        overlay_vars = set(three_d_vars[three_d_vars.index("cldavel"):])

    temp_zm = None
    if "temp" in ds and set(ds["temp"].dims) == {"plev", "lat", "lon"}:
        daT = ds["temp"]
        fvT = _get_fill_value(daT)
        daT2 = daT.where(daT != fvT) if fvT is not None else daT
        temp_zm = daT2.mean(dim="lon", skipna=True).values  # (plev, lat)

    pdf_path = os.path.join(out_dir, f"all_map3d_{token}.pdf")
    pdf = PdfPages(pdf_path) if len(three_d_vars) > 1 else None

    for vname, da in ds.data_vars.items():
        if set(da.dims) != {"plev", "lat", "lon"}:
            continue

        fv = _get_fill_value(da)
        da2 = da.where(da != fv) if fv is not None else da

        zm = da2.mean(dim="lon", skipna=True)
        z = zm.values

        spec = get_plot_spec(vname, z)

        fig = plt.figure(figsize=(8.0, 4.8))
        ax = fig.add_subplot(111)

        pcm = ax.pcolormesh(
            lat, plev, z,
            shading="auto",
            cmap=spec["cmap"],
            vmin=spec["vmin"],
            vmax=spec["vmax"]
        )

        cb = fig.colorbar(pcm, ax=ax, extend=spec["extend"])

        long_name = da.attrs.get("long_name", vname)
        units_raw = da.attrs.get("units", "")
        units = format_units(units_raw)
        cb_label = f"{long_name} [{units}]" if units else f"{long_name}"
        cb.set_label(cb_label, fontsize=12)

        if spec["ticks"] is not None:
            cb.set_ticks(spec["ticks"])
            cb.update_ticks()
        cb.ax.tick_params(labelsize=11)

        ax.set_xlim(-90, 90)
        ax.set_xticks(lat_ticks)
        ax.set_xticklabels(lat_ticklabels)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelsize=12)

        ax.set_ylabel("Pressure (hPa)", fontsize=13)
        ax.set_yscale("linear")
        ax.set_ylim(1000, 0)
        ax.set_yticks(np.arange(0, 1001, 200))
        ax.tick_params(axis="y", labelsize=12)

        ax.set_title(long_name, fontsize=13)

        # ---- Overlay isotherms: both solid, thinner, subtle halo ----
        if (vname in overlay_vars) and (temp_zm is not None):
            cs273 = draw_isotherm_subtle_halo(
                ax, lat, plev, temp_zm,
                level=273.0,
                color="#0B3D91",   # dark blue
                lw_color=1.05,
                lw_halo=1.9,
                alpha_halo=0.70
            )
            cs238 = draw_isotherm_subtle_halo(
                ax, lat, plev, temp_zm,
                level=238.0,
                color="#67B7FF",   # light blue
                lw_color=1.05,
                lw_halo=1.9,
                alpha_halo=0.70
            )

            lbl273 = ax.clabel(cs273, fmt={273.0: "273 K"}, inline=True, fontsize=9)
            lbl238 = ax.clabel(cs238, fmt={238.0: "238 K"}, inline=True, fontsize=9)
            for t in (lbl273 + lbl238):
                t.set_path_effects([pe.Stroke(linewidth=2.2, foreground="white"), pe.Normal()])

        fig.tight_layout()

        stem = f"map3d_{token}_{vname}_zonalmean"
        out_pdf = os.path.join(out_dir, f"{stem}.pdf")
        fig.savefig(out_pdf, dpi=dpi)
        outputs = [out_pdf]
        if write_png:
            out_png = os.path.join(out_dir, f"{stem}.png")
            fig.savefig(out_png, dpi=dpi)
            outputs.append(out_png)
        if pdf is not None:
            pdf.savefig(fig)

        plt.close(fig)
        print(f"[OK] {vname} -> {', '.join(outputs)}")

    if pdf is not None:
        pdf.close()
        print(f"[OK] PDF saved -> {pdf_path}")

    ds.close()


def main():
    ap = argparse.ArgumentParser(
        description="Plot zonal-mean (lon-mean) latitude-pressure sections for all 3D vars (plev,lat,lon)."
    )
    ap.add_argument("nc", help="Input NetCDF file (e.g., ecare_map3d_202512.nc)")
    ap.add_argument("-o", "--outdir", default="zonal_mean_figs", help="Output directory")
    ap.add_argument("--png", action="store_true",
                    help="Also save PNG files (PDF is the standard output)")
    ap.add_argument("--pdf", action="store_true",
                    help="Deprecated compatibility option: save PDF and PNG")
    ap.add_argument("--dpi", type=int, default=200, help="Output dpi (default: 200)")
    args = ap.parse_args()

    if not os.path.isfile(args.nc):
        ap.error(f"input NetCDF not found: {args.nc}")
    if args.dpi <= 0:
        ap.error("--dpi must be positive")
    plot_zonal_means(args.nc, args.outdir, write_png=(args.png or args.pdf), dpi=args.dpi)


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)

