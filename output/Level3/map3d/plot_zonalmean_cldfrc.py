#!/usr/bin/env python3
#
# USAGE: python3.11 plot_zonalmean_cld3d.py [ecare_cld3d_YYYYMM.nc] -o figs --pdf

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


def nice_ticks(vmin, vmax, n=6):
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

    # Derived products are all percentages [%]
#    if vname == "cldfrc3d":
#        spec.update(vmin=0.0, vmax=100.0, cmap="turbo", extend="neither")
#        spec["ticks"] = nice_ticks(0.0, 100.0, n=6)
    # Derived products are all percentages [%]
    if vname == "cldfrc3d":
        spec.update(vmin=0.0, vmax=60.0, cmap="turbo", extend="max")
        spec["ticks"] = nice_ticks(0.0, 60.0, n=7)

    elif vname.startswith("typefrc3d_"):
        spec.update(vmin=0.0, vmax=100.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 100.0, n=6)

    else:
        # Fallback
        if np.isfinite(zmin) and np.isfinite(zmax) and (zmin < 0.0) and (zmax > 0.0):
            vmax = max(abs(zmin), abs(zmax))
            spec.update(vmin=-vmax, vmax=vmax, cmap="RdBu_r", extend="neither")
            spec["ticks"] = nice_ticks(spec["vmin"], spec["vmax"], n=7)
        elif np.isfinite(zmin) and np.isfinite(zmax):
            spec.update(vmin=zmin, vmax=zmax, cmap="turbo", extend="neither")
            spec["ticks"] = nice_ticks(spec["vmin"], spec["vmax"], n=7)

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


def _sanitize_count_da(da):
    """Replace FillValue with NaN and cast to float."""
    fv = _get_fill_value(da)
    da2 = da.astype("float64")
    if fv is not None:
        da2 = da2.where(da2 != fv)
    return da2


def build_derived_dataset(ds):
    """
    Build a dataset containing:
      - cldfrc3d = n_cld3d_0 / n_obs3d * 100
      - typefrc3d_i = n_cld3d_i / n_cld3d_0 * 100  (i=1..14)
    """
    required_coords = ("lat", "lon", "plev")
    for k in required_coords:
        if k not in ds:
            raise ValueError(f"Dataset must contain coordinate '{k}'.")

    if "n_obs3d" not in ds:
        raise ValueError("Dataset must contain 'n_obs3d'.")
    if "n_cld3d_0" not in ds:
        raise ValueError("Dataset must contain 'n_cld3d_0'.")

    n_obs = _sanitize_count_da(ds["n_obs3d"])
    n_cld0 = _sanitize_count_da(ds["n_cld3d_0"])

    # mask invalid denominators
    n_obs_valid = n_obs.where(n_obs > 0)
    n_cld0_valid = n_cld0.where(n_cld0 > 0)

    out = xr.Dataset(coords={"lat": ds["lat"], "lon": ds["lon"], "plev": ds["plev"]})

    # Cloud Fraction [%]
    cldfrc = (n_cld0_valid / n_obs_valid) * 100.0
    cldfrc.attrs["long_name"] = "Cloud Fraction"
    cldfrc.attrs["units"] = "%"
    out["cldfrc3d"] = cldfrc

    # Type fractions within cloudy scenes [%]
    type_name_map = {
        1: "warm water",
        2: "supercooled water",
        3: "3D ice",
        4: "2D plate ice",
        5: "3D+2D ice mixture",
        6: "liquid drizzle",
        7: "mixed-phase drizzle",
        8: "rain",
        9: "snow",
        10: "water and drizzle",
        11: "water and rain",
        12: "mixed phase",
        13: "unknown",
        14: "melting layer",
    }

    for i in range(1, 15):
        vsrc = f"n_cld3d_{i}"
        if vsrc not in ds:
            continue
        num = _sanitize_count_da(ds[vsrc])
        frac = (num / n_cld0_valid) * 100.0
        vdst = f"typefrc3d_{i:02d}"
        frac.attrs["long_name"] = f"Cloud Type Fraction ({type_name_map.get(i, f'type {i}')})"
        frac.attrs["units"] = "%"
        out[vdst] = frac

    return out


def _get_temp_zonalmean(ds_raw):
    """
    Return zonal-mean temperature (plev, lat) for isotherm overlay if available.
    """
    if "temp" not in ds_raw:
        return None

    daT = ds_raw["temp"]
    if not (daT.ndim == 3 and all(d in daT.dims for d in ("plev", "lat", "lon"))):
        return None

    fvT = _get_fill_value(daT)
    daT2 = daT.astype("float64")
    if fvT is not None:
        daT2 = daT2.where(daT2 != fvT)

    temp_zm = daT2.mean(dim="lon", skipna=True).transpose("plev", "lat").values
    return temp_zm


def period_token(nc_path):
    stem = os.path.splitext(os.path.basename(nc_path))[0]
    prefix = "ecare_cld3d_"
    if not stem.startswith(prefix):
        raise ValueError(f"cannot identify period token from input filename: {stem}")
    token = stem[len(prefix):]
    if token != "all" and re.fullmatch(r"\d{4}(?:\d{2})?", token) is None:
        raise ValueError(f"unsupported period token in input filename: {token}")
    return token


def plot_zonal_means(nc_path: str, out_dir: str, write_png: bool = False,
                     dpi: int = 200):
    os.makedirs(out_dir, exist_ok=True)
    token = period_token(nc_path)

    ds_raw = xr.open_dataset(nc_path, decode_times=False)

    # Ensure plev is descending (1000 -> 0) for consistent plotting
    if "plev" not in ds_raw:
        raise ValueError("Dataset must contain coordinate 'plev'.")
    plev_raw = ds_raw["plev"].values
    if np.all(np.diff(plev_raw) > 0):  # ascending
        ds_raw = ds_raw.sortby("plev", ascending=False)

    # Derived fraction dataset
    ds = build_derived_dataset(ds_raw)

    # Temperature overlay (optional)
    temp_zm = _get_temp_zonalmean(ds_raw)

    lat = ds["lat"].values
    plev = ds["plev"].values  # hPa

    lat_ticks = [-90, -60, -30, 0, 30, 60, 90]
    lat_ticklabels = ["90S", "60S", "30S", "Eq", "30N", "60N", "90N"]

    variables = [(name, da) for name, da in ds.data_vars.items()
                 if da.ndim == 3 and all(d in da.dims for d in ("plev", "lat", "lon"))]
    pdf_path = os.path.join(out_dir, f"all_cld3d_{token}.pdf")
    pdf = PdfPages(pdf_path) if len(variables) > 1 else None

    for vname, da in variables:

        # Zonal mean, force dimension order
        zm = da.mean(dim="lon", skipna=True).transpose("plev", "lat")
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

        # ---- Remove colorbar label (as requested) ----
        if spec["ticks"] is not None:
            cb.set_ticks(spec["ticks"])
            cb.update_ticks()
        cb.ax.tick_params(labelsize=11)

        # ---- Title includes units ----
        long_name = da.attrs.get("long_name", vname)
        units_raw = da.attrs.get("units", "")
        units = format_units(units_raw)
        title_txt = f"{long_name} [{units}]" if units else f"{long_name}"

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

        ax.set_title(title_txt, fontsize=13)

        # ---- Overlay isotherms (if temp exists) ----
        if temp_zm is not None:
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

        stem = f"cld3d_{token}_{vname}_zonalmean"
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

    ds_raw.close()
    ds.close()


def main():
    ap = argparse.ArgumentParser(
        description="Plot zonal-mean (lon-mean) latitude-pressure sections for EarthCARE cloud-type fractions (ecare_cld3d)."
    )
    ap.add_argument("nc", help="Input NetCDF file (e.g., ecare_cld3d_202512.nc)")
    ap.add_argument("-o", "--outdir", default="zonal_mean_figs", help="Output directory")
    ap.add_argument("--png", action="store_true",
                    help="Also save PNG files (PDF is the standard output)")
    ap.add_argument("--pdf", action="store_true",
                    help="Deprecated compatibility option: save PDF and PNG")
    ap.add_argument("--dpi", type=int, default=200,
                    help="Output resolution in dpi (default: 200)")
    args = ap.parse_args()

    if not os.path.isfile(args.nc):
        ap.error(f"input NetCDF not found: {args.nc}")
    if args.dpi <= 0:
        ap.error("--dpi must be positive")
    plot_zonal_means(args.nc, args.outdir,
                     write_png=(args.png or args.pdf), dpi=args.dpi)


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
