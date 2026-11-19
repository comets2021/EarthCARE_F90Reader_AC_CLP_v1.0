#!/usr/bin/env python3
#
# USAGE: python3.11 plot_map2d.py [filename.nc] -o figs --pdf

import os
import argparse
import re
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as mticker

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except Exception as e:
    raise ImportError(
        "This script requires 'cartopy'. Install e.g.:\n"
        "  pip install cartopy\n"
        "or (conda):\n"
        "  conda install -c conda-forge cartopy"
    ) from e


def _get_fill_value(da):
    fv = da.attrs.get("_FillValue", None)
    if fv is None:
        fv = da.encoding.get("_FillValue", None)
    return fv


def format_units(u: str) -> str:
    """Convert common unit strings to mathtext superscripts."""
    if not u:
        return ""
    s = u.strip()
    s = s.replace("g m-2", r"g m$^{-2}$").replace("g m -2", r"g m$^{-2}$")
    s = s.replace("m s-1", r"m s$^{-1}$").replace("m s -1", r"m s$^{-1}$")
    s = s.replace("kg kg-1", r"kg kg$^{-1}$").replace("kg kg -1", r"kg kg$^{-1}$")
    s = s.replace("# cm-3", r"# cm$^{-3}$").replace("# cm -3", r"# cm$^{-3}$")
    s = s.replace("g m-3", r"g m$^{-3}$").replace("g m -3", r"g m$^{-3}$")
    if s.lower() == "none":
        return ""
    return s


def nice_ticks(vmin, vmax, n=6):
    if vmin is None or vmax is None:
        return None
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        return None
    return np.linspace(vmin, vmax, n)


def area_weighted_global_mean(z2d, lat1d):
    """Area-weighted mean using cos(lat). NaNs ignored."""
    w = np.cos(np.deg2rad(lat1d)).astype(float)
    w2 = w[:, None]
    mask = np.isfinite(z2d)
    if not np.any(mask):
        return np.nan
    z_fill = np.where(mask, z2d, 0.0)
    num = np.sum(z_fill * w2)
    den = np.sum(w2 * mask)
    return num / den if den != 0 else np.nan


def get_map_spec(vname: str, z: np.ndarray):
    """Choose colormap / limits / extend / ticks."""
    spec = {"vmin": None, "vmax": None, "cmap": "turbo", "extend": "neither", "ticks": None}

    finite = np.isfinite(z)
    zmin = float(np.nanmin(z)) if np.any(finite) else np.nan
    zmax = float(np.nanmax(z)) if np.any(finite) else np.nan

    if vname in ("pdf_clr", "pdf_cld", "pdf_drz", "pdf_prc"):
        spec.update(vmin=0.0, vmax=100.0, cmap="turbo", extend="neither")
        spec["ticks"] = nice_ticks(0.0, 100.0, n=6)

    elif vname == "cot":
        spec.update(vmin=0.0, vmax=20.0, cmap="turbo", extend="max")
        spec["ticks"] = nice_ticks(0.0, 20.0, n=6)

    elif vname in ("lwp", "iwp"):
        spec.update(vmin=0.0, vmax=150.0, cmap="turbo", extend="max")
        spec["ticks"] = nice_ticks(0.0, 150.0, n=6)

    elif vname == "cth":
        spec.update(vmin=0.0, vmax=15000.0, cmap="turbo", extend="max")
        spec["ticks"] = nice_ticks(0.0, 15000.0, n=6)

    elif vname == "n_obs2d":
        if np.isfinite(zmin) and np.isfinite(zmax):
            spec.update(vmin=zmin, vmax=zmax, cmap="turbo", extend="neither")
            spec["ticks"] = nice_ticks(zmin, zmax, n=6)

    else:
        if np.isfinite(zmin) and np.isfinite(zmax) and (zmin < 0.0) and (zmax > 0.0):
            vmax = max(abs(zmin), abs(zmax))
            spec.update(vmin=-vmax, vmax=vmax, cmap="RdBu_r", extend="neither")
            spec["ticks"] = nice_ticks(-vmax, vmax, n=7)
        else:
            if np.isfinite(zmin) and np.isfinite(zmax):
                spec.update(vmin=zmin, vmax=zmax, cmap="turbo", extend="neither")
                spec["ticks"] = nice_ticks(zmin, zmax, n=7)

    return spec


def lat_label_eq_ns(y, pos=None):
    yi = int(round(y))
    if yi == 0:
        return "Eq"
    return f"{abs(yi)}N" if yi > 0 else f"{abs(yi)}S"


def _get_gl_xlabel_artists(gl):
    """Cartopy version-safe getter for gridliner x-label Text artists."""
    for attr in ("xlabel_artists", "_xlabel_artists"):
        if hasattr(gl, attr):
            arts = getattr(gl, attr)
            if arts is not None:
                return list(arts)
    return []


def add_right_edge_zero_aligned(fig, ax, gl, data_crs, *,
                                eps=1e-3,
                                lat_for_x=-88.0,
                                ref_text="180",
                                y_shift_px=0,
                                fontsize=10):
    """
    Place extra '0' at the RIGHT edge, with y aligned to Cartopy bottom lon labels.

    We align using the *bottom of the rendered bbox* (y0) of a reference bottom label
    (default: '180'). Then apply a small pixel shift (negative -> down).
    """
    # Ensure gridliner has created label artists
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # --- reference label bbox y0 (display pixels) ---
    xlabels = _get_gl_xlabel_artists(gl)

    ref_bbox_y0 = None
    # 1) try exact match
    for t in xlabels:
        if (t is None) or (not t.get_visible()):
            continue
        if t.get_text().strip() == ref_text:
            bb = t.get_window_extent(renderer=renderer)
            ref_bbox_y0 = bb.y0
            break

    # 2) fallback: median y0 among visible bottom labels
    if ref_bbox_y0 is None:
        y0s = []
        for t in xlabels:
            if (t is None) or (not t.get_visible()):
                continue
            s = t.get_text().strip()
            if not s:
                continue
            bb = t.get_window_extent(renderer=renderer)
            y0s.append(bb.y0)
        if len(y0s) > 0:
            ref_bbox_y0 = float(np.median(y0s))
        else:
            ax_bb = ax.get_window_extent(renderer=renderer)
            ref_bbox_y0 = ax_bb.y0 - 10.0

    y_ref = float(ref_bbox_y0) + float(y_shift_px)

    # --- choose seam side that projects to the RIGHT edge ---
    candidates = (-eps, eps)
    xdisp = []
    for lon in candidates:
        xp, yp = ax.projection.transform_point(lon, lat_for_x, data_crs)
        xd, yd = ax.transData.transform((xp, yp))
        xdisp.append(xd)
    lon_right = candidates[int(np.argmax(xdisp))]

    # --- compute x display coord near right boundary ---
    xp, yp = ax.projection.transform_point(lon_right, lat_for_x, data_crs)
    x_ref, _ = ax.transData.transform((xp, yp))

    # --- place text in figure fraction coords ---
    xf, yf = fig.transFigure.inverted().transform((x_ref, y_ref))
    t0 = fig.text(
        xf, yf, "0",
        transform=fig.transFigure,
        ha="center",
        va="bottom",
        fontsize=fontsize
    )
    t0.set_in_layout(False)  # prevent bbox_inches='tight' from expanding

    return t0, lon_right


def add_mean_at_top_right(fig, ax, data_crs, lon_right, gmean, *,
                          lat_top=89.0,
                          x_offset_px=20,
                          fontsize=10):
    """
    Place mean text near the map's top edge (near 90N) on the right side,
    slightly to the right of the right-edge seam (lon_right).
    """
    if not np.isfinite(gmean):
        return None

    xp, yp = ax.projection.transform_point(lon_right, lat_top, data_crs)
    x_disp, y_disp = ax.transData.transform((xp, yp))

    x_disp2 = x_disp + float(x_offset_px)
    xf, yf = fig.transFigure.inverted().transform((x_disp2, y_disp))
    t = fig.text(
        xf, yf, f"mean = {gmean:.1f}",
        transform=fig.transFigure,
        ha="left",
        va="center",
        fontsize=fontsize
    )
    t.set_in_layout(False)
    return t


def period_token(nc_path):
    stem = os.path.splitext(os.path.basename(nc_path))[0]
    prefix = "ecare_map2d_"
    if not stem.startswith(prefix):
        raise ValueError(f"cannot identify period token from input filename: {stem}")
    token = stem[len(prefix):]
    if not re.fullmatch(r"(?:[0-9]{6}|[0-9]{4}|all)", token):
        raise ValueError(f"unsupported period token '{token}' in input filename")
    return token


def plot_maps(nc_path: str, out_dir: str, write_png: bool = False, dpi: int = 200):
    os.makedirs(out_dir, exist_ok=True)
    token = period_token(nc_path)

    ds = xr.open_dataset(nc_path, decode_times=False)
    for k in ("lat", "lon"):
        if k not in ds:
            raise ValueError(f"Dataset must contain coordinate '{k}'.")

    lat = ds["lat"].values
    lon = ds["lon"].values
    Lon, Lat = np.meshgrid(lon, lat)

    # 0 at edges, 180 at center
    proj = ccrs.Robinson(central_longitude=180.0)
    data_crs = ccrs.PlateCarree()

    variables = [(name, da) for name, da in ds.data_vars.items()
                 if set(da.dims) == {"lat", "lon"}]
    pdf_path = os.path.join(out_dir, f"all_map2d_{token}.pdf")
    pdf = PdfPages(pdf_path) if len(variables) > 1 else None

    # Gridlines every 30°
    grid_xlocs = np.arange(-180, 181, 30)
    grid_ylocs = [-60, -30, 0, 30, 60]

    # Only label these longitudes (90° spacing)
    label_lons = {-180, -90, 0, 90, 180}

    def lon_formatter_selective(x, pos=None):
        """
        x comes from grid_xlocs (PlateCarree). Normalize to [-180,180].
        Label only 0, ±90, ±180 as 0 / 90E / 180 / 90W.
        """
        xr_ = int(round(((x + 180) % 360) - 180))  # [-180,180]
        if xr_ not in label_lons:
            return ""
        if xr_ == 0:
            return "0"
        if abs(xr_) == 180:
            return "180"
        if xr_ == 90:
            return "90E"
        if xr_ == -90:
            return "90W"
        return ""

    for vname, da in variables:

        fv = _get_fill_value(da)
        da2 = da.where(da != fv) if fv is not None else da
        z = da2.values

        spec = get_map_spec(vname, z)
        gmean = area_weighted_global_mean(z, lat)

        fig = plt.figure(figsize=(9.2, 4.9))
        ax = plt.axes(projection=proj)
        ax.set_global()

        # Coastlines only (no borders)
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.6)

        gl = ax.gridlines(
            crs=data_crs,
            draw_labels=True,
            linewidth=0.4,
            color="gray",
            alpha=0.55,
            linestyle="--",
        )
        gl.top_labels = False
        gl.bottom_labels = True
        gl.left_labels = True
        gl.right_labels = True

        gl.xlocator = mticker.FixedLocator(grid_xlocs)
        gl.ylocator = mticker.FixedLocator(grid_ylocs)
        gl.xformatter = mticker.FuncFormatter(lon_formatter_selective)
        gl.yformatter = mticker.FuncFormatter(lat_label_eq_ns)

        gl.xlabel_style = {"size": 10}
        gl.ylabel_style = {"size": 10}
        gl.xpadding = 3
        gl.ypadding = 3

        pcm = ax.pcolormesh(
            Lon, Lat, z,
            transform=data_crs,
            shading="auto",
            cmap=spec["cmap"],
            vmin=spec["vmin"],
            vmax=spec["vmax"],
        )

        long_name = da.attrs.get("long_name", vname)
        units = format_units(da.attrs.get("units", ""))

        ax.set_title(long_name, fontsize=13)

        # Slim colorbar
        cb = fig.colorbar(
            pcm, ax=ax,
            orientation="horizontal",
            pad=0.06,
            fraction=0.035,
            aspect=45,
            extend=spec["extend"]
        )
        if spec["ticks"] is not None:
            cb.set_ticks(spec["ticks"])
            cb.update_ticks()
        cb.ax.tick_params(labelsize=11)

        # Unit label
        if units:
            xpad = 1.06 if spec["extend"] != "neither" else 1.02
            tunit = cb.ax.text(
                xpad, 0.5, units,
                transform=cb.ax.transAxes,
                va="center",
                ha="left",
                fontsize=11
            )
            tunit.set_in_layout(False)

        # ---- Extra "0" on the right edge (aligned with other lon labels) ----
        _, lon_right = add_right_edge_zero_aligned(
            fig, ax, gl, data_crs,
            eps=1e-3,
            lat_for_x=-88.0,
            ref_text="180",
            y_shift_px=0,   # if still slightly high/low, tweak -1..-3
            fontsize=10
        )

        # ---- Mean text near top boundary on the right side ----
        add_mean_at_top_right(
            fig, ax, data_crs,
            lon_right=lon_right,
            gmean=gmean,
            lat_top=89.0,
            x_offset_px=20,
            fontsize=10
        )

        stem = f"map2d_{token}_{vname}"
        out_pdf = os.path.join(out_dir, f"{stem}.pdf")
        fig.savefig(out_pdf, dpi=dpi, bbox_inches="tight")
        outputs = [out_pdf]
        if write_png:
            out_png = os.path.join(out_dir, f"{stem}.png")
            fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
            outputs.append(out_png)
        if pdf is not None:
            pdf.savefig(fig, bbox_inches="tight")

        plt.close(fig)
        print(f"[OK] {vname} -> {', '.join(outputs)}")

    if pdf is not None:
        pdf.close()
        print(f"[OK] PDF saved -> {pdf_path}")

    ds.close()


def main():
    ap = argparse.ArgumentParser(
        description="Plot EarthCARE L3 2D maps (lat,lon) on Robinson projection for all 2D vars."
    )
    ap.add_argument("nc", help="Input NetCDF file (e.g., ecare_map2d_202512.nc)")
    ap.add_argument("-o", "--outdir", default="figs", help="Output directory (default: figs)")
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
    plot_maps(args.nc, args.outdir, write_png=(args.png or args.pdf), dpi=args.dpi)


if __name__ == "__main__":
    try:
        main()
    except (KeyError, OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
