"""Make the figures of the documentation from a real scan.

    python docs/make_figures.py --scan sub-01_T1w.nii.gz [--cache FILE]

- ``t1sim-scheme.png``: the path from the scan to the simulated image and
  its ground truth.
- ``t1sim-frontend.png``: what :func:`t1sim.frontend.segment` returns.
- ``t1sim-examples.png``: one detail of the anatomy under every artefact, a
  bias field on the whole slice, and an atrophy.

Every image is one axial slice of the working grid, and everything shown is
made by t1sim from the scan.  ``--cache`` keeps the slices, so that a change
of the layout does not run the simulation again (about four minutes).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from t1sim import anatomy, atlas, frontend, simulate  # noqa: E402

# One colour per tissue wherever a tissue is shown by colour, and the inks of
# the text.  The three hues stay apart for colour-blind readers as well.
TISSUE_COLOR = {"CSF": "#257de0", "GM": "#f47230", "WM": "#00b894"}
INK, INK_SOFT, SURFACE = "#191919", "#666666", "#ffffff"
DONE, DONE_EDGE = "#e3effc", "#257de0"       # a step t1sim takes
OPEN, OPEN_EDGE = "#f3f3f3", "#9c9c9c"       # a step that is MATLAB only

#: Height in mm of the slice above the origin of the template (MNI), and the
#: corner and size in mm of the enlarged detail.
SLICE_MM = 18.0
ZOOM_MM = (-52.0, 4.0, 44.0, 44.0)

#: Region of the Hammers atlas whose cortex is thinned for the examples (the
#: left middle frontal gyrus, in which the detail lies), and by how much.
ATROPHY = ("Hammers", 28, -0.25)


def _slice_from_scan(volume, mat, seg, k):
    """Slice *k* of the working grid, sampled from a volume on the grid of the scan."""
    grid = seg["grid"]
    nx, ny = grid.shape[:2]
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    voxels = np.stack([i, j, np.full_like(i, k), np.ones_like(i)]).reshape(4, -1)
    # grid voxel -> template mm -> scan mm -> voxel of the scan
    to_scan = np.linalg.inv(mat) @ np.linalg.inv(seg["world"]) @ grid.affine
    coords = (to_scan @ voxels)[:3]
    return map_coordinates(volume, coords, order=1, mode="constant").reshape(nx, ny)


def collect(scan):
    """The slices the figures show, as a dict of 2-D arrays and a few numbers."""
    seg = frontend.segment(scan, verbose=False, deformation=True)
    grid = seg["grid"]
    zooms = np.sqrt((grid.affine[:3, :3] ** 2).sum(0))
    k = int(round((SLICE_MM - grid.affine[2, 3]) / zooms[2]))

    raw = nib.as_closest_canonical(nib.load(scan))
    if not np.allclose(raw.affine, seg["native_mat"], atol=1e-4):
        raise RuntimeError("the scan is not on the grid the front end reports")
    out = {
        "scan": _slice_from_scan(raw.get_fdata(dtype=np.float32), raw.affine, seg, k),
        "corrected": seg["corrected"][:, :, k],
        "mask": seg["mask"][:, :, k],
        "label": seg["label"][:, :, k],
        "fractions": seg["fractions"][:, :, k],
        "means": np.array([seg["means"][t] for t in frontend.TISSUES]),
        "zooms": zooms,
        "origin": grid.affine[:3, 3],
    }

    truth = simulate.ground_truth(scan, seg=seg)
    out["sim_fractions"] = truth["fractions"][:, :, k]
    out["sim_ideal"] = truth["ideal"][:, :, k]
    plain = dict(snr=0, bias_percent=0, seed=1)
    coarse = dict(plain, resolution=1.0, upsample=True)
    examples = {
        "sim_image": dict(snr=30, bias_percent=30, seed=1),
        "ex_noise": dict(snr=12, bias_percent=0, seed=1),
        "ex_bias": dict(snr=0, bias_percent=-60, bias_type=(3, 2), seed=1),
        "ex_contrast_low": dict(plain, contrast=0.5),
        "ex_contrast_high": dict(plain, contrast=1.5),
        "ex_resolution": coarse,
        "ex_motion": dict(coarse, motion=2),
        "ex_ringing": dict(coarse, ringing=2),
        "ex_gibbs": dict(coarse, ringing={"strength": 2, "type": "gibbs"}),
    }
    for key, options in examples.items():
        out[key] = simulate.acquire(truth, **options)["image"][:, :, k]
    name, region, change = ATROPHY
    thin = simulate.ground_truth(scan, seg=seg, atrophy={"atlas": name, "regions": [region], "change": change})
    out["ex_atrophy"] = thin["ideal"][:, :, k]
    out["ex_atrophy_fractions"] = thin["fractions"][:, :, k]
    out["ex_atrophy_reached"] = np.array([thin["parameters"]["Atrophy"][0]["Achieved"]])
    out["roi"] = anatomy.roi_from_atlas(atlas.to_grid(name, seg["iy"], seg["template_mat"]), [region])[:, :, k]
    return out


# ------------------------------------------------------------------ drawing
def _box(d):
    """Bounding box of the head in the slice, as slices, with a small margin."""
    head = d["corrected"] > 0.08 * np.percentile(d["corrected"], 99.5)
    xs, ys = np.where(head.any(1))[0], np.where(head.any(0))[0]
    pad = 6
    return (slice(max(xs[0] - pad, 0), xs[-1] + pad), slice(max(ys[0] - pad, 0), ys[-1] + pad))


def _zoom(d):
    """The enlarged detail, as slices."""
    x0, y0, w, h = ZOOM_MM
    i0 = int(round((x0 - d["origin"][0]) / d["zooms"][0]))
    j0 = int(round((y0 - d["origin"][1]) / d["zooms"][1]))
    return (slice(i0, i0 + int(round(w / d["zooms"][0]))),
            slice(j0, j0 + int(round(h / d["zooms"][1]))))


def _show(ax, image, box, extent=None, **kwargs):
    """An axial slice with anterior up and the right side of the head on the right."""
    return ax.imshow(image[box].swapaxes(0, 1), origin="lower", extent=extent,
                     interpolation="nearest", **kwargs)


def _gray(image, d, reference="corrected"):
    """Display range of an intensity image: black to a little above white matter."""
    wm = np.median(d[reference][d["label"] > 2.9]) if reference else None
    if wm is None:
        wm = np.percentile(image, 99.0) / 1.15
    return {"cmap": "gray", "vmin": 0, "vmax": 1.25 * wm}


def _tissue_cmap(tissue):
    return LinearSegmentedColormap.from_list(tissue, [SURFACE, TISSUE_COLOR[tissue]])


def _composite(fractions):
    """The three fractions as one colour image on the surface colour."""
    colors = np.array([matplotlib.colors.to_rgb(TISSUE_COLOR[t]) for t in frontend.TISSUES])
    rgb = fractions @ colors
    return rgb + (1 - fractions.sum(-1, keepdims=True)) * np.array(matplotlib.colors.to_rgb(SURFACE))


def _sides(ax):
    """Mark the sides of the head in a panel."""
    for x, side in ((0.035, "L"), (0.965, "R")):
        ax.text(x, 0.035, side, transform=ax.transAxes, ha="center", va="bottom", fontsize=9,
                color="#ffffff", weight="bold")


def figure_frontend(d, path):
    """What ``frontend.segment`` returns: the whole slice above, a detail below."""
    box, zoom = _box(d), _zoom(d)
    scan_range = {"cmap": "gray", "vmin": 0, "vmax": 1.25 * np.median(d["scan"][d["label"] > 2.9])}
    panels = [
        ("T1w scan", "the input, shown on the working grid", d["scan"], scan_range),
        ('seg["corrected"]', "denoised and bias-corrected", d["corrected"], _gray(d["corrected"], d)),
        ('seg["label"]', "the label map", d["label"], {"cmap": "gray", "vmin": 0, "vmax": 3}),
    ] + [
        (f'seg["fractions"][..., {n}]', f"fraction of {t}", d["fractions"][..., n],
         {"cmap": _tissue_cmap(t), "vmin": 0, "vmax": 1})
        for n, t in enumerate(frontend.TISSUES)
    ]
    aspect = (box[1].stop - box[1].start) / (box[0].stop - box[0].start)
    width = 2.6
    fig, axes = plt.subplots(
        2, len(panels), figsize=(width * len(panels), width * (aspect + 1) + 1.85),
        gridspec_kw={"height_ratios": [aspect, 1], "hspace": 0.05, "wspace": 0.06},
        facecolor=SURFACE)
    fig.subplots_adjust(left=0.035, right=0.995, top=0.865, bottom=0.125)
    images = []
    for col, (title, sub, image, style) in enumerate(panels):
        top, low = axes[0, col], axes[1, col]
        _show(top, image, box, **style)
        images.append(_show(low, image, zoom, **style))
        code = title.startswith("seg")
        top.text(0.5, 1.115, title, transform=top.transAxes, ha="center", va="bottom", color=INK,
                 fontsize=9.6 if code else 12, family="monospace" if code else None,
                 weight=None if code else "bold")
        top.text(0.5, 1.025, sub, transform=top.transAxes, ha="center", va="bottom",
                 fontsize=10.5, color=INK_SOFT)
        top.add_patch(Rectangle((zoom[0].start - box[0].start - 0.5, zoom[1].start - box[1].start - 0.5),
                                zoom[0].stop - zoom[0].start, zoom[1].stop - zoom[1].start,
                                fill=False, edgecolor="#ffd23f" if col < 3 else INK, linewidth=1.3))
        for ax in (top, low):
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_edgecolor("#c8c8c8")
                spine.set_linewidth(0.8)
    _sides(axes[0, 0])
    axes[0, 0].set_ylabel("axial slice", fontsize=10.5, color=INK_SOFT)
    axes[1, 0].set_ylabel(f"detail, {ZOOM_MM[2]:.0f} mm wide", fontsize=10.5, color=INK_SOFT)

    def bar(col, ticks, labels, text):
        pos = axes[1, col].get_position()
        cax = fig.add_axes([pos.x0 + 0.08 * pos.width, 0.088, 0.84 * pos.width, 0.014])
        cb = fig.colorbar(images[col], cax=cax, orientation="horizontal")
        cb.set_ticks(ticks)
        cb.set_ticklabels(labels)
        cb.ax.tick_params(labelsize=9, colors=INK_SOFT, length=2.5)
        cb.outline.set_edgecolor("#c8c8c8")
        cb.ax.set_xlabel(text, fontsize=9.5, color=INK_SOFT, labelpad=3)

    bar(2, [0, 1, 2, 3], ["0", "CSF", "GM", "WM"], "in between: partial volumes")
    for col, tissue in enumerate(frontend.TISSUES, start=3):
        bar(col, [0, 0.5, 1], ["0", "0.5", "1"], f"share of {tissue} in the voxel")
    means = ",  ".join(f"{t} {v:.0f}" for t, v in zip(frontend.TISSUES, d["means"]))
    left = axes[1, 0].get_position().x0
    fig.text(left, 0.085, 'seg["means"]', fontsize=9.6, color=INK, family="monospace", va="center")
    fig.text(left, 0.040, f"intensity of the pure tissues in the corrected\nscan, from the mixture of "
             f"Gaussians:\n{means}", fontsize=9.5, color=INK_SOFT, va="center", linespacing=1.35)
    fig.suptitle("What t1sim.frontend.segment() makes of a real scan", fontsize=15, color=INK,
                 x=0.035, ha="left", y=0.985, weight="bold")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def figure_scheme(d, path):
    """From the scan to the simulated image: the products, and the steps between them."""
    box = _box(d)
    aspect = (box[1].stop - box[1].start) / (box[0].stop - box[0].start)
    tw = 2.15                      # width of an image
    th = tw * aspect
    sw = 2.25                      # width of a step
    gap = 0.34                     # between an image and a step
    left = 0.50
    pitch = tw + sw + 2 * gap
    width = left * 2 + 4 * tw + 3 * (sw + 2 * gap)
    head = 1.78                    # room above the images of a row: titles of row and images
    foot = 0.50                    # room below them
    foot_low = 1.62                # and below those of the second row, for the cleanup
    between = 1.25                 # between the two rows
    y_low = 0.22 + foot_low
    y_top = y_low + th + head + between + foot
    height = y_top + th + head + 0.22
    fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.axis("off")

    def panel(y, foot, face, edge, title, note):
        ax.add_patch(FancyBboxPatch((0.16, y - foot), width - 0.32, th + head + foot,
                                    boxstyle="round,pad=0,rounding_size=0.22",
                                    facecolor=face, edgecolor=edge, linewidth=1.2, zorder=0))
        ax.text(0.46, y + th + head - 0.36, title, fontsize=14.5, weight="bold", color=INK, va="center")
        ax.text(0.46, y + th + head - 0.70, note, fontsize=10.5, color=INK_SOFT, va="center")

    def image(col, y, data, title, sub, **style):
        x = left + col * pitch
        if data is not None:
            _show(ax, data, box, extent=(x, x + tw, y, y + th), zorder=2, **style)
            ax.add_patch(Rectangle((x, y), tw, th, fill=False, edgecolor="#b0b0b0", linewidth=0.8, zorder=3))
        else:
            ax.add_patch(Rectangle((x, y), tw, th, facecolor="#fafafa", edgecolor="#b0b0b0",
                                   linewidth=0.8, linestyle=(0, (4, 3)), zorder=2))
            ax.text(x + tw / 2, y + th / 2, "not made\nby t1sim yet", ha="center", va="center",
                    fontsize=10, color=INK_SOFT)
        ax.text(x + tw / 2, y + th + 0.36, title, ha="center", va="bottom", fontsize=12, weight="bold", color=INK)
        ax.text(x + tw / 2, y + th + 0.10, sub, ha="center", va="bottom", fontsize=9.6, color=INK_SOFT)
        return x

    def arrow(points, color=INK_SOFT):
        """A line through *points* that ends in an arrow head."""
        xs, ys = zip(*points)
        ax.plot(xs[:-1], ys[:-1], color=color, linewidth=1.5, solid_capstyle="butt", zorder=4)
        ax.annotate("", xy=points[-1], xytext=points[-2], zorder=4,
                    arrowprops=dict(arrowstyle="-|>", color=color, linewidth=1.5,
                                    shrinkA=0, shrinkB=0, mutation_scale=15))

    def step_box(x, top, w, head_text, body, done):
        lines = 1 + body.count("\n")
        h = 0.60 + 0.166 * lines
        ax.add_patch(FancyBboxPatch((x, top - h), w, h, boxstyle="round,pad=0,rounding_size=0.09",
                                    facecolor=DONE if done else OPEN,
                                    edgecolor=DONE_EDGE if done else OPEN_EDGE, linewidth=1.3,
                                    linestyle="-" if done else (0, (4, 3)), zorder=5))
        ax.text(x + w / 2, top - 0.13, head_text, ha="center", va="top", fontsize=10.8, weight="bold",
                color=INK, zorder=6)
        ax.text(x + w / 2, top - 0.47, body, ha="center", va="top", fontsize=9.1, color=INK,
                linespacing=1.3, zorder=6)
        return h

    def steps(col, y, items, done=True):
        """The steps between image *col* and the next, stacked, and the arrows into and out of them."""
        x = left + col * pitch + tw + gap
        mid = y + th / 2
        hs = [0.60 + 0.166 * (1 + body.count("\n")) for _, body in items]
        top = mid + (sum(hs) + 0.16 * (len(items) - 1)) / 2
        for (head_text, body), h in zip(items, hs):
            step_box(x, top, sw, head_text, body, done)
            top -= h + 0.16
        arrow([(x - gap + 0.05, mid), (x - 0.03, mid)])
        arrow([(x + sw + 0.03, mid), (x + sw + gap - 0.05, mid)])

    # ------------------------------------------------ first row: t1sim
    y = y_top
    panel(y, foot, "#f6faff", DONE_EDGE, "1  The tissue fractions of a real scan",
          "t1sim.frontend.segment()")
    scan_style = {"cmap": "gray", "vmin": 0, "vmax": 1.25 * np.median(d["scan"][d["label"] > 2.9])}
    x = image(0, y, d["scan"], "T1w scan", "a real scan of good quality", **scan_style)
    for xs, side in ((x + 0.10, "L"), (x + tw - 0.10, "R")):
        ax.text(xs, y + 0.07, side, ha="center", va="bottom", fontsize=9, color="#ffffff", weight="bold", zorder=4)
    steps(0, y, [("SANLM denoising", "at the resolution\nof the scan"),
                 ("Unified segmentation", "six-class templates,\naffine registration,\nbias field,\nmixture of Gaussians")])
    x = image(1, y, d["corrected"], "Corrected scan", "bias-corrected; brain mask outlined",
              **_gray(d["corrected"], d))
    ax.contour(np.linspace(x, x + tw, d["mask"][box].shape[0]),
               np.linspace(y, y + th, d["mask"][box].shape[1]),
               d["mask"][box].swapaxes(0, 1).astype(float), levels=[0.5], colors=["#ffd23f"],
               linewidths=1.1, zorder=3)
    steps(1, y, [("LAS", "local intensity\nnormalisation of the brain;\npure GM is the GM mean\nof the mixture")])
    image(2, y, d["label"], "Label map", "CSF 1, GM 2, WM 3, mixtures in between",
          cmap="gray", vmin=0, vmax=3)
    steps(2, y, [("Decomposition", "every voxel as a mixture\nof the two tissues\nnext to its label")])
    x = image(3, y, _composite(d["fractions"]), "Tissue fractions", "what the image is made from")
    for n, tissue in enumerate(frontend.TISSUES):
        ax.add_patch(Rectangle((x + 0.10 + n * 0.68, y - 0.33), 0.17, 0.17, color=TISSUE_COLOR[tissue], zorder=3))
        ax.text(x + 0.32 + n * 0.68, y - 0.245, tissue, fontsize=9.6, color=INK, va="center")

    # ------------------------------------------------ second row: the image and its ground truth
    y2 = y_low
    panel(y2, foot_low, "#f6faff", DONE_EDGE, "2  From the fractions to the simulated image and its ground truth",
          "t1sim.simulate: ground_truth() up to the synthesised image and the ground truth, "
          "acquire() for every simulated image")
    x0 = image(0, y2, _composite(d["fractions"]), "Tissue fractions", "as segmented, or changed")
    steps(0, y2, [("Synthesis", "brain: fractions x\ntissue intensities;\nrest of the head:\nthe corrected scan")])
    image(1, y2, d["sim_ideal"], "Synthesised T1w", "free of noise and bias field",
          **_gray(d["sim_ideal"], d))
    steps(1, y2, [("Acquisition", "contrast, bias field,\nvoxel size, movement,\nringing, noise")])
    x2 = image(2, y2, d["sim_image"], "Simulated T1w", "here with a bias field and noise",
               **_gray(d["sim_image"], d))
    x3 = image(3, y2, _composite(d["sim_fractions"]), "Ground truth", "the fractions after the cleanup")
    ax.text((x2 + tw + x3) / 2, y2 + th / 2, "belongs\nto", ha="center", va="center", fontsize=10.6,
            color=INK_SOFT, linespacing=1.3)
    for xa, xb in ((x2 + tw + 0.10, (x2 + tw + x3) / 2 - 0.42), (x3 - 0.10, (x2 + tw + x3) / 2 + 0.42)):
        ax.plot([xa, xb], [y2 + th / 2] * 2, color=INK_SOFT, linewidth=1.5, linestyle=(0, (2, 2)), zorder=4)

    # the cleanup: the same fractions, on their way to the ground truth
    line = y2 - 0.82
    bw = 6.3
    bx = (x0 + tw / 2 + x3 + tw / 2) / 2 - bw / 2
    h = step_box(bx, line + 0.465, bw, "Cleanup",
                 "blood vessels, dura and the rims of the ventricles look like grey matter:\n"
                 "they are corrected in the ground truth, and the image keeps them", done=True)
    del h
    arrow([(x0 + tw / 2, y2 - 0.04), (x0 + tw / 2, line), (bx - 0.03, line)])
    arrow([(bx + bw + 0.03, line), (x3 + tw / 2, line), (x3 + tw / 2, y2 - 0.04)])

    # ------------------------------------------------ from the first row to the second
    y_gap = y2 + th + head + between / 2
    start_x = left + tw / 2
    end_x = left + 3 * pitch + tw / 2
    arrow([(end_x, y - foot), (end_x, y_gap), (start_x, y_gap), (start_x, y2 + th + head - 0.02)])
    body = "cortical thickness: a band of a chosen thickness around the white matter   |   atrophy: a region deformed"
    bw = 8.1
    step_box(width / 2 - bw / 2, y_gap + 0.385, bw, "Change the anatomy (optional)", body, done=True)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def figure_examples(d, path):
    """What an acquisition and an atrophy do: a bias field on the whole slice, everything else on one detail."""
    zoom, box = _zoom(d), _box(d)
    gray = {"cmap": "gray", "vmin": 0, "vmax": 1.25 * float(d["means"][2])}     # one window for all images
    gm_before, gm_after = d["sim_fractions"][..., 1], d["ex_atrophy_fractions"][..., 1]
    change = LinearSegmentedColormap.from_list("change", [TISSUE_COLOR["CSF"], SURFACE, TISSUE_COLOR["GM"]])
    reached = abs(100 * float(d["ex_atrophy_reached"][0]))
    rows = [
        ("On the 0.5 mm grid of the ground truth", [
            ("Reference", "no artefact", d["sim_ideal"], gray),
            ("Noise", "SNR 12 in white matter", d["ex_noise"], gray),
            ("Contrast", "exponent 0.5: less", d["ex_contrast_low"], gray),
            ("Contrast", "exponent 1.5: more", d["ex_contrast_high"], gray)]),
        ("Acquired at 1 mm and brought back onto the 0.5 mm grid", [
            ("Voxel size", "1 mm, no other artefact", d["ex_resolution"], gray),
            ("Movement", "severity 2", d["ex_motion"], gray),
            ("Ringing", 'strength 2, type "notch"', d["ex_ringing"], gray),
            ("Ringing", 'strength 2, type "gibbs"', d["ex_gibbs"], gray)]),
        (f"Atrophy: the cortex of the left middle frontal gyrus thinned until it lost {reached:.0f}% of its volume", [
            ("Image before", "the reference", d["sim_ideal"], gray),
            ("Image after", "rendered from the deformed fractions", d["ex_atrophy"], gray),
            ("Ground truth after", "fraction of grey matter", gm_after,
             {"cmap": _tissue_cmap("GM"), "vmin": 0, "vmax": 1}),
            ("Change of the ground truth", "after minus before; region outlined", gm_after - gm_before,
             {"cmap": change, "vmin": -0.6, "vmax": 0.6})]),
    ]

    aspect = (box[1].stop - box[1].start) / (box[0].stop - box[0].start)
    fw, fh = 15.6, 13.3
    fig = plt.figure(figsize=(fw, fh), facecolor=SURFACE)

    def frame(ax):
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor("#c8c8c8")
            spine.set_linewidth(0.8)

    def titles(ax, title, sub, size=12.5):
        ax.text(0.5, 1.09, title, transform=ax.transAxes, ha="center", fontsize=size, weight="bold", color=INK)
        ax.text(0.5, 1.025, sub, transform=ax.transAxes, ha="center", fontsize=9.8, color=INK_SOFT)

    # the whole slice, twice: where the detail lies, and the one artefact a detail cannot show
    wide = 0.205
    tall = wide * aspect * fw / fh
    for n, (title, sub, image) in enumerate((
            ("Synthesised image", "no artefact; the detail is outlined", d["sim_ideal"]),
            ("Bias field", "60%, inverted, strength 3", d["ex_bias"]))):
        ax = fig.add_axes([0.012, 0.878 - tall - n * (tall + 0.075), wide, tall])
        _show(ax, image, box, **gray)
        if n == 0:
            ax.add_patch(Rectangle((zoom[0].start - box[0].start - 0.5, zoom[1].start - box[1].start - 0.5),
                                   zoom[0].stop - zoom[0].start, zoom[1].stop - zoom[1].start,
                                   fill=False, edgecolor="#ffd23f", linewidth=1.5))
            _sides(ax)
        frame(ax)
        titles(ax, title, sub)

    size = 0.182
    height = size * fw / fh
    x0, pitch_x, pitch_y = 0.238, size + 0.008, height + 0.094
    for r, (heading, panels) in enumerate(rows):
        top = 0.878 - r * pitch_y
        fig.text(x0, top + 0.043, heading, fontsize=11.5, color=INK, weight="bold", va="bottom")
        for c, (title, sub, image, style) in enumerate(panels):
            ax = fig.add_axes([x0 + c * pitch_x, top - height, size, height])
            shown = _show(ax, image, zoom, **style)
            if title.startswith("Change"):
                ax.contour(d["roi"][zoom].swapaxes(0, 1).astype(float), levels=[0.5], colors=[INK],
                           linewidths=0.9)
            frame(ax)
            titles(ax, title, sub, size=11.5)
            if style is not gray:
                cax = fig.add_axes([x0 + c * pitch_x + 0.15 * size, top - height - 0.022, 0.7 * size, 0.008])
                bar = fig.colorbar(shown, cax=cax, orientation="horizontal")
                ticks = [0, 0.5, 1] if style["vmin"] == 0 else [-0.6, 0, 0.6]
                bar.set_ticks(ticks)
                bar.set_ticklabels([f"{t:g}" if style["vmin"] == 0 else f"{t:+g}" if t else "0" for t in ticks])
                bar.ax.tick_params(labelsize=8.5, colors=INK_SOFT, length=2.5)
                bar.outline.set_edgecolor("#c8c8c8")
    fig.text(0.012, 0.878 - 2 * tall - 0.105, f"The detail is {ZOOM_MM[2]:.0f} mm wide.\nAll images share one grey scale,\n"
             "from 0 to 1.25 times the intensity\nof white matter.", fontsize=10, color=INK_SOFT, va="top",
             linespacing=1.4)
    fig.suptitle("What t1sim.simulate.acquire() and a change of the anatomy do", fontsize=15,
                 color=INK, x=0.012, ha="left", y=0.985, weight="bold")
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--scan", required=True, help="a T1w scan (NIfTI)")
    parser.add_argument("--cache", default=None, help="keep the slices in this .npz file")
    parser.add_argument("--out", default=str(Path(__file__).resolve().parent), help="where the figures go")
    args = parser.parse_args(argv)
    if args.cache and Path(args.cache).is_file():
        with np.load(args.cache) as cached:
            data = {key: cached[key] for key in cached.files}
    else:
        data = collect(args.scan)
        if args.cache:
            np.savez_compressed(args.cache, **data)
    out = Path(args.out)
    figure_frontend(data, out / "t1sim-frontend.png")
    figure_scheme(data, out / "t1sim-scheme.png")
    figure_examples(data, out / "t1sim-examples.png")
    print(f"wrote t1sim-frontend.png, t1sim-scheme.png and t1sim-examples.png to {out}")


if __name__ == "__main__":
    main()
