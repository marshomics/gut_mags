"""
plotting_theme.py -- one place that controls how every Python figure looks and
how it is saved. Guarantees the two things the brief requires:

  * every figure is written as BOTH .png (raster, dpi from config) and .svg;
  * SVG text stays as text (svg.fonttype = 'none'), so labels are editable in
    Illustrator / Inkscape rather than being converted to vector outlines.

Import and call `apply_theme(cfg)` at the top of each figure script, then
`save(fig, path_without_extension, cfg)` to emit both formats.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def apply_theme(cfg):
    f = cfg["figures"]
    fam = f.get("font_family", "Arial")
    plt.rcParams.update({
        "svg.fonttype": "none",          # editable text in SVG
        "pdf.fonttype": 42,              # editable text if PDF ever used
        "font.family": "sans-serif",     # resolve via the fallback list below
        "font.sans-serif": [fam, "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": f.get("base_font_size", 7),
        "axes.titlesize": f.get("base_font_size", 7) + 1,
        "axes.labelsize": f.get("base_font_size", 7),
        "xtick.labelsize": f.get("base_font_size", 7) - 1,
        "ytick.labelsize": f.get("base_font_size", 7) - 1,
        "legend.fontsize": f.get("base_font_size", 7) - 1,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "figure.dpi": 120,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def niche_palette(cfg):
    return dict(cfg["figures"]["palette"])


def save(fig, path_no_ext, cfg):
    dpi = cfg["figures"].get("dpi", 400)
    for fmt in cfg["figures"].get("formats", ["png", "svg"]):
        fig.savefig(f"{path_no_ext}.{fmt}", format=fmt,
                    dpi=dpi if fmt == "png" else None)
    plt.close(fig)


# millimetre -> inch helper for journal-sized panels
def mm(x):
    return x / 25.4
