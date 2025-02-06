#!/usr/bin/env python3

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def apply_mpl_settings():
    plt.rc("text", usetex=False)
    plt.rc("font", size=10, family="cmu serif")
    plt.rc("mathtext", fontset="cm")
    plt.rc("axes", unicode_minus=False)
    plt.rc("xtick", direction="in", top=True)
    plt.rc("ytick", direction="in", right=True)
    plt.ioff()
apply_mpl_settings()


def savefig(outname, dpi=300, overwrite=True):
    outpath = Path(outname)
    if outpath.exists() and not overwrite:
        print(f"Figure exists, continuing: {outpath}")
    else:
        for ext in ("pdf", "png"):
            filen = str(outpath) + f".{ext}"
            plt.savefig(filen, dpi=dpi)
        print(f"-- Figure saved to: {outpath}")
        plt.close("all")


def set_minor_ticks(ax, x=True, y=True):
    if x:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    if y:
        ax.yaxis.set_minor_locator(AutoMinorLocator())


def set_grid(ax):
    ax.grid(linestyle="dashed", color="0.3", linewidth=0.3)


