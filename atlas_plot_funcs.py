import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import ScalarFormatter
from collections import defaultdict
from coloraide import Color
from scipy.interpolate import interp1d
import os
import glob
from PIL import Image
import copy as cp
import pandas as pd


def formatPlots(
    show=False, figSize=(4, 3), usetex=False, dpi=600, Transparent_Background=False
):
    plt.style.use("fast")
    plt.style.use("ggplot")
    # mpl.use("pgf")
    # mpl.rcParams["pgf.texsystem"] = "pdflatex"
    if usetex:
        rc("text", usetex=True)
    else:
        rc("text", usetex=False)
    rc("font", family="sans-serif")
    mpl.rcParams["text.latex.preamble"] = " ".join(
        [
            r"\usepackage{amsmath}",
            r"\usepackage{siunitx}",
            r"\sisetup{detect-all}",
            r"\sisetup{math-micro=\text{µ},text-micro=µ}",
            r"\sisetup{math-ohm=\Omega}",
            r"\usepackage{amssymb}",
            r"\usepackage{sansmath}",
            r"\sansmath",
        ]
    )
    mpl.rcParams["lines.linewidth"] = 1
    mpl.rcParams["axes.titlesize"] = 11
    mpl.rcParams["axes.labelsize"] = 8
    mpl.rcParams["axes.linewidth"] = 0.8
    mpl.rcParams["axes.labelcolor"] = (0, 0, 0)
    mpl.rcParams["axes.facecolor"] = [1, 1, 1, 0.9]  # Make area inside axis white
    mpl.rcParams["axes.edgecolor"] = "black"
    mpl.rcParams["axes.formatter.use_mathtext"] = True  # offset is \times not e
    mpl.rcParams["grid.color"] = (0.92, 0.92, 0.92)
    mpl.rcParams["lines.markersize"] = 6
    mpl.rcParams["xtick.labelsize"] = 6
    mpl.rcParams["xtick.bottom"] = True
    mpl.rcParams["xtick.color"] = (0, 0, 0)
    mpl.rcParams["ytick.color"] = (0, 0, 0)
    mpl.rcParams["xtick.major.size"] = 2
    mpl.rcParams["ytick.labelsize"] = 6
    mpl.rcParams["ytick.left"] = True
    mpl.rcParams["ytick.major.size"] = 2
    mpl.rcParams["legend.fontsize"] = 8
    mpl.rcParams["figure.autolayout"] = True
    mpl.rcParams["figure.figsize"] = figSize
    mpl.rcParams["figure.facecolor"] = [1, 1, 1, 0]
    mpl.rcParams["figure.dpi"] = dpi
    mpl.rcParams["savefig.facecolor"] = [
        1,
        1,
        1,
        0,
    ]  # Save fig with tranparency outside
    # mpl.rcParams['savefig.transparent'] = Transparent_Background
    mpl.rcParams["savefig.dpi"] = dpi
    mpl.rcParams["path.simplify"] = True  # Set to true to increase rendering speed
    mpl.rcParams["path.simplify_threshold"] = 0.05  # Below this, points are removed
    if show:
        mpl.use("Qt5Agg")  # Interactive Backend
    else:
        mpl.use("Agg")  # Don't show plot


def plot_2d(
    infile,
    xrange,
    yrange,
    clabel,
    crange,
    saveDir,
    copies=2,
    title=None,
    half_sim=True,
    c_scale="linear",
    fig_size=None,
):
    """
    Plot and save an image from tonyplot with an improved axis.

    Parameters
    ----------
    infile : str
        The file path of the input image.
    xrange : tuple
        A tuple specifying the range of the x-axis.
    yrange : tuple
        A tuple specifying the range of the y-axis.
    clabel : str
        Label for the color bar.
    crange : tuple
        A tuple specifying the range of the color bar.
    saveDir : str
        The directory where the plot will be saved.
    copies : int, optional
        Number of copies to add horizontally. Defaults to 2.
    title : str, optional
        Title for the plot. Defaults to None.
    half_sim : bool, optional
        Whether to generate a half-symmetric plot. Defaults to True.
    c_scale : str, optional
        Scale for the color bar, either "linear" or "log". Defaults to "linear".
    fig_size : tuple, optional
        Size of the figure in inches (width, height). Defaults to None.

    Returns
    -------
    tuple
        A tuple containing the figure and the axis objects.
    """

    path, file = os.path.split(infile)
    file = file.strip(".png")
    cmap = rainbow30_cmap()

    with Image.open(infile) as im:
        bw = im.convert(mode="L")
        # Indexing is backwards as y,x based on row x collumn
        px = np.array(bw.getdata()).reshape(bw.height, bw.width)

        # find left edge
        left_idx = (
            np.argmin(
                np.sum(px[int(0.25 * bw.height) : int(0.75 * bw.height), :], axis=0)
            )
            + 1
        )
        # find right edge of frame
        right_frame = bw.width - np.argmin(
            np.flip(
                np.sum(
                    px[int(0.25 * bw.height) : int(0.75 * bw.height), left_idx:-1],
                    axis=0,
                )
            )
        )
        # find right edge of drawing
        right_idx = (
            right_frame
            - 12
            - np.argmin(
                np.flip(
                    np.sum(
                        px[
                            int(0.05 * bw.width) : int(0.95 * bw.width),
                            left_idx : right_frame - 10,
                        ],
                        axis=0,
                    )
                )
            )
        )
        # find top edge
        top_idx = (
            np.argmin(
                np.sum(px[:, int(0.25 * bw.width) : int(0.75 * bw.width)], axis=1)
            )
            + 1
        )
        # find bottom edge
        bottom_idx = top_idx + np.argmin(
            np.sum(px[top_idx:-1, int(0.25 * bw.width) : int(0.75 * bw.width)], axis=1)
        )
        # Crop colored image
        px = np.array(im.getdata()).reshape(im.height, im.width, -1)
        px = px[top_idx:bottom_idx, left_idx:right_idx]
        # Make full image by mirroring
        if half_sim:
            px = np.concatenate([px, np.flip(px, axis=1)], axis=1)
            copy = cp.copy(px)
        # add copies
        for i in range(copies - 1):
            px = np.concatenate([px, copy], axis=1)

        fig, axL = plt.subplots(nrows=1, ncols=1)
        if fig_size:
            fig.set_size_inches(fig_size[0], fig_size[1])
        axL.imshow(px, extent=[xrange[0], xrange[1] * 2 * copies, yrange[0], yrange[1]])
        # ADD A COLORBAR
        # get a normalizing object from the minimum to maximum value
        if c_scale == "log":
            normalize = mpl.colors.LogNorm(crange[0], crange[1])
        else:
            normalize = mpl.colors.Normalize(crange[0], crange[1])

        sm = ScalarMappable(cmap=cmap, norm=normalize)
        cb = fig.colorbar(sm, fraction=0.2, aspect=20, shrink=0.85)
        cb.set_label(clabel)  # Set the label
        sm.set_clim(crange[0], crange[1])  # set the colorbar limits
        cb.ax.tick_params(labelsize=8)  # set the font size
        cb.ax.yaxis.set_offset_position("right")

        axL.grid(False)
        if title:
            axL.set_title(title)
        print(saveDir + file + "_full.png")
        fig.savefig(
            saveDir + file + "_full.png", bbox_inches="tight", pad_inches=0.025, dpi=600
        )
    return fig, axL


def plot_cutline(
    file,
    name,
    x_label,
    y_label,
    saveDir,
    x_scale="linear",
    y_scale="linear",
    filetypes=["png"],
    xlim=None,
    ylim=None,
    yfactor=1,
):
    fig, axL = plt.subplots(nrows=1, ncols=1)
    axL.xaxis.set_major_locator(MaxNLocator(integer=True))
    axL.yaxis.set_major_locator(MaxNLocator(integer=True))

    data = pd.read_csv(file, sep=" ", header=None, names=["x", "y"], skiprows=4)

    axL.plot(data["x"], data["y"] * yfactor, linestyle="-", zorder=1)

    axL.set_xlabel(x_label)
    axL.set_ylabel(y_label)
    axL.set_xscale(x_scale)
    axL.set_yscale(y_scale)
    if xlim:
        axL.set_xlim(xlim)
    if ylim:
        axL.set_ylim(ylim)
    axL.grid(True, "major", lw=0.6)
    axL.margins(x=0, y=0)

    for ext in filetypes:
        plot_name = saveDir + name + "." + ext
        fig.savefig(plot_name, bbox_inches="tight", pad_inches=0.025, dpi=300)

    return fig, axL


def plot_IV(
    file,
    name,
    x_range,
    x_label,
    y_label,
    saveDir,
    x_scale="linear",
    y_scale="linear",
    filetypes=["png"],
    xlim=None,
    ylim=None,
):
    """
    Plot and save current-voltage (IV) characteristics from a given file.

    Parameters
    ----------
    file : str
        The file path containing the data to be plotted.
    name : str
        The name of the plot to be saved.
    x_range : tuple
        A tuple specifying the range of the x-axis.
    x_label : str
        Label for the x-axis.
    y_label : str
        Label for the y-axis.
    saveDir : str
        The directory where the plot will be saved.
    x_scale : str, optional
        Scale for the x-axis, either "linear" or "log". Defaults to "linear".
    y_scale : str, optional
        Scale for the y-axis, either "linear" or "log". Defaults to "linear".
    filetypes : list, optional
        A list of file types to save the plot (e.g., ["png", "pdf"]). Defaults to ["png"].
    xlim : tuple, optional
        Limits for the x-axis. Defaults to None.
    ylim : tuple, optional
        Limits for the y-axis. Defaults to None.

    Returns
    -------
    tuple
        A tuple containing the figure and the axis objects.
    """
    print("saved", name)
    fig, axL = plt.subplots(nrows=1, ncols=1)
    axL.xaxis.set_major_locator(MaxNLocator(integer=True))
    axL.yaxis.set_major_locator(MaxNLocator(integer=True))

    data = pd.read_csv(file, sep=" ", header=None, names=["V", "I"], skiprows=4)
    scale = 1 / x_range[1] * 1e8 * 1e-3  ## kA/cm^2

    axL.plot(data["V"], data["I"] * scale, linestyle="-", zorder=1)

    axL.set_xlabel(x_label)
    axL.set_ylabel(y_label)
    axL.set_xscale(x_scale)
    axL.set_yscale(y_scale)
    if xlim:
        axL.set_xlim(xlim)
    if ylim:
        axL.set_ylim(ylim)
    axL.grid(True, "major", lw=0.6)
    axL.margins(x=0, y=0)

    name, _ = name.split('.')
    for ext in filetypes:
        plot_name = saveDir + name + "." + ext
        fig.savefig(plot_name, bbox_inches="tight", pad_inches=0.025, dpi=300)

    return fig, axL



def plot_extract_voltage_sweep(
    files,
    saveName,
    ylim=[],
    xlim=[],
    yscale="linear",
    xscale="linear",
    xlabel="",
    ylabel="",
    yfactor=1,
    xfactor=1,
):
    """
    Plot and save the extracted sweep data from given files.

    Parameters
    ----------
    files : list
        A list of file paths containing the data to be plotted. The file name must be of the form "ID_contact_voltage_sweeptype".
    saveName : str
        The file name (with path) to save the plot.
    ylim : list, optional
        Limits for the y-axis. Defaults to [].
    xlim : list, optional
        Limits for the x-axis. Defaults to [].
    yscale : str, optional
        Scale for the y-axis, either "linear" or "log". Defaults to "linear".
    xscale : str, optional
        Scale for the x-axis, either "linear" or "log". Defaults to "linear".
    xlabel : str, optional
        Label for the x-axis. Defaults to "".
    ylabel : str, optional
        Label for the y-axis. Defaults to "".
    yfactor : float, optional
        Factor to scale the y-axis values. Defaults to 1.
    xfactor : float, optional
        Factor to scale the x-axis values. Defaults to 1.

    Returns
    -------
    tuple
        A tuple containing the figure and the axis objects.
    """
    data = defaultdict(list)
    for file in files:
        head, tail = os.path.split(file)
        ID, contact, V, sweep = tail.split("_")
        cut_data = pd.read_csv(file, sep=" ", header=None, names=["x", "y"], skiprows=4)
        data["V"].append(float(V))
        data["cut_data"].append(cut_data)
    data = pd.DataFrame.from_dict(data)
    data.sort_values("V", inplace=True, ignore_index=True)

    colors = Color.interpolate(
        [Color("rebeccapurple"), Color("orange")], method="natural", space="lch"
    )
    colors = [colors(x / len(data.V)).to_string(hex=True) for x in range(len(data.V))]
    cmap = mpl.colors.ListedColormap(colors)

    fig, axL = plt.subplots(nrows=1, ncols=1)
    for index, row in data.iterrows():
        axL.plot(
            row.cut_data.x * xfactor, row.cut_data.y * yfactor, color=colors[index]
        )

    axL.set_xscale(xscale)
    axL.set_yscale(yscale)
    axL.set_xlabel(xlabel)
    axL.set_ylabel(ylabel)
    if xlim:
        axL.set_xlim(xlim)
    if ylim:
        axL.set_ylim(ylim)

    # get a normalizing object from the minimum to maximum value
    # LogNorm for log scale, Normalize for linear
    normalize = mpl.colors.Normalize(data.V.max(), data.V.min())
    colormap = cmap  # get the colormap
    sm = ScalarMappable(cmap=colormap, norm=normalize)
    cb = plt.colorbar(sm)  # Add it to the plot
    cb.set_label("Voltage [V]")  # Set the label
    sm.set_clim(data.V.max(), data.V.min())  # set the colorbar limits
    cb.ax.tick_params(labelsize=8)

    fig.savefig(saveName, bbox_inches="tight", pad_inches=0.025, dpi=600)

    return fig, axL


def rainbow30_cmap():
    '''
        Returns the colormap to match atlas/tonyplot
    '''
    cmap = (
        np.array(
            [
                (255, 0, 255),
                (229, 0, 255),
                (203, 0, 255),
                (178, 0, 255),
                (152, 0, 255),
                (127, 0, 255),
                (102, 0, 255),
                (0, 0, 255),
                (0, 101, 255),
                (0, 127, 255),
                (0, 153, 255),
                (0, 178, 255),
                (0, 203, 255),
                (0, 229, 255),
                (0, 255, 255),
                (0, 255, 229),
                (0, 255, 203),
                (0, 255, 178),
                (0, 255, 153),
                (0, 255, 127),
                (50, 255, 0),
                (127, 255, 0),
                (178, 255, 0),
                (204, 255, 0),
                (229, 255, 0),
                (255, 255, 0),
                (255, 229, 0),
                (255, 204, 0),
                (255, 178, 0),
                (255, 153, 0),
                (255, 127, 0),
                (255, 102, 0),
                (255, 50, 0),
                (255, 0, 0),
            ]
        )
        / 255
    )
    cmap = mpl.colors.ListedColormap(cmap)  # get the colormap
    return cmap


def sci_fmt(x):
    """
        Takes a string, x given in 'e' format and formats in math scientific notation
        Note: Likely not needed with newer versions of python
    """
    a = r"{:.1g}".format(x)
    if "e" in a:
        a, b = a.split("e")
        a = int(a)
        b = int(b)
        if abs(b) >= 1:
            if a == 1:
                sci_str = r"$10^{" + f"{b:d}" + "}$"
            else:
                sci_str = r"$" + f"{a:.2g}" + r"\times 10^{" + f"{b:d}" + r"}$"
        else:
            sci_str = "$" + f"{a:.2g}" + "$"
    else:
        sci_str = a
    print(sci_str)
    return sci_str
