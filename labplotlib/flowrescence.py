"""
For flow fluorescence data
"""
from argparse import ArgumentError
import pandas as pd
import numpy as np
import os
import re
from shapely.geometry import Polygon, Point, LineString
import warnings
from matplotlib import pyplot as plt
from seaborn.utils import adjust_legend_subtitles

def read_csv(folder):
    """
    subfolders in `folder` are "specimens". csv files in subfolders are "samples".
    returns concatenated dataframes where index is (specimen, sample)
    """
    df = {}
    for folder_specimen in os.listdir(folder):
        subfolder = os.path.join(folder, folder_specimen)
        if not os.path.isdir(subfolder):
            continue
        for file in os.listdir(subfolder):
            sample_name = re.search("[A-H]1?[0-9]_(.*)\.csv", file).group(1)
            df[(folder_specimen, sample_name)] = pd.read_csv(os.path.join(subfolder, file)).add_prefix("data_")
    return pd.concat(df, names=("fname_specimen", "fname_sample"))

def parse_index(df, level, regex=None, apply=None, names=None, replace=None, default=None, droplevel=False):
    """
    maps `level` to new level `names` using mapping:
        - `regex` where group names should correspond to new levels
        - `apply` is a function, output can be tuple
        - `replace` is a dict, output can be tuple
    """
    if names is None:
        names = tuple()
    elif isinstance(names, str):
        names = (names,)
    if default is None and regex is not None:
        default = dict()
    newdf = {}
    for idx, group in df.groupby(level):
        if regex is not None:
            match = re.search(regex, idx)
            if not match:
                warnings.warn(f"dropping {idx} because no regex match")
                continue
            ind = tuple(match[x] if match[x] else default.get(x) for x in names)
        elif apply is not None:
            ind = apply(idx)
        elif replace is not None:
            ind = replace.get(idx, default)
        toappend = group
        if droplevel:
            toappend = toappend.droplevel(level=level)
        if ind in newdf:
            newdf[ind] = pd.concat([newdf[ind], toappend])
        else:
            newdf[ind] = toappend
    return pd.concat(newdf, names=names)

def split_sample(df, splits, level_time="sample_time", level_sample="fname_sample"):
    newdf = {}
    for (fname_sample, sample_time), group in df.groupby([level_sample, level_time]):
        split = splits.get(fname_sample)
        if split is None or split["split_on"] < sample_time:
            newdf[(fname_sample, sample_time)] = group.droplevel([level_sample, level_time])
            continue
        # TODO, make this recursive i.e. look what else this should get split into further down the line
        for split_to in split["split_to"]:
            newdf[(split_to, sample_time)] = group.droplevel([level_sample, level_time])
    return pd.concat(newdf, names=(level_sample, level_time))

def squeeze_index(df):
    return df.droplevel([idx for idx in df.index.names if df.index.unique(idx).size <= 1])

def log10(series):
    return np.log10(series.clip(lower=1))

FLUOR_CHANNELS = (
    "data_PE-Texas Red-H", "data_PE-Texas Red-A",
    "data_FITC-H", "data_FITC-A",
    )

def gate2d(df, boundary, interior=True):
    if isinstance(boundary, Polygon):
        return df.apply(Point, axis=1).apply(boundary.covers) == interior
    elif isinstance(boundary, LineString):
        (x1, y1), (x2, y2) = boundary.coords
        x, y = df.iloc[:,  0], df.iloc[:, 1]
        if x.name in FLUOR_CHANNELS:
            x1, x2 = np.log10((x1, x2))
            x = log10(x)
        if y.name in FLUOR_CHANNELS:
            y1, y2 = np.log10((y1, y2))
            y = log10(y)
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
        return (y <= m * x + b) == interior
    else:
        raise ArgumentError()

def norm_day0(df):
    return

def _plot_hist_default_cmap():
    return plt.cm.winter

def _plot_hist_default_norm(df, hue_level):
    level_values = df.index.unique(hue_level)
    return plt.Normalize(vmin=level_values.min(), vmax=level_values.max())

def plot_hist(ax, df, hue_level, data_column, bins=40, cmap=None, norm=None, hide_ticks=True):
    if cmap is None:
        cmap = _plot_hist_default_cmap()
    if norm is None:
        norm = _plot_hist_default_norm(df, hue_level)
    if data_column[:5] != "data_":
        data_column = f"data_{data_column}"
    assert data_column in FLUOR_CHANNELS
    df = squeeze_index(df)
    for level in df.index.names:
        if level not in (hue_level, "fname_sample", "fname_specimen", None):
            warnings.warn(f"combining experiment conditions: {level} = {df.index.unique(level)}")
    if isinstance(bins, int):
        bins = np.logspace(np.log10(df[data_column].min()), np.log10(df[data_column].max()), bins + 1)
    log_bins = np.log10(bins)
    bins_ctr = np.power(10, log_bins[:-1] + 0.5 * np.diff(log_bins))
    for hue_val, group in df.groupby(hue_level):
        bin_counts, _ = np.histogram(group[data_column], bins=bins)
        bin_mass = bin_counts / bin_counts.sum()
        ax.fill_between(x=bins_ctr, y1=0, y2=bin_mass, alpha=0.4, color=cmap(norm(hue_val)))
    ax.set_ylim(bottom=0)
    if hide_ticks:
        ax.set_yticks([])
    ax.set_xscale("log")

def plot_hist_legend(ax, df=None, hue_level=None, cmap=None, norm=None, title=None, captions=None):
    if cmap is None:
        cmap = _plot_hist_default_cmap()
    if norm is None:
        norm = _plot_hist_default_norm(df, hue_level)
    if title:
        ax.plot([], [], label=title, visible=False)
    if captions is None:
        captions_func = lambda x: x
    elif isinstance(captions, dict):
        captions_func = lambda x: captions.get(x, x)
    elif callable(captions):
        captions_func = captions
    else:
        raise ArgumentError()
    for hue_val in sorted(df.index.unique(hue_level), key=norm):
        ax.fill_between(x=[], y1=[], y2=[], alpha=0.4, label=captions_func(hue_val), color=cmap(norm(hue_val)))
    ax.axis("off")
    legend = ax.legend(loc="center")
    adjust_legend_subtitles(legend)

def counts(df, groupby=None):
    if groupby is None:
        groupby = df.index.names[:-1]
    return df.iloc[:, 0].groupby(groupby).count().rename("counts")
