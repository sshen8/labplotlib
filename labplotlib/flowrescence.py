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
from sklearn.mixture import GaussianMixture

# TODO: https://github.com/matplotlib/matplotlib/issues/6321

def read_csv(folder, well_id=True):
    """
    subfolders in `folder` are "specimens". csv files in subfolders are "samples".
    returns concatenated dataframes where index is (specimen, sample)
    """
    regex_pattern = "[A-H]1?[0-9]_(.*)\.csv" if well_id else "(.*)\.csv"
    df = {}
    for folder_specimen in os.listdir(folder):
        subfolder = os.path.join(folder, folder_specimen)
        if not os.path.isdir(subfolder):
            continue
        for file in os.listdir(subfolder):
            sample_name = re.search(regex_pattern, file).group(1)
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
            ind = tuple(match[x] if match[x] else default.get(x, "") for x in names)
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
    return np.log10(np.maximum(series, 1))

def geomean(series, axis=None):
    return np.power(10, np.mean(log10(series), axis=axis))

def geomstd(series, axis=None):
    logstd = np.std(log10(series), axis=axis)
    logmean = np.mean(log10(series), axis=axis)
    return (
          geomean(series, axis=axis) - np.power(10, logmean - logstd),
        - geomean(series, axis=axis) + np.power(10, logmean + logstd),
        )

FLUOR_CHANNELS = (
    "data_PE-Texas Red-H", "data_PE-Texas Red-A",
    "data_FITC-H", "data_FITC-A",
    "data_mCherry-H", "data_mCherry-A",
    "data_APC-H", "data_APC-A",
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

def gmm(data, groupby):
    if type(data) is not pd.Series:
        raise ArgumentError(f"expecting pd.Series but got {type(data).__name__}")
    def _calc_gmm(subdata):
        vals = np.expand_dims(subdata, axis=1)
        gm = GaussianMixture(n_components=2, random_state=2020).fit(vals)
        # assume transduced cells have higher fluorescence than untransduced cells
        transd_idx = gm.means_.squeeze().argmax()
        untransd_idx = gm.means_.squeeze().argmin()
        return {
            "transd": gm.weights_[transd_idx].item(),
            "mean": gm.means_[transd_idx].item(),
            "var": gm.covariances_[transd_idx].item(),
            "mean_untransd": gm.means_[untransd_idx].item(),
            "var_untransd": gm.covariances_[untransd_idx].item(),
        }
    # TODO: consider untransduced cells together
    return data.groupby(groupby).apply(_calc_gmm).unstack()

def norm_day0(df):
    return

def _plot_hist_default_cmap():
    return plt.cm.winter

def _plot_hist_default_norm(df, hue_level):
    if np.issubdtype(df.index.dtypes[hue_level], np.number):
        level_values = df.index.unique(hue_level)
        return plt.Normalize(vmin=level_values.min(), vmax=level_values.max())
    else: # categorical data
        return lambda x: x

def _plot_process_color_args(df, hue_level, cmap, norm):
    if cmap is None:
        cmap = _plot_hist_default_cmap()
    elif isinstance(cmap, dict):
        cmap_dict = dict(**cmap)
        cmap = lambda x: cmap_dict[x]
    if norm is None:
        norm = _plot_hist_default_norm(df, hue_level)    
    return cmap, norm

def _plot_process_args(df, x_level=None, hue_level=None, data_column=None):
    assert hue_level is not None
    if data_column:
        if data_column[:5] != "data_":
            data_column = f"data_{data_column}"
        assert data_column in FLUOR_CHANNELS
    for level in df.index.names:
        if level not in (x_level, hue_level, "fname_sample", "fname_specimen", None) and len(df.index.unique(level)) > 1:
            warnings.warn(f"combining experiment conditions: {level} = {df.index.unique(level)}")
    return df, x_level, hue_level, data_column

def plot_hist(ax, df, hue_level, data_column, bins=40, cmap=None, norm=None, hide_ticks=True):
    cmap, norm = _plot_process_color_args(df, hue_level, cmap, norm)
    df, _, hue_level, data_column = _plot_process_args(df, None, hue_level, data_column)
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

def plot_hist_legend(ax, df=None, hue_level=None, cmap=None, norm=None, title=None, captions=None, sort_key=None, alpha=0.4):
    cmap, norm = _plot_process_color_args(df, hue_level, cmap, norm)
    if sort_key is None:
        sort_key = norm
    elif isinstance(sort_key, list):
        sort_key = sort_key.index
    vals = df.index.unique(hue_level)
    if sort_key:
        vals = sorted(vals, key=sort_key)
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
    for hue_val in vals:
        ax.fill_between(x=[], y1=[], y2=[], alpha=alpha, label=captions_func(hue_val), color=cmap(norm(hue_val)))
    ax.axis("off")
    legend = ax.legend(loc="center")
    adjust_legend_subtitles(legend)

def plot_bar(ax, df, x_level, hue_level, cmap=None, norm=None, mean_column="mean_log10", std_column="std_log10"):
    cmap, norm = _plot_process_color_args(df, hue_level, cmap, norm)
    df, x_level, hue_level, _ = _plot_process_args(df, x_level, hue_level, None)
    x_pos, x_width = np.linspace(-0.45, 0.45, num=len(df.index.unique(hue_level)) + 1, retstep=True)
    x_pos = x_pos + 0.5 * x_width
    x_pos = x_pos[:-1]
    x_width *= 0.9
    x_ticks = []
    for i, (x_val, group) in enumerate(df.groupby(x_level)):
        ax.bar(x_pos + i, np.power(10, group[mean_column]), width=x_width, color=[cmap(norm(x)) for x in group.index.get_level_values(hue_level)])
        if std_column:
            ax.errorbar(x=x_pos + i, y=np.power(10, group[mean_column]),
                yerr=(
                    (- np.power(10, group[mean_column] - group[std_column]) + np.power(10, group[mean_column])),
                    (  np.power(10, group[mean_column] + group[std_column]) - np.power(10, group[mean_column]))),
                color="k", fmt="none",
                elinewidth=0.4, capsize=2, capthick=0.4)
        x_ticks.append(x_val)
    ax.set_yscale("log")
    return x_ticks

def plot_bar_legend(*args, **kwargs):
    return plot_hist_legend(*args, **kwargs, alpha=1)

def counts(df, groupby=None):
    if groupby is None:
        groupby = df.index.names[:-1]
    return df.iloc[:, 0].groupby(groupby).count().rename("counts")
