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
            if match:
                match = match.groupdict()
            else:
                warnings.warn(f"dropping {idx} because no regex match")
                continue
            ind = tuple(match.get(x, default.get(x)) for x in names)
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
        if split is None or split["split_on"] > sample_time:
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
