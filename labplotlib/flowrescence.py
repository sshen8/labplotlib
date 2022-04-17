"""
For flow fluorescence data
"""
from argparse import ArgumentError
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from shapely.geometry import Polygon, Point, LineString

def read_csv(timepoints, regex):
    """
    `timepoints`: dict of folders containing csv data mapped to datetime.datetime of data collection
    """
    return pd.concat(
        {(sample_time, re.search(regex, file)["sample_name"]): pd.read_csv(os.path.join(folder_timepoint, file)).add_prefix("data_")
                for folder_timepoint, sample_time in timepoints.items()
                for file in os.listdir(folder_timepoint)},
        names=("sample_time", "sample_name")
    )

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
