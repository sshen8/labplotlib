from labplotlib.blood import blood_cell_plot
import pandas as pd
import numpy as np
from labplotlib.blood.blood import CELL_TO_ID
import matplotlib.pyplot as plt

def test_plot():
    expression_across_cells = pd.Series(np.random.rand(len(CELL_TO_ID)), index=list(CELL_TO_ID.keys()), name="CD33")
    fig, ax = plt.subplots()
    blood_cell_plot(expression_across_cells, ax=ax)
    plt.close(fig)
