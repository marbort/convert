import argparse
import os
import numpy
import pandas as pd
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def extract_data(path):
    data=pd.read_excel(path,sheet_name="NEDA", skiprows=17)
    return(data)

data=extract_data('/home/marco/OneDrive_UiO/2024/ML_Schlenk/EDA_summary.ods')
print(data)
