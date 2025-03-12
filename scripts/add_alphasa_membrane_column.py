#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

Time(1)
Query("ALTER TABLE `alphasa` ADD COLUMN `membrane` CHAR(1) NULL DEFAULT NULL AFTER `sec`, ADD INDEX `Membrane` (`membrane` ASC) VISIBLE")
Time(1)

print("\nDone!")
