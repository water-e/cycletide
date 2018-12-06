import numpy as np
from read_ts import *

ts = read_ts("freeport_flow.rdb")
with open("freeport_flow.csv","w") as g:
    for el in ts:
        val = "NA" if np.isnan(el.value) else el.value
        g.write("{}\n".format(val))
