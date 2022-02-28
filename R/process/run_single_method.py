import decoupler as dc
import pandas as pd
import numpy as np
from methods_params import opts_list
import sys


# Parse args
args = sys.argv[1:]
mat_path = args[0]
net_path = args[1]
method = args[2]

# Methods dictionary
methods = {
  'aucell' : dc.run_aucell,
  'udt' : dc.run_udt,
  'mdt' : dc.run_mdt,
  'wmean' : dc.run_wmean,
  'ulm' : dc.run_ulm,
  'mlm' : dc.run_mlm,
  'wsum' : dc.run_wsum,
  'viper' : dc.run_viper,
  'gsva' : dc.run_gsva,
  'ora' : dc.run_ora,
  'gsea' : dc.run_gsea
}

# Open
mat = pd.read_csv(mat_path).T
mat.columns = (np.arange(mat.shape[1]) + 1).astype(str)
mat.index = (np.arange(mat.shape[0]) + 1).astype(str)
net = pd.read_csv(net_path)
net['source'] = net['source'].astype(str)
net['target'] = net['target'].astype(str)
net.columns = ['source','target','weight']

# Run
f = methods[method]
params = opts_list[method]
res = f(mat, net, **params)
