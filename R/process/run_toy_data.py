import decoupler as dc
import pandas as pd

# Read data
mat = pd.read_csv('data/raw/toy_mat.csv', index_col=0).T
net = pd.read_csv('data/raw/toy_net.csv')
net.columns = ['source', 'target', 'weight']

# Args
seed = 42
n_jobs = 8
opts_list = {
  'udt' : {'min_leaf':3, 'seed':seed},
  'mdt' : {'trees':10, 'min_leaf':3, 'n_jobs':n_jobs, 'seed':seed},
  'aucell' : {'n_up':3, 'seed':seed},
  'wmean' : {'times':100, 'seed':seed},
  'wsum' : {'times':100, 'seed':seed},
  'ulm' : {},
  'mlm' : {},
  'viper' : {'pleiotropy':True},
  'gsva' : {'kcdf':True},
  'ora' : {'n_up':3, 'n_bottom':0, 'n_background':20000, 'seed':seed},
  'gsea' : {'times':100, 'seed':seed}
}

# Run
res = dc.decouple(mat, net, min_n=0, verbose=False, methods='all', args=opts_list)

# Melt and sort
res = dc.melt(res).sort_values(['method', 'source', 'sample'])

# Store
res.to_csv('data/prc/toy_acts_Python.csv', index=None)
