# Random seed
seed = 42
n_jobs = 8

# List of options for each method
opts_list = {
  'udt' : {'min_leaf':20, 'seed':seed},
  'mdt' : {'trees':10, 'min_leaf':20, 'n_jobs':n_jobs, 'seed':seed},
  'aucell' : {'seed':seed},
  'wmean' : {'times':100, 'seed':seed},
  'wsum' : {'times':100, 'seed':seed},
  'ulm' : {},
  'mlm' : {},
  'viper' : {'pleiotropy':True},
  'gsva' : {'kcdf':True},
  'ora' : {'n_up':300, 'n_bottom':300, 'n_background':20000, 'seed':seed},
  'gsea' : {'times':100, 'seed':seed}
}
