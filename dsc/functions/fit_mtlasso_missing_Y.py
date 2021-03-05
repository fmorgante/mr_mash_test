def fit_sparse_multi_task_lasso_missing_Y(X, Y, standardize, nfolds, B_init, grid_limits, grid_length):
  
  import numpy as np
  import mtlasso
  import pandas as pd
  import sklearn.model_selection as skms
  import time
  
  if grid_limits is None:
    grid = np.geomspace(.1, 1, grid_length) * X.shape[0]
  else:
    grid = np.linspace(start=grid_limits[1], stop=grid_limits[0], num=grid_length)
    
  t = time.process_time()
  
  Ym=np.ma.fix_invalid(Y, fill_value=0)
    
  cv_scores = mtlasso.lasso.sparse_multi_task_lasso_cv(
  X,
  Ym,
  cv=skms.KFold(n_splits=nfolds),
  lambda1=grid,
  lambda2=grid,
  max_iter=5000,
  standardize=standardize,
  verbose=False)
  cv_scores = pd.DataFrame(cv_scores)
  cv_scores.columns = ['fold', 'lambda1', 'lambda2', 'mse']
  lambda1, lambda2 = cv_scores.groupby(['lambda1', 'lambda2'])['mse'].agg(np.mean).idxmin()
  Bhat, B0hat = mtlasso.lasso.sparse_multi_task_lasso(X, Ym, init=B_init, standardize=standardize, lambda1=lambda1, lambda2=lambda2)
  
  elapsed_time = time.process_time() - t
  
  return Bhat, B0hat, elapsed_time
