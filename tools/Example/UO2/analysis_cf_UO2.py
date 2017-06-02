import h5py
import numpy as np
from scipy.linalg import eigvalsh

Ryd_eV = 13.605698066

f = h5py.File("glog.h5", 'r')
cf_bare = f["/Impurity_1/One_body_local"]
cf_norm = f["/Impurity_1/GA_La"]

cf_bare_eval = []
evals = eigvalsh(cf_bare[:2, 0:2])
for eval in evals:
  cf_bare_eval.append(eval)
evals = eigvalsh(cf_bare[4:6, 4:6])
for eval in evals:
  cf_bare_eval.append(eval)
cf_bare_eval.append(cf_bare[12,12])
cf_bare_eval = np.array(cf_bare_eval)

cf_norm_eval = []
evals = eigvalsh(cf_norm[:2, 0:2])
for eval in evals:
  cf_norm_eval.append(eval)
evals = eigvalsh(cf_norm[4:6, 4:6])
for eval in evals:
  cf_norm_eval.append(eval)
cf_norm_eval.append(cf_norm[12,12])
cf_norm_eval = np.array(cf_norm_eval)

print ' Bare crystal field (eV): \n', ((cf_bare_eval - cf_bare_eval[0])*Ryd_eV).real
print ' Renormalized crystal field (eV): \n', ((cf_norm_eval - cf_norm_eval[0])*Ryd_eV).real


