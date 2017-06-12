def create_GOMP_file():
  '''
  Create GOMP.INP file based on GMPI_X.INP for openMP execution.
  '''
  import numpy as np
  with open('GMPI_0.INP', 'r') as f:
    lines = f.readlines()
    line = (lines[0]).split()
    num_procs = int(line[1])
    if num_procs == 1:
      return
  idx_vector = []
  for iproc in range(num_procs):
    with open('GMPI_'+str(iproc)+'.INP', 'r') as f:
      lines = f.readlines()
      line = (lines[0]).split()
      n_vec = int(line[4])
      for iline in range(n_vec):
        idx_vector.append(lines[1+iline])
  nk_list = np.array([int((line.split())[1]) for line in idx_vector])
  idx = nk_list.argsort()
  idx_vector_sort = []
  for i in idx[::-1]:
    idx_vector_sort.append(idx_vector[i])
  with open('GOMP.INP', 'w') as f:
    f.write(' 0 '+str(num_procs)+' 0  T  '+str(len(idx_vector))+'\n')
    for line in idx_vector_sort:
      f.write(line)
    f.write(" MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR,VECTORS")

if __name__=='__main__':
  create_GOMP_file()
