import numpy as np

def get_data_array_from_string(str):
  '''
  Convert str to data list.
  e.g., str= 1 0
             2 1
  it outputs np.array([[1,0],[2,1]])

  '''
  res = []
  lines = str.split('\n')
  for line in lines:
    line_ = line.split()
    if len(line_) == 0:
      continue
    try:
      res.append(map(float, line_))
    except:
      print " Problematic line: ", line_
      pass
  return np.array(res)
