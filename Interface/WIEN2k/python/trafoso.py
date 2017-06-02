import numpy as np

def trafoso(L):
  '''
  trafoso provides transformation matrices from
  |L,1/2,mL,mS> (L=0,1,2,3, mS=-1/2,1/2) basis to
  basis |J,L,S,mJ>, J=L-1/2, L+1/2
  H. Watanabe 'Operator Methods in Ligand Field Theory'
  Prentice Hall, 1966, Table 1.8-1.
  ordering because of the convention used in WIEN is:
                     mS=1/2        mS=-1/2
                   -L .......L  -L ...... L     (2*(2L+1) columns)
          -(L-1/2)
             .
  J=L-1/2    .
             .
           (L-1/2)
           -L-1/2
             .
  J=L+1/2    .
             .
            L+1/2
  '''
  cf=np.zeros((2*(2*L+1),2*(2*L+1)))
  if L==0:
    cf[0,1]=1.0
    cf[1,0]=1.0
  else:
    k1=-1
    for ms in range(-1,2,2):
      ams=-ms/2.
      for ml in range(-L,L+1):
        k1=k1+1
        k2=-1
        for mj in range(-2*L+1,2*L,2):  # L-1/2 states
          amj=mj/2.
          k2=k2+1
          d=amj-ml-ams
          if abs(d)<0.0001:
            if ms==1:
              cf[k2,k1]=-np.sqrt((L+0.5+amj)/(2*L+1))
            else:
              cf[k2,k1]= np.sqrt((L+0.5-amj)/(2*L+1))
        for mj in range(-2*L-1,2*L+2,2): # L+1/2 states
          amj=mj/2.
          k2=k2+1
          d=amj-ml-ams
          if abs(d)<0.0001:
            if ms==1:
              cf[k2,k1]= np.sqrt((L+0.5-amj)/(2*L+1))
            else:
              cf[k2,k1]= np.sqrt((L+0.5+amj)/(2*L+1))
  return cf

if __name__=="__main__":
  print trafoso(2)
