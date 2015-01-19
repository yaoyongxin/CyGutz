


class material_class():
   'Class for a material in Wien2k+GA'

   # rb = 0.529177249      # Bohr radius A
   # convE = 13.605698066  # 1 Ry = 13.605698066 eV

   def __init__(self, folder): # %%%
       import numpy as np
       import os

       self.rb = 0.529177249      # Bohr radius A
       self.convE = 13.605698066  # 1 Ry = 13.605698066 eV

       # Name files where data is taken
       endfolder = folder.split("/")[-1]
       self.struct = folder + "/" + endfolder + ".struct"
       self.outputkgen = folder + "/" + endfolder + ".outputkgen"
       self.indmfl = folder + "/" + endfolder + ".indmfl"
       self.output = folder + "/" + endfolder + ".outputdmf1.0"
       self.scf = folder + "/" + endfolder + ".scf"
       self.GL = folder + "/" + "GL.INP"
       self.GLX = folder + "/" + "GLX.OUT"

       # Read default name material as first line structure file
       with open(self.struct, 'r') as f:
          name_material = f.readline().split()[0]
          self.name_material = name_material

       # Set properties of the material by reading data from LDA files
       self.set_volumecell()
       self.set_structure()
       self.n_atoms = np.sum(np.array(self.mult_unequiv_atoms))
       self.volume = self.volumecell/(self.n_atoms *1.0)

       # If initialization already done read corresponding information
       if os.path.exists(self.indmfl): 
          self.set_corr_atoms_init() # reads from self.indmf atoms information


   # **********************************************************
   # MAIN METHODS (ACTIONS THAT ONE CAN TAKE AS A USER) # %%%
   # **********************************************************
   
   def print_info(self):
      import os
      # Read information from DFT calculation
      print "Material:", self.name_material
      print
      print "Total number of atoms:", self.n_atoms
      print "Number unequivalent atoms:", self.n_unequiv_atoms
      if not os.path.exists(self.indmfl):
         print
         print "DFT+GA calculation not still initialized"
         print
         print "List groups unequivalent atoms:"
         for i,a in enumerate(self.unequiv_atoms):
            print i, a, "; multiplicity =", self.mult_unequiv_atoms[i]
      # Read information from ga_init_dmft.py initialization (if done)
      if os.path.exists(self.indmfl):
         print
         print "DFT+GA calculation already initialized"
         print
         print "List groups unequivalent atoms:"
         c = 0
         for i,a in enumerate(self.unequiv_atoms):
            s = self.atoms_corrstatus[i]
            if s == "uncorrelated":
               print i, a, "("+s+")", "; multiplicity =", self.mult_unequiv_atoms[i]
            if s == "correlated":
               atom = self.corr_atoms_instances[c]
               print i, atom.name, "("+ s + " ["+str(c)+"] " + ": " + atom.orb_label+")", \
                   "; multiplicity =", atom.mult, "; $\Sigma$ =", atom.ind_sigma_labels
               c = c+1

   def analyze_GA(self): # Analyze GA calculation
      self.set_energycell()
      self.energy = self.energycell/(self.n_atoms *1.0)
      self.set_corr_atoms_instances_results()
   #
   # **********************************************************

   #----------------------------------------
   # METHODS TO SET INFORMATION GIVEN BY USER
   # ABOUT GA SOLVER AND CREATE INPUT FILES # %%%

   def set_parametrization(self, parametrization = "Slater"):
      # Slater, Kanamori, Fock
      self.parametrization = parametrization

   def set_intParameters(self, U_list,J_list):
      # Just provide list for all interaction-parameters for all correlated atoms
      assert(len(U_list) == self.n_corr_atoms)
      assert(len(J_list) == self.n_corr_atoms)
      self.U_list = U_list
      self.J_list = J_list

   def set_localHilbert(self, Hilbert_list):
      # List of lists: for each k give portion of Hilbert space considere
      assert(len(Hilbert_list) == self.n_corr_atoms)
      self.Hilbert_list = Hilbert_list

   def set_nf0(self, nf0_list):
      # List of occupations $n_f^0$ defining fix double-counting potentials
      assert(len(nf0_list) == self.n_corr_atoms)
      self.nf0_list = nf0_list

   #----------------------------------------
   # METHODS TO READ AUTOMATICALLY INFORMATION
   # FROM DFT CALCULATION # %%%

   def set_name_material(self,name_material):
      self.name_material = name_material # (e.g., "FeO")

   def set_structure(self): # names unequivalent atoms and respective multiplicities
       import numpy as np
       import os
       #
       # Name and number of unequivalent atoms #
       command = self.grep_command("RMT=", self.struct, grepfile="ff")
       exec command
       unequiv_atoms=[]
       for line in ff:
           line = line[:-1]
           line = line.split()
           unequiv_atoms.append(line[0])
       n_unequiv_atoms = len(unequiv_atoms)
       #
       self.unequiv_atoms = unequiv_atoms
       self.n_unequiv_atoms = n_unequiv_atoms
       #
       # Multiplicity of each group of equivalent atoms
       command = self.grep_command("MULT=", self.struct, grepfile="ff")
       exec command
       #
       mult_unequiv_atoms=[]
       for line in ff:
           line = line[:-1]
           line = line.split()
           mult_unequiv_atoms.append(int(line[1]))
       # 
       self.mult_unequiv_atoms = mult_unequiv_atoms

   def set_volumecell(self): # vol unit cell (A^3)
       import numpy as np
       import os
       #
       command = self.grep_command("R1 = ", self.outputkgen, grepfile="f")
       exec command
       line = f.readline().split()
       R1 = np.array([float(line[2]),float(line[3]),float(line[4])]) * self.rb
       #
       command = self.grep_command("R2 = ", self.outputkgen, grepfile="f")
       exec command
       line = f.readline().split()
       R2 = np.array([float(line[2]),float(line[3]),float(line[4])]) * self.rb
       #
       command = self.grep_command("R3 = ", self.outputkgen, grepfile="f")
       exec command
       line = f.readline().split()
       R3 = np.array([float(line[2]),float(line[3]),float(line[4])]) * self.rb
       #
       M =  np.zeros(shape=(3,3))
       M[:,0]=R1
       M[:,1]=R2
       M[:,2]=R3
       # 
       self.volumecell = np.linalg.det(M)
 
   #----------------------------------------
   # METHODS TO READ AUTOMATICALLY INFORMATION GIVEN IN
   # GA INITIALIZATION FILES AFTER RUNNING ga_init_dmft.py # %%%

   def read_label_corr(self):  # reads labels correlated atoms from self.indmfl
      import os
      import numpy as np
      #
      iatom_list = []
      #
      command = self.grep_command("# iatom,", self.indmfl, grepfile="f")
      exec command
      for line in f:
         line = line[:-1]
         line = line.split()
         iatom_list.append(int(line[0])-1)
      #
      return iatom_list  # list all atomic labels

   def read_struct_sigma(self):  # reads structure sigma from self.indmfl
      import numpy as np
      #
      cix_num_list = []
      dim_list = []
      struct_sigma_list = []
      readsig = "off"
      #
      command = self.open_command(self.indmfl, openfile="f")
      exec command
      for line in f:
         line = line[:-1]
         line = line.split()
         if len(line) >= 4+1:
            if line[4] == "cix-num,":
               cix = int(line[0])
               cix_num_list.append(cix)
               dim = int(line[1])
               dim_list.append(dim)
               sigmat = np.zeros((dim,dim))
         if len(line) >= 2+1:
            if line[2] == "Sigind":
               readsig = "wait one more"
               count = 0
         if readsig == "wait one more": 
            readsig = "on"
         elif readsig == "on":
            ll = []
            for l in line:
               ll.append(int(l))
            ll = np.array(ll)
            sigmat[count,:] = ll
            count = count + 1
            if count == dim:
               readsig = "off"
               struct_sigma_list.append(sigmat)
      #
      return struct_sigma_list  # list all structure-sigma matrices

   def read_sigma_labels(self): # reads labels independent components sigma
      import numpy as np
      #
      labels_list = []
      readlabels = "off"
      #
      command = self.open_command(self.indmfl, openfile="f")
      exec command
      for line in f:
         line = line[:-1]
         line = line.split()
         if "Independent" in line:
            readlabels = "wait one more"
         if readlabels == "wait one more":
            readlabels = "on"
         elif readlabels == "on":
            labels = line
            # 
            for i,l in enumerate(labels):
               for char in l:
                  if char in "'":
                     labels[i] = labels[i].replace(char,"")
            #
            labels_list.append(labels)
            readlabels = "off"
      #
      return labels_list # independent components sigma (e.g., ["5/2","7/2"])

   def set_corr_atoms_init(self):
      # Assuming that equivalent atoms in self.struct are
      # treaded as equivalent also in GA
      import numpy as np
      iatom_list = self.read_label_corr()
      struct_sigma_list = self.read_struct_sigma()
      labels_list = self.read_sigma_labels()
      #
      atoms_corrstatus = []
      corr_atoms_instances = []
      ind = 0; indc=0
      for i,atom_name in enumerate(self.unequiv_atoms):
         c = "uncorrelated"
         atom = self.corr_atom_class(atom_name)
         # Instances only for the correlated unequivalent atoms
         if ind in iatom_list:
            c = "correlated"
            struct_sigma = np.copy(struct_sigma_list[indc])
            atom.set_struct_sigma(struct_sigma)
            atom.set_ind_sigma_labels(labels_list[indc])
            atom.set_order(i)
            atom.set_order2(ind)
            atom.set_mult(self.mult_unequiv_atoms[i])
            corr_atoms_instances.append(atom)
            indc=indc+atom.mult
         # 
         for rep in range(self.mult_unequiv_atoms[i]):
            ind = ind + 1
         #
         atoms_corrstatus.append(c)
         #
      self.nt_corr_atom=indc
      self.atoms_corrstatus = atoms_corrstatus
      self.corr_atoms_instances = corr_atoms_instances




#   def set_orblabels(self, orblabels_list):
#      # "s","p","d","f", or "uncorrelated"
#      assert(len(orblabels_list) == len(self.unequiv_atoms))
#      self.orblabels_list = orblabels_list
#      #
#      n_corr_atoms = np.copy(len(self.unequiv_atoms))
#      NA_list = []
#      Hilbert_list = []
#      for o in self.orblabels_list:
#         if o == "uncorrelated":
#            n_corr_atoms = n_corr_atoms - 1
#            if o == "s": NA = 2*(2*0+1)
#            if o == "p": NA = 2*(2*1+1)
#            if o == "d": NA = 2*(2*2+1)
#            if o == "f": NA = 2*(2*3+1)
#            NA_list.append(NA)
#            Hilbert_list.append(range(NA)) # by default it uses all of the local Hilbert space
#      self.n_corr_atoms = n_corr_atoms
#      self.NA_list = NA_list
#      self.Hilbert_list = Hilbert_list

#   def set_struct_sigma(self, struct_sigma_list):
#      # List of struct-sigma matrices (as in corr_atom_class below)
#      assert(len(struct_sigma_list) == self.n_corr_atoms)
#      for i,s in enumerate(struct_sigma_list):
#         assert(s.shape[0] == s.shape[1])
#         assert(s.shape[0] == NA_list[i])
#      self.struct_sigma_list = struct_sigma_list

   def get_INPUTFILES(self):
       GL_INP=''
       LSCF =self.get_LSCF()
       GL_INP=GL_INP+str(LSCF)+"       ! GL%LSCF: 1--updating VS; 4--generating VS; 6--fixed VS"+'\n'
       parametrization=self.get_parametrization()
       self.set_parametrization(parametrization)
       LGPRJ=self.get_LGPRJ()
       GL_INP=GL_INP+str(LGPRJ)+"      ! GL%LGPRJ: 0--FOCK PROJECTOR; 1--N,J,E; 3--N,S,E; 4--N,E,E,Sz"+'\n'
       LHUB =self.get_LHUB()
       GL_INP=GL_INP+str(LHUB) +"      ! GL%LHUB: 1--Slater; 2--Kanamori"+'\n'
       LDC=self.get_LDC()
       GL_INP=GL_INP+str(LDC)  +"      ! GL%LDC: 1--std dc; 2-fixed dc"+'\n'
       GL_INP=GL_INP+ "1           ! GL%LDIAT: 0--FULL T; 1--DIAGONAL T"+'\n'
       GL_INP=GL_INP+ "-10         ! GL%LEIGS: 0--FULL P; -N--RED P WITH ECUT=10^-N"+'\n'
       GL_INP=GL_INP+ "1000        ! GL%MAX_ITER"+'\n'
       GL_INP=GL_INP+ "0           ! GL%LPSICOMP: Calculate <Psi|LO>?"+'\n'
       GL_INP=GL_INP+ "0           ! GL%LS2"+'\n'
       GL_INP=GL_INP+ "0           ! GL%LL2"+'\n'
       GL_INP=GL_INP+ "0           ! GL%LJ2"+'\n'
       GL_INP=GL_INP+ "0           ! GL%LEBANDCOMP : Band energy component analysis"+'\n'
       GL_INP=GL_INP+ "0           ! GL%LENTANGLEMENT : Entanglement entropy"+'\n'

       corr_atom_type,type_1atom=self.get_corr_atom_type()
       n_corr_atom_type=len(corr_atom_type)
       GL_INP=GL_INP+str(n_corr_atom_type)+"  "+str(self.nt_corr_atom)+"  ! NTYP, NIONS"+'\n'
       GL_INP=GL_INP+"CORRELATED ATOM TYPE INFO"+'\n'
       for i,type in enumerate(corr_atom_type):
           GL_INP=GL_INP+"NT="+str(i+1)+'\n'
           prompt="Please enter U J parameters (eV) for correlated atom type %s" %(type)+'\n'
           userin = raw_input(prompt).split()
           U=userin[0]; J=userin[1]
           print "You have entered U=%s(eV) J=%s(eV)" %(U,J)+'\n'
           GL_INP=GL_INP+U+'  '+J+"  ! U,J"+'\n'
           NA2=self.corr_atoms_instances[type_1atom[i]].NA
           GL_INP=GL_INP+str(NA2/2)+"   ! NA"+'\n'
           prompt="Please enter Hilbert space occupation range Nf1 Nf2 (%2d =<Nf1<=Nf2<= %2d) " %(0,NA2)+'\n'
           userin = raw_input(prompt).split()
           Nf1=userin[0]; Nf2=userin[1]
           print "You have entered Nf1=%s Nf2=%s" %(Nf1,Nf2)+'\n'
           Nf_DC='  '
           if LDC==2:
                prompt="Please enter fixed Nf_DC for fixed DC (%s < Nf_DC < %s):"%(Nf1,Nf2)+'\n'
                Nf_DC  = raw_input(prompt).strip()
                print "You have enetered fixed Nf0=%s" %(Nf_DC)+'\n'
           GL_INP=GL_INP+Nf1+'  '+Nf2+'  '+Nf_DC+" ! OCC, Nf_DC"+'\n'
       GL_INP=GL_INP+"CORRELATED ATOM INFO"+'\n'
       icatm=0
       for i,atom in enumerate(self.corr_atoms_instances):
           NT=corr_atom_type.index(atom.nameL)+1
           for j in range(atom.mult):
               icatm=icatm+1
               GL_INP=GL_INP+"NI=%2d"%(icatm)+'\n'
               GL_INP=GL_INP+"%2d  %2d  %2d"%(NT,atom.order2+j+1,icatm-j)+" ! NT, NI0, NIMAP"+'\n'
               GL_INP=GL_INP+''.join("%2d  "%(i) for i in atom.struct_sigma.diagonal())+" ! SYMID"+'\n'
       with open(self.GL,'w') as f:
           print >>f, GL_INP

   def get_parametrization(self):
       prompt="Choose parametrization method for srcreened Coulomb integrals: "+'\n'+ \
              " 0 -- Fock, density-density type;  "+'\n'+\
              " 1 -- Slater, F0, F2, F4 ... with fixed ratio;  "+'\n'+\
              " 2 -- Kanamori, U, U', J.  "+'\n'
       userin = raw_input(prompt).strip()
       if userin=='0': parametrization='Fock'
       elif userin=='2': parametrization='Kanamori'
       else: parametrization='Slater'
       print "You have chosen %s type parametrization." %(parametrization)+'\n'
       return parametrization

   def get_corr_atom_type(self): # for instance, Ce1 and Ce2 belong to same atomic type Ce
       corr_atom_type=[]; type_1atom=[]
       for i,a in enumerate(self.corr_atoms_instances):
           if not a.nameL in corr_atom_type:
               corr_atom_type.append(a.nameL)
               type_1atom.append(i)
       return corr_atom_type,type_1atom

   def get_LSCF(self):
       LSCF=6
       prompt="Choose SCF type: "+'\n'+ \
              " 1 -- updating variational setup run;  "+'\n'+\
              " 4 -- generating initial variational setup;  "+'\n'+\
              " 6 -- fixed variational setup run.  "+'\n'
       userin = raw_input(prompt).strip()
       if userin=='1': LSCF=1
       elif userin=='4': LSCF=4
       print "You have chosen LSCF = %d" %(LSCF)+'\n'
       return LSCF

   def get_LGPRJ(self):
       if self.parametrization.lower()=='fock':
           LGPRJ=0
       else:
           prompt="Choose quantum labels for the diagonal Gutzwiller projector: "+'\n'+ \
              " 1 -- occupation N, total angular momentum J^2, atomic eigen energy E;  "+'\n'+\
              " 3 -- occupation N, total spin S^2, atomic eigen energy E;  "+'\n'+\
              " 4 -- occupation N, total spin S^2, atomic eigen energy E, S_z (spin symmetry broken).  "+'\n'
           userin = raw_input(prompt).strip()
           if userin=='3': LGPRJ=3
           elif userin=='4': LGPRJ=4
           else: LGPRJ=1
           print "You have chosen LGPRJ = %d" %(LGPRJ)+'\n'
       return LGPRJ

   def get_LHUB(self):
       if self.parametrization.lower()=='fock' or self.parametrization.lower()=='kanamori':
           LHUB=2
           print "You have chosen Kanamori type Coulomb interactions."+'\n'
       else: 
           LHUB=1
           print "You have chosen Slater type Coulomb interactions."+'\n'
       return LHUB

   def get_LDC(self):
       prompt="Choose double couting type: "+'\n'+ \
              " 1 -- standard double counting (LDA+U);  "+'\n'+\
              " 2 -- fixed double couting through fixed Nf to be provided.  "+'\n'
       userin = raw_input(prompt).strip()
       if userin=='2': LDC=2
       else: LDC=1
       print "You have chosen %d" %(LDC)+'\n'
       return LDC

   #----------------------------------------
   # METHODS TO READ RESULTS GA CALCULATION ONCE IT'S ALL DONE # %%%

   def set_energycell(self): # total energy cell (eV)
       import numpy as np
       import os
       #
       command = self.grep_command(":ENE", self.scf, grepfile="f")
       exec command
       #
       for line in f:
           line = line[:-1]
           line = line.split()
           energycell = float(line[8]) * self.convE  # Total energy given in Ry!!
       # 
       self.energycell = energycell

   def get_W(self,icatm):
       lread=0
       W=[]; n_block=[]
       with open(self.output,'r') as f:
           for line in f.readlines():
               if 'OCC_CONFIG_WEIGHT' in line:
                   if lread==1:
                       break
                   line=line.split()
                   NI=int(line[1])
                   if NI==icatm:
                       lread=1
               elif lread==1:
                   if 'OCC=' in line:
                       line=line.split()
                       n_block.append(int(line[1]))
                       W.append(float(line[3]))
                   else:
                       break
       return n_block,W

   def get_list_X(self,struct_sigma,icatm,STR):
      import numpy as np
      lread=0
      with open(self.output,'r') as f:
         for line in f.readlines():
            if STR in line:
               lread=1
            elif lread==1:
               line=line.split()
            if line[0]=='NI=' and int(line[1])==icatm:
               lread=2
            elif lread==2:
               line=line.split()
               list=[]
               for i in line:
                  list.append(float(i))
                  lread=0
         ind_sigma = np.amax(struct_sigma)
         list_X = []
         for i in range(ind_sigma):
            list_X.append(0.0)
         NA = struct_sigma.shape[0]
         for i in range(NA):
            ii=struct_sigma[i][i]
            list_X[ii-1]=list[i]

   #----------------------------------------
   # CLASS FOR THE CORRELATED ATOMS (used to analyze GA results) # %%%
   class corr_atom_class():

      def __init__(self, name):
         import numpy as np
         self.name  = name
         self.nameL = ''.join([i for i in name if not i.isdigit()])
       
      # INFORMATION AVAILABLE ALREADY AFTER INITIALIZATION WITH ga_init_dmft.py

      def set_order(self,order): # Wien2k label of groups of atoms in material
         self.order = order
      #
      def set_order2(self,order2): # Wien2k label of atoms in material
         self.order2=order2
      #
      def set_mult(self, mult): # multiplicity of correlated atom in material
         self.mult = mult
      #
      def set_positions(self, positions):
         assert(len(positions) == mult)
         self.positions = positions # (list of vectors, each one a 3d numpy array)

      def set_struct_sigma(self, struct_sigma): # 2d numpy array of zeros, ones, and so on
         # (assuming 0,1,2,3... convention, important that components are INTEGERS)
         import numpy as np
         assert(struct_sigma.shape[0] == struct_sigma.shape[1])
         self.struct_sigma = struct_sigma # (two-labels numpy array)
         self.NA = self.struct_sigma.shape[0] # size self-energy matrix
         maxlabel=int(np.amax(self.struct_sigma))
         ind_sigma_list = []
         for i in range(self.NA):
            for j in range(self.NA):
               s = int(self.struct_sigma[i,j])
               if s >= 1:
                  ind_sigma_list.append(s)
         isl = np.array(ind_sigma_list)
         self.ind_sigma = max(isl)-min(isl) + 1 # number of independent components self-energy 
         #
         dims = np.zeros(self.ind_sigma)
         for comp in range(self.ind_sigma):
            for i in range(self.NA):
               for j in range(self.NA):
                  if struct_sigma[i,j] == comp+1:
                     dims[comp]=dims[comp]+1
         self.dims = dims.tolist()
         #      
         if self.NA == 2*(2*0+1):
            self.orb_label = "s"
         if self.NA == 2*(2*1+1):
            self.orb_label = "p"
         if self.NA == 2*(2*2+1):
            self.orb_label = "d"
         if self.NA == 2*(2*3+1):
            self.orb_label = "f"

      def set_ind_sigma_labels(self, ind_sigma_labels): # labels of self-energy components
         assert(len(ind_sigma_labels) == self.ind_sigma)
         self.ind_sigma_labels = ind_sigma_labels # (list, e.g., ["5/2","7/2"])

      # (The next two methods are convenience routines to convert: lists <---> matrices)
      #
      def list_to_mat(self, l): # creating matrix M from list l
         import numpy as np
         M = np.zeros(shape=(self.NA,self.NA))
         for i in range(self.NA):
            for j in range(self.NA):
               if self.struct_sigma[i,j] >= 1:
                  M[i,j] = l[self.struct_sigma[i,j]-1]
         return M
      #
      def mat_to_list(self, M): # creating list from matrix M
         import numpy as np
         l = np.zeros(self.ind_sigma)
         for comp in range(1,self.ind_sigma+1):
            for i in range(self.NA):
               for j in range(self.NA):
                  if self.struct_sigma[i,j] == comp:
                     l[comp-1] = M[i,j]
         return l.tolist()
      
      # ----------------------------------------------------------------------
      # INFORMATION AVAILABLE ONLY AFTER THAT THE GA CALCULATION IS DONE
      # (See Eq. 30 draft DMFT-GA for the meaning of the quantities below)
                  
      def set_R(self, list_R): # renormalization-matrix non-zero values
         # (assuming convention of struct_sigma for the order!)
         import numpy as np
         assert(len(list_R) == self.ind_sigma)
         self.list_R = list_R # (list, e.g., [0.7, 0.3])
         #
         self.R = self.list_to_mat(list_R)
         self.Z = np.array( np.matrix(self.R) * np.matrix(self.R).H )
         self.list_Z = self.mat_to_list(self.Z)

      def set_Lambda(self, list_Lambda): 
         # (assuming convention of struct_sigma for the order!)
         import numpy as np
         assert(len(list_Lambda) == self.ind_sigma)
         self.list_Lambda = list_Lambda # (list, e.g., [0.7, 0.3])
         #
         self.Lambda = self.list_to_mat(self.list_Lambda)

      def set_Eta(self, list_Eta = "zero"):
         # (assuming convention of struct_sigma for the order!) 
         import numpy as np
         if list_Eta == "zero": list_Eta = np.zeros(self.ind_sigma).tolist()
         #
         assert(len(list_Eta) == self.ind_sigma)
         self.list_Eta = list_Eta # (list, e.g., [0.7, 0.3])
         #
         self.Eta = self.list_to_mat(self.list_Eta)
      
      def set_mu(self, EF):
         self.Mu = EF  # Fermi-level

      def set_onsite_energies(self, list_cf):
         # independent quadratic part of on-site Hamiltonian
         # (assuming convention of struct_sigma for the order!)
         self.list_cf = list_cf
         self.cf = self.list_to_mat(self.list_cf)

      def get_Sigma(self):
         import numpy as np
         from numpy.linalg import inv
         #
         identity = np.eye(self.NA)
         full_S0 = np.array( np.matrix(inv(self.R)) * \
                                np.matrix(self.Lambda + self.Eta) * \
                                np.matrix(inv(self.R)).H ) - \
                                self.Mu * np.array( np.matrix(identity - self.Z) * np.matrix(inv(self.Z)) )
         # Here Sigma is extrapolated according to the standard DMFT convention
         # (the on-site energies are not included in the definition of Sigma0
         self.Sigma0 = full_S0 - self.cf  # Sigma at the Fermi-level (note that I subtract the crystal field)
         self.Sigma1 = -np.array( np.matrix(identity-self.Z) * np.matrix(inv(self.Z)) ) # First order expansion

      # -----------------------------------------------------------------------

      def set_W(self, n_block, W): # information configuration weights
         self.n_block = n_block  # list not-neglected local eigenspaces of number-operator
         n_symbols = []
         for n in n_block:
            n_symbols.append("$" + self.orb_label + "^{" + str(n) + "}$")
         self.n_symbols = n_symbols
         self.W = W  # configuration weights ordered as n_block

      def set_n(self,list_n): # populations-matrix non-zero values
         import numpy as np
         assert(len(list_n) == self.ind_sigma)
         self.list_n = list_n # (list, e.g., [0.4,0.6])
         self.n = self.list_to_mat(list_n)  # creating matrix from list

      # ------------------------------------------------------------------------
 
#      def set_order(self,order): # Wien2k label of groups of atoms in material
#         self.order = order

      def set_Hloc(self, Hint_parametrization, U,J, list_HlocKS): # local Hamiltonian
         self.Hint_parametrization = Hint_parametrization # (e.g., "Slater")
         self.U,self.J = U,J
         assert(len(list_HlocKS) == self.ind_sigma)  # (list, e.g., [4.5,0.6])
         self.list_HlocKS = list_HlocKS
         self.HlocKS = self.list_to_mat(list_HlocKS)  # creating matrix from list

      def set_ladder(self, f,fdagger): # ladder operators stored as sparse matrices
         self.f = f  # (order consistent with self-energy matrix)
         self.fdagger = fdagger
   #    
   # ------------------------------------------------------------------------
   # ------------------------------------------------------------------------


   # To move above together with the methods to read the data
   # after that the calculation is done. 
   # Note that the instances are already created if the calculation
   # has been already performed, so the atoms instances already 
   # exist, and one should just set the various remainging things

   def set_corr_atoms_instances_results(self):      
      icatm=0 # correlated atom index including the equivalent ones
      for i,a in self.corr_atoms_instances:

         list_R = self.get_list_X(struct_sigma,icatm,'RAA')
         a.set_R(list_R)

         list_n = self.get_list_X(struct_sigma,icatm,'NKS')
         a.set_n(list_n)

         n_block,W = self.get_W(icatm)
         a.set_W(n_block, W)
 
         icatm=icatm+a.mult
   
   

   # --------- (convenience methods) ---------- #
   #
   def grep_command(self,string,filename, grepfile="f",option=""):
       import os
       command = grepfile + " = os.popen(\'grep \"" + option + string + "\" " + filename + "\')"
       return command
   #
   def open_command(self,filename, openfile="f"):
       import os
       command = openfile + "= open(\'" + filename + "\')"
       return command
   # -----------------



       

#######################################################
#######################################################



# TESTING #

import numpy as np

folder = "/home/ykent/WIEN_GUTZ/EXAMPLE/SmB6/SmB6"
#folder = "/home/ykent/WIEN_GUTZ/EXAMPLE/TESTSET/Ce_V28.5/Ce"
#folder = "./../new-Plutonium/ALPHA-Pu/GA-8tenth/alpha_Pu_3"
#folder = "./../new-Plutonium/ALPHA-Pu/LDA/alpha_Pu_3"
##folder = "./f_Pr_18"
##folder = "../../ALPHA_U-Pr/GA-6/a_Pr_18"

#folder = "./SmB6"

material = material_class(folder)
material.print_info()
material.get_INPUTFILES()



#material.analyze_GA() # Analyze GA calculation
#quit()
##material.set_name_material("XYZ")
#print "Name material:", material.name_material
#print
#print "Paths files:"
#print material.struct
#print material.outputkgen
#print material.indmfl
#print material.output
#print material.scf
#print material.GL
#print material.GLX

#print
#print "List groups unequivalent atoms:", material.unequiv_atoms
#print "Multiplicity groups unequivalent atoms:", material.mult_unequiv_atoms
#print "Total number of atoms:", material.n_atoms
#print "Volume/atom:", material.volume
#print "Total energy:", material.energy
#print

#corr_atoms = material.corr_atoms_instances
#print "CORRELATED ATOMS:"
#print
#for i in range(len(corr_atoms)):
#   a=corr_atoms[i]
#   print 'dims:', a.dims
#   print 'Atom name:', a.name, a.orb_label
#   print 'Renormalization matrix:'
#   print a.R
##   print a.struct_sigma
#   print 'Independent components Z matrix:', a.list_Z
##   print a.ind_sigma
#   print a.Z #np.diag(a.Z)
#   print 'ETA:'
#   a.set_Eta()
#   print a.Eta
##   print 'Orbital populations:'
##   print a.n
#   print 'Independent components orbital populations:', a.list_n
#   print "Configuration weights:"
#   for i,w in enumerate(a.W):
#      print a.n_block[i], a.n_symbols[i], w
   
#   print "Orbital labels:", a.ind_sigma_labels
#   print #a.mult

