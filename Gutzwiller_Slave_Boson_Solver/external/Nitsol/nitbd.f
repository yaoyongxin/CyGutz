      block data nitbd

c  Purpose- Initialize nitparam common block with default values.

      include 'nitdflts.h'
      include 'nitparam.h'
      include 'nitprint.h'
  
      data choice1_exp  /DFLT_CHOICE1_EXP/
      data choice2_exp  /DFLT_CHOICE2_EXP/
      data choice2_coef /DFLT_CHOICE2_COEF/
      data eta_cutoff   /DFLT_ETA_CUTOFF/
      data etamax       /DFLT_ETA_MAX/
      data thmin        /DFLT_THMIN/
      data thmax        /DFLT_THMAX/
      data etafixed     /DFLT_ETA_FIXED/

      data iplvl        /DFLT_PRLVL/
      data ipunit       /STDOUT/

      end
