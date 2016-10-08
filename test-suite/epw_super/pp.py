#
# Post-processing script QE --> EPW 
# 14/07/2015 - Samuel Ponce
#
import os

# Enter the number of irr. q-points 
user_input = raw_input('Enter the prefix used for PH calculations (e.g. diam)\n')
prefix = str(user_input)

# Enter the number of irr. q-points 
user_input = raw_input('Enter the number of irreducible q-points\n')
nqpt = user_input
try:
  nqpt = int(user_input)
except ValueError:
  raise Exception('The value you enter is not an integer!')

os.system('mkdir save')

for iqpt in range(nqpt+1):
  label = str(iqpt)

  os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
  if (iqpt == 1):
    os.system('cp _ph0/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
    os.system('cp -r _ph0/'+prefix+'.phsave save/')
  else:
    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
    os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )

