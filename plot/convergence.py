import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path
import subprocess

# Python script to plot the outputs
# directory
graphdir ='../graphs/'
datadir  ='../data/'
pardir   ='../par/'
figformat = 'png'

# Program to be run
program = "./main"

# test case
tc = 3

# advection scheme
hords = (8, 8, 0, 0)

# departure point scheme
dps = (1, 2, 1, 2)

# N
N = 48
Ns=[]

# time step for N
if tc==1 or tc==2 or tc==3:
  dt = 14400
elif tc==4:
  dt = 7200
else:
  print('invalid TC')
  exit()

# number of grids
ngrids = 3

# aux routine
def replace_line(filename, content, line_number):
    import re
    if os.path.exists(filename): # The file exists
        # Open the grid file
        file  = open(filename, "r")
        lines = file.readlines()

        # new content
        lines[line_number-1] = content+'\n'

        # Close the file
        file.close()

        # Write new file
        with open(filename, 'w') as file:
            for line in lines:
                file.write(line)

    else:   # The file does not exist
        print("ERROR in edit_file_line: file "+filename+" not found in /par.")
        exit()
 
# error arrays
errors_linf = np.zeros((ngrids, len(dps)))
errors_l1   = np.zeros((ngrids, len(dps)))
errors_l2   = np.zeros((ngrids, len(dps)))
# compile the code
subprocess.run('cd .. ; make', shell=True)
parfile = pardir+'input.par'
for n in range(0, ngrids):
   for k in range(0, len(dps)):
      dp   = dps[k]
      hord = hords[k]
      # update parameters
      replace_line(parfile, str(tc)  , 3)
      replace_line(parfile, str(N)   , 5)
      replace_line(parfile, str(dt)  , 7)
      replace_line(parfile, str(hord), 9)
      replace_line(parfile, str(dp)  , 11)

      # Run the program
      subprocess.run('cd .. ; ./main ', shell=True)

      # error filename
      filename = "tc"+str(tc)+"_N"+str(N)+"_hord"+str(hord)+"_dp"+str(dp)+"_errors.txt"
   
      # load the errors
      errors = np.loadtxt(datadir+filename)
      errors_linf[n,k] = errors[0]
      errors_l1[n,k]   = errors[1]
      errors_l2[n,k]   = errors[2]

   # update N and dt
   Ns.append(N)
   N = 2*N
   dt = dt*0.5
   
 
# Plotting parameters
colors = ('green', 'green', 'red', 'red')
markers = ('*','+','x','*', '+', 'x', '*', '+')
lines_style = ('-','--','-','--')

# Plot error graph 
errors = [errors_linf, errors_l1, errors_l2]
names = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']
enames = ['linf','l1','l2']
for e in range(0, len(errors)):
   error = errors[e]
   for m in range(0, len(dps)):
      hord = hords[m]
      dp = dps[m]

      # convergence rate
      n = len(Ns)-1
      CR = (np.log(error[n-1,m])-np.log(error[n,m]))/np.log(2.0)
      CR = str("{:2.1f}".format(CR))

      # plot
      plt.loglog(Ns, error[:,m], lines_style[m], color=colors[m], marker=markers[m],\
      label = 'hord'+str(hord)+'-dp'+str(dp)+" - order "+CR)

   # Label
   title =names[e]+" error - TC"+str(tc)
   figname =  graphdir+'tc'+str(tc)+'_'+enames[e]+"_error"
   plt.xlabel('$N$')
   plt.ylabel('Error')
   plt.title(title)
   plt.legend()
   plt.grid(True, which="both")
   plt.savefig(figname+'.'+figformat, format=figformat)
   #plt.show()
   plt.close()
