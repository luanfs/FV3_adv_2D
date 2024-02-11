import numpy as np
import matplotlib.pyplot as plt

# Python script to plot the outputs

# directory
graphdir ='../graphs/' # must exist
datadir  ='../data/' # must exist
figformat = 'png'

# some constants
N    = 96 # number of cells
tc   = 4  # test case
hord = 0  # 1d adv scheme
adv   = 1  #2d adv scheme
nplots = 12


# x axis points for plotting
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
x, y = np.meshgrid(x,y)

if tc == 1 or tc == 2 or tc == 3:
   qmin = -0.1
   qmax =  1.2
elif tc == 4:
   qmin = -0.1
   qmax =  3.3

for t in range(0, nplots+1):
   # basename for plotting
   basename = "tc"+str(tc)+"_N"+str(N)+"_hord"+str(hord)+"_adv"+str(adv)+"_t"
   input_name  = datadir+basename+str(t)+'.txt'
   output_name = graphdir+basename+str(t)+'.'+figformat
   data = np.loadtxt(input_name)

   # plot the graph
   z = data[3:]
   z = np.reshape(z, (N,N))
   print(np.amin(z), np.amax(z))
   plt.contourf(x, y, z, cmap='jet', levels=np.linspace(qmin,qmax,20))
   plt.colorbar(orientation='vertical', fraction=0.046, pad=0.04, format='%.1e')

   time = data[0]
   time = str("{:.2e}".format(time))

   massvar = data[1]
   massvar = str("{:.2e}".format(massvar))

   cfl = data[2]
   cfl = str("{:.2e}".format(cfl))

   # Label
   plt.xlabel('$x$')
   plt.ylabel('$y$')
   plt.grid(True, which="both")
   title = "N="+str(N)+", time = "+time+" days, CFL="+cfl#+", mass variation="+massvar  
   plt.title(title)
   plt.savefig(output_name+'.'+figformat, format=figformat)
   plt.close()
   print("Plotted "+ output_name)
