import numpy as np
import matplotlib.pyplot as plt

# Python script to plot the outputs

# directory
graphdir ='../graphs/' # must exist
datadir  ='../data/' # must exist
figformat = 'png'

# some constants
N    = 192 # number of cells
tc   = 3  # test case
hord = 8 # 1d adv scheme
dp   = 1  #2d adv scheme
iadv = 1
nplots = 12


# Domain size
erad = 6371.0 # earth radius (km)
Lx = 0.5*np.pi*erad
Lxo2 = Lx*0.5

# x axis points for plotting
x = np.linspace(-Lxo2, Lxo2, N)
y = np.linspace(-Lxo2, Lxo2, N)
x, y = np.meshgrid(x,y)

if tc == 1 or tc == 2 or tc == 3:
   qmin =  0.0
   qmax =  1.2
elif tc == 4:
   qmin =  0.0
   qmax =  3.5

for t in range(0, nplots+1):
   # basename for plotting
   basename = "tc"+str(tc)+"_N"+str(N)+"_hord"+str(hord)+"_iadv"+str(iadv)+"_dp"+str(dp)+"_t"
   input_name  = datadir+basename+str(t)+'.txt'
   output_name = graphdir+'adv2d_'+basename+str(t)+'.'+figformat
   data_info = np.loadtxt(input_name)

   # get binary data
   input_name_bin  = datadir+basename+str(t)+'.dat'
   z = open(input_name_bin, 'rb')
   z = np.fromfile(z, dtype=np.float64)
   z = np.reshape(z, (N,N))

   # plot the graph
   print(np.amin(z), np.amax(z))
   plt.contourf(x, y, z, cmap='jet', levels=np.linspace(qmin,qmax,20))
   plt.colorbar(orientation='vertical', fraction=0.046, pad=0.04, format='%.1e')

   time = data_info[0]
   time = str("{:.2e}".format(time))

   massvar = data_info[1]
   massvar = str("{:.2e}".format(massvar))

   cfl = data_info[2]
   cfl = str("{:.2e}".format(cfl))

   if iadv==1:
      sp = 'PL'
   elif iadv==2:
      sp = 'LT'
   # Label
   plt.xlabel('$x$ (km)')
   plt.ylabel('$y$ (km)')
   #plt.xlim(-2000, 2000)
   #plt.ylim(-2000, 2000)
   plt.grid(True, which="both")
   title = "N="+str(N)+", time = "+time+" days, CFL="+cfl+ '\n'+ \
   sp +'.hord'+str(hord)+'.dp'+str(dp)
   #+", mass variation="+massvar  
   plt.title(title)
   plt.savefig(output_name, format=figformat)
   plt.close()
   print("Plotted "+ output_name)
