#compare the free run between hete and homo

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from matplotlib import cm
from netCDF4 import Dataset

# Read coordinates
coordinate_file='/p/home/jusers/tang1/jureca/dev_hgs-pdaf/postprocessing/write_coordinate/coordinates.dat'
coord_data=np.loadtxt(coordinate_file)

# model simulation
homo_ol_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_homo_ol_new/001/hgs-flow.da.nc'
hete_ol_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_hete_ol_new/001/hgs-flow.da.nc'

homo_h_h_file='/p/home/jusers/tang1/jureca/dev_hgs-pdaf/run_scripts/jureca/work/n100_homo_h_005/001/hgs-flow.da.nc'
hete_h_h_file='/p/home/jusers/tang1/jureca/dev_hgs-pdaf/run_scripts/jureca/work/n100_hete_h_005/001/hgs-flow.da.nc'

homo_hk_h_file='/p/home/jusers/tang1/jureca/dev_hgs-pdaf/run_scripts/jureca/work/n100_homo_hk_005/001/hgs-flow.da.nc'
hete_hk_h_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_hete_hk_h_damp0.1/001/hgs-flow.da.nc'

homo_hs_h_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_homo_hs_h/001/hgs-flow.da.nc'
hete_hs_h_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_hete_hs_h/001/hgs-flow.da.nc'

homo_hs_s_file='/p/scratch/icei-prace-2023-0004/tang1/output/hgs-pdaf/n100_homo_hs_s/001/hgs-flow.da.nc'
hete_hs_s_file='/p/scratch/icei-prace-2023-0004/tang1/output/hgs-pdaf/n100_hete_hs_s/001/hgs-flow.da.nc'

homo_hs_hs_file='/p/scratch/icei-prace-2023-0004/tang1/output/hgs-pdaf/n100_homo_hs_hs/001/hgs-flow.da.nc'
hete_hs_hs_file='/p/scratch/icei-prace-2023-0004/tang1/output/hgs-pdaf/n100_hete_hs_hs/001/hgs-flow.da.nc'

homo_hsk_h_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_homo_hsk_h/001/hgs-flow.da.nc'
hete_hsk_h_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_hete_hsk_h_damp0.02/001/hgs-flow.da.nc'

homo_hsk_s_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_homo_hsk_s/001/hgs-flow.da.nc'
hete_hsk_s_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_hete_hsk_s/001/hgs-flow.da.nc'

homo_hsk_hs_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_homo_hsk_hs/001/hgs-flow.da.nc'
hete_hsk_hs_file='/p/home/jusers/tang1/jusuf/dev_hgs-pdaf/run_scripts/jusuf/work/n100_hete_hsk_hs_damp0.02/001/hgs-flow.da.nc'


homo_ol_data=Dataset(homo_ol_file)
hete_ol_data=Dataset(hete_ol_file)

homo_h_h_data=Dataset(homo_h_h_file)
hete_h_h_data=Dataset(hete_h_h_file)

homo_hk_h_data=Dataset(homo_hk_h_file)
hete_hk_h_data=Dataset(hete_hk_h_file)

homo_hs_h_data=Dataset(homo_hs_h_file)
hete_hs_h_data=Dataset(hete_hs_h_file)

homo_hs_s_data=Dataset(homo_hs_s_file)
hete_hs_s_data=Dataset(hete_hs_s_file)

homo_hs_hs_data=Dataset(homo_hs_hs_file)
hete_hs_hs_data=Dataset(hete_hs_hs_file)

homo_hsk_h_data=Dataset(homo_hsk_h_file)
hete_hsk_h_data=Dataset(hete_hsk_h_file)

homo_hsk_s_data=Dataset(homo_hsk_s_file)
hete_hsk_s_data=Dataset(hete_hsk_s_file)

homo_hsk_hs_data=Dataset(homo_hsk_hs_file)
hete_hsk_hs_data=Dataset(hete_hsk_hs_file)

# observations
obs_file='/p/project/icei-prace-2023-0004/tang1/input_HGS-PDAF/observation/obs_HEAD.nc'
obs_data=Dataset(obs_file)

# select observation point
# read observation file
# read the number of observations
nobs=obs_data.dimensions['n_obs'].size
x=obs_data.variables['x'][:]
y=obs_data.variables['y'][:]
z=obs_data.variables['z'][:]
id_obs=obs_data.variables['obs_id'][:]
nsteps=95

diff_homo_ol=np.zeros((nsteps,nobs))
diff_hete_ol=np.zeros((nsteps,nobs))
diff_homo_h_h=np.zeros((nsteps,nobs))
diff_hete_h_h=np.zeros((nsteps,nobs))
diff_homo_hk_h=np.zeros((nsteps,nobs))
diff_hete_hk_h=np.zeros((nsteps,nobs))
diff_homo_hs_h=np.zeros((nsteps,nobs))
diff_hete_hs_h=np.zeros((nsteps,nobs))
diff_homo_hs_s=np.zeros((nsteps,nobs))
diff_hete_hs_s=np.zeros((nsteps,nobs))
diff_homo_hs_hs=np.zeros((nsteps,nobs))
diff_hete_hs_hs=np.zeros((nsteps,nobs))
diff_homo_hsk_h=np.zeros((nsteps,nobs))
diff_hete_hsk_h=np.zeros((nsteps,nobs))
diff_homo_hsk_s=np.zeros((nsteps,nobs))
diff_hete_hsk_s=np.zeros((nsteps,nobs))
diff_homo_hsk_hs=np.zeros((nsteps,nobs))
diff_hete_hsk_hs=np.zeros((nsteps,nobs))

for j in range(0,nobs):
  title='avg_abs(forecast-obs)'
  obs_head=obs_data.variables['Head'][0:nsteps,j]
  head_homo_ol=homo_ol_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_ol=hete_ol_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_h_h=homo_h_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_h_h=hete_h_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_hk_h=homo_hk_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_hk_h=hete_hk_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_hs_h=homo_hs_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_hs_h=hete_hs_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_hs_s=homo_hs_s_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_hs_s=hete_hs_s_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_hs_hs=homo_hs_hs_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_hs_hs=hete_hs_hs_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_hsk_h=homo_hsk_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_hsk_h=hete_hsk_h_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_hsk_s=homo_hsk_s_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_hsk_s=hete_hsk_s_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_homo_hsk_hs=homo_hsk_hs_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  head_hete_hsk_hs=hete_hsk_hs_data.variables['head_f'][0:nsteps,id_obs[j]-1]
  diff_homo_ol[:,j]=np.abs(head_homo_ol-obs_head)
  diff_hete_ol[:,j]=np.abs(head_hete_ol-obs_head)
  diff_homo_h_h[:,j]=np.abs(head_homo_h_h-obs_head)
  diff_hete_h_h[:,j]=np.abs(head_hete_h_h-obs_head)
  diff_homo_hk_h[:,j]=np.abs(head_homo_hk_h-obs_head)
  diff_hete_hk_h[:,j]=np.abs(head_hete_hk_h-obs_head)
  diff_homo_hs_h[:,j]=np.abs(head_homo_hs_h-obs_head)
  diff_hete_hs_h[:,j]=np.abs(head_hete_hs_h-obs_head)
  diff_homo_hs_s[:,j]=np.abs(head_homo_hs_s-obs_head)
  diff_hete_hs_s[:,j]=np.abs(head_hete_hs_s-obs_head)
  diff_homo_hs_hs[:,j]=np.abs(head_homo_hs_hs-obs_head)
  diff_hete_hs_hs[:,j]=np.abs(head_hete_hs_hs-obs_head)
  diff_homo_hsk_h[:,j]=np.abs(head_homo_hsk_h-obs_head)
  diff_hete_hsk_h[:,j]=np.abs(head_hete_hsk_h-obs_head)
  diff_homo_hsk_s[:,j]=np.abs(head_homo_hsk_s-obs_head)
  diff_hete_hsk_s[:,j]=np.abs(head_hete_hsk_s-obs_head)
  diff_homo_hsk_hs[:,j]=np.abs(head_homo_hsk_hs-obs_head)
  diff_hete_hsk_hs[:,j]=np.abs(head_hete_hsk_hs-obs_head)

#yverage over nobs
avg_diff_homo_ol=np.average(diff_homo_ol,1)
avg_diff_hete_ol=np.average(diff_hete_ol,1)
avg_diff_homo_h_h=np.average(diff_homo_h_h,1)
avg_diff_hete_h_h=np.average(diff_hete_h_h,1)
avg_diff_homo_hk_h=np.average(diff_homo_hk_h,1)
avg_diff_hete_hk_h=np.average(diff_hete_hk_h,1)
avg_diff_homo_hs_h=np.average(diff_homo_hs_h,1)
avg_diff_hete_hs_h=np.average(diff_hete_hs_h,1)
avg_diff_homo_hs_s=np.average(diff_homo_hs_s,1)
avg_diff_hete_hs_s=np.average(diff_hete_hs_s,1)
avg_diff_homo_hs_hs=np.average(diff_homo_hs_hs,1)
avg_diff_hete_hs_hs=np.average(diff_hete_hs_hs,1)
avg_diff_homo_hsk_h=np.average(diff_homo_hsk_h,1)
avg_diff_hete_hsk_h=np.average(diff_hete_hsk_h,1)
avg_diff_homo_hsk_s=np.average(diff_homo_hsk_s,1)
avg_diff_hete_hsk_s=np.average(diff_hete_hsk_s,1)
avg_diff_homo_hsk_hs=np.average(diff_homo_hsk_hs,1)
avg_diff_hete_hsk_hs=np.average(diff_hete_hsk_hs,1)

t=np.arange(0,nsteps,1)
# we divide these scenarios into two groups:
# group 1: heterogeneous
# group 2: homogeneous

fig,[ax1, ax2]=plt.subplots(1,2,sharey=True)
ax2.plot(t,avg_diff_homo_ol,label='ol')
ax1.plot(t,avg_diff_hete_ol,label='ol')
ax2.plot(t,avg_diff_homo_h_h,label='DA_h_h')
ax1.plot(t,avg_diff_hete_h_h,label='DA_h_h')
ax2.plot(t,avg_diff_homo_hk_h,label='DA_hk_h')
ax1.plot(t,avg_diff_hete_hk_h,label='DA_hk_h')
ax2.plot(t,avg_diff_homo_hs_h,label='DA_hs_h')
ax1.plot(t,avg_diff_hete_hs_h,label='DA_hs_h')
ax2.plot(t,avg_diff_homo_hs_s,label='DA_hs_s')
ax1.plot(t,avg_diff_hete_hs_s,label='DA_hs_s')
ax2.plot(t,avg_diff_homo_hs_hs,label='DA_hs_hs')
ax1.plot(t,avg_diff_hete_hs_hs,label='DA_hs_hs')
ax2.plot(t,avg_diff_homo_hsk_h,label='DA_hsk_h')
ax1.plot(t,avg_diff_hete_hsk_h,label='DA_hsk_h')
ax2.plot(t,avg_diff_homo_hsk_s,label='DA_hsk_s')
ax1.plot(t,avg_diff_hete_hsk_s,label='DA_hsk_s')
ax2.plot(t,avg_diff_homo_hsk_hs,label='DA_hsk_hs')
ax1.plot(t,avg_diff_hete_hsk_hs,label='DA_hsk_hs')

ax1.set_title('(a) heterogeneous K')
ax2.set_title('(b) homogeneous K')

ax1.set_xlabel('time (day)')
ax2.set_xlabel('time (day)')
ax1.set_ylabel('average absolute head difference (m)')
ax1.legend()
plt.tight_layout()
#plt.show()
plt.savefig("abs_head_diff",dpi=800,bbox_inches='tight')
