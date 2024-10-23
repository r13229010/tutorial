import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta, datetime

# ---------------------------

# (k1, k2, j1, j2, i1, i2)
domain_range=(0, 50, 0, 128, 0, 128)
time = np.arange(0,721)
time_steps = 721
start_time = datetime(2024, 9, 18, 0, 0)  # Start at 00:00
time_datetime = [start_time + timedelta(minutes=2*i) for i in range(time_steps)]  # 2-minute intervals
lat=np.arange(0,25600,200)
lon=np.arange(0,25600,200)
exp_name = 'HW6a'

class VVMplotter:
    
    def plot_contourf_with_PBL(self, variable, title, time, z=np.arange(0,2000,40), PBL1=None, PBL2=None, cmap='seismic', cmin=None, cmax=None):
        plt.figure(figsize=(10, 6))
        plt.contourf(time, z, variable.T, levels=np.linspace(-0.06,0.06,25), cmap=cmap, extend='both')
        # plt.contourf(time, z, variable.T, cmap=cmap, extend='both')
        plt.colorbar()
        plt.contour(time, z, variable.T, linewidths = 1, levels=[0], colors = 'black')
        plt.plot(time, PBL1, color='red', label='PBL (th = surface_th + 0.5K)', linewidth=2)
        plt.plot(time, PBL2, color='lightgreen', label='PBL (dth/dz max)', linewidth=2)
        plt.title(title)
        plt.xlabel('Time')
        plt.ylabel('Height (m)')
        
        plt.legend(loc='upper right')
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=1))
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(exp_name + '_' +title + '.png', dpi=300)
        # plt.show()
        plt.close()
    
    def plot_horizontal(self, tracer, title, hours=np.arange(0,12,2), z=160, cmap='seismic', cmin=-0.06, cmax=0.06):
        i=0
        for t in hours*30:
            time_str = str(int(hours[i])).zfill(2)
            lev = int(z/40)
            
        
            fig, ax = plt.subplots(figsize=(7, 6))
            plotter = ax.contourf(lon,lat,tracer[t,lev],levels=np.linspace(cmin,cmax,25),cmap=cmap,extend='both')
            # plotter = ax.contourf(lon,lat,tracer[num,0],cmap=cmap,extend='both')
            ax.contour(lon, lat, tracer[t,lev], linewidths = 1, levels=[0], colors = 'black')
            cbar = fig.colorbar(plotter)
            # cbar.set_label("[ppb]")
            ax.set_xlabel('x [m]')
            ax.set_ylabel('y [m]')
            ax.set_aspect('equal')
            ax.set_title(title + str(z) + 'm ' + ' at ' + time_str + ':00')
            plt.tight_layout()
            plt.savefig(exp_name + '_'+str(z) + 'm_'+time_str + '.png',dpi=300)
            plt.close()
            i=i+1
            
    def plot_profile(self, u,w,wth,title, hours=np.arange(0,12,2), z=np.arange(0,2000,40), cmap='seismic', cmin=None, cmax=None):
        space = 2
        for t in hours*30:
            time_str = str(int(t/30)).zfill(2)
            wth_plot = np.nanmean(wth[t],axis=1)
            fig, ax = plt.subplots(figsize=(12, 5))
            plotter = ax.contourf(lon,z,wth_plot,levels=np.linspace(cmin,cmax,25),cmap=cmap,extend='both')
            wind_plot = ax.quiver(lon[::space],z[::space],np.nanmean(u[t],axis=1)[::space,::space],np.nanmean(w[t],axis=1)[::space,::space],color='black',scale = 10)
            ax.contour(lon, z, wth_plot, linewidths = 1, levels=[0], colors = 'black')
            fig.colorbar(plotter)
            # cbar.set_label("[ppb]")
            ax.set_xlabel('x [m]')
            ax.set_ylim(0,1000)
            ax.set_ylabel('z [m]')
            # ax.set_aspect('equal')
            ax.set_title(title + ' at ' + time_str + ':00')
            plt.tight_layout()
            plt.savefig(exp_name + ' wind profile_' + time_str + '.png', dpi=300)
            plt.close()


wth = np.load(exp_name + '_' + 'wth.npy')
PBL1_1, PBL2_1=np.load(exp_name + '_' +'PBL 1.npy')
PBL1_2, PBL2_2=np.load(exp_name + '_' +'PBL 2.npy')
PBL1_12, PBL2_12=np.load(exp_name + '_' +'PBL 1+2.npy')
Turbulance_1 = np.load(exp_name + '_' +'Turbulance 1.npy')
Turbulance_2 = np.load(exp_name + '_' +'Turbulance 2.npy')
Turbulance_12 = np.load(exp_name + '_' +'Turbulance 1+2.npy')
u = np.load(exp_name + '_' +'u.npy')
w = np.load(exp_name + '_' +'w.npy')
plotter = VVMplotter()
plotter.plot_profile(u=u, w=w, wth=wth, title='Sounding 1 ', cmin=-0.06, cmax=0.06)
plotter.plot_horizontal(wth, title='Sounding 1', hours=np.arange(0,12,2), z=160, cmap='jet', cmin=-0.01, cmax=0.06)

plotter.plot_contourf_with_PBL(np.array(Turbulance_1),'Sounding 1 pasture',time = time_datetime, PBL1=PBL1_1, PBL2=PBL2_1)

plotter.plot_contourf_with_PBL(np.array(Turbulance_2),'Sounding 1 evergreen needle',time = time_datetime, PBL1=PBL1_2, PBL2=PBL2_2)

plotter.plot_contourf_with_PBL(np.array(Turbulance_12),'Sounding 1 whole domain',time = time_datetime, PBL1=PBL1_12, PBL2=PBL2_12)

