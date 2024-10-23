import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta, datetime
# class Person:
#     def __init__(self, Age, name):
#         self.age = Age # other function can use 'age'
#         pass
#     def greet(self):
#         print(self.age)
# Person(Age=20, name='BeiBei')
# # 繼承
# class Student(Person):
#     def __init__(self, Age, Name, ID):
#         super().__init__(Age=Age, name=Name)
#         self.ID=ID
# mog = Student(Age=20, Name='mog', ID='R13229069')

# print(mog.ID)


# ---------------------------
from vvmtools import VVMTools

def wthbar(w,th):
    
    # calculate x-y mean (axis=0)
    w_mean = np.mean(w, axis=(2,3))        # shape: (t, z)
    th_mean = np.mean(th, axis=(2,3))  # shape: (t, z)
    # 計算擾動項
    w_prime = w - w_mean[:, :,np.newaxis, np.newaxis]  # shape: (time, z, y, x)
    th_prime = th - th_mean[:, :,np.newaxis, np.newaxis]  # shape: (time, z, y, x)
    
    # 計算擾動項乘積的平均 (時間方向)
    # print((w_prime * th_prime).shape)
    w_th_prime_mean = np.mean(w_prime * th_prime, axis=(2,3))  # shape: (y, x)
    # print('mean = ',(w_th_prime_mean))
    return w_th_prime_mean
      # 應該是 (y, x)
def wth(w,th):
    # calculate x-y mean (axis=0)
    w_mean = np.mean(w, axis=(2,3))        # shape: (t, z)
    th_mean = np.mean(th, axis=(2,3))  # shape: (t, z)
    # 計算擾動項
    w_prime = w - w_mean[:, :,np.newaxis, np.newaxis]  # shape: (time, z, y, x)
    th_prime = th - th_mean[:, :,np.newaxis, np.newaxis]  # shape: (time, z, y, x)
    return w_prime * th_prime
    
# Calculate PBL Heights
def calc_PBL_heights(th, z):
    PBL1 = np.zeros(th.shape[0])  # Definition 1
    PBL2 = np.zeros(th.shape[0])  # Definition 2
    
    for t in range(th.shape[0]):
        # PBL Definition 1: th = surface_th + 0.5K, exclude the lowest layer
        surface_th = th[t, 0]
        target_th = surface_th + 0.5
        # Only search above the first layer (index 1 and above)
        PBL1[t] = z[np.argmin(np.abs(th[t, :] - target_th))]  # Find height where th is closest to target_th, excluding the lowest level
        
        # PBL Definition 2: dth/dz max (second derivative of th), exclude the lowest layer
        dth_dz = np.gradient(th[t,:], z)
        # Only search above the first layer (index 1 and above)
        PBL2[t] = z[np.argmax(dth_dz)]  # Find height where second derivative is maximum, excluding the lowest level
    
    return PBL1, PBL2

class newVVMTools(VVMTools):
    def __init__(self, case_path):
        super().__init__(case_path)
        
    def TKE(self, t):
        u = self.get_var('u', t, numpy=True)
        v = self.get_var('v', t, numpy=True)
        w = self.get_var('w', t, numpy=True)
        
        return np.nanmean(1/2*(u**2 + v**2 + w**2),axis = (1,2))
    
    def Enstrophy(self, t):
        
        zeta = self.get_var('zeta', t, numpy=True)
        eta = self.get_var('eta_2', t, numpy=True)
        if zeta.shape != eta:
            eta = self.get_var('eta', t, numpy=True)
        xi = self.get_var('xi', t, numpy=True)
        return np.nanmean((eta**2 + zeta**2 + xi**2),axis = (1,2))

    def Turbulance(self, time):
        w_list = []
        th_list = []
    
        # 遍歷每個時間步
        for t in time:
            w = self.get_var('w', t, numpy=True)
            th = self.get_var('th', t, numpy=True)
            w_list.append(w)
            th_list.append(th)
    
        # 將列表轉換為 numpy 陣列，shape: (time, z, y, x)
        w_all = np.stack(w_list, axis=0)
        th_all = np.stack(th_list, axis=0)
        
        return w_all, th_all
    
    def UW(self, t):
        w_list = []
        u_list = []
    
        # 遍歷每個時間步
        for t in time:
            w = self.get_var('w', t, numpy=True)
            u = self.get_var('u', t, numpy=True)
            w_list.append(w)
            u_list.append(u)
    
        # 將列表轉換為 numpy 陣列，shape: (time, z, y, x)
        w_all = np.stack(w_list, axis=0)
        u_all = np.stack(u_list, axis=0)
        return u_all,w_all

def plot_contourf_with_PBL(variable, title, time, z=np.arange(0,2000,40), PBL1=None, PBL2=None, cmap='jet', cmin=None, cmax=None):
    plt.figure(figsize=(10, 6))
    # plt.contourf(time, z, variable.T, levels=np.linspace(cmin,cmax,100), cmap=cmap, extend='both')
    plt.contourf(time, z, variable.T, cmap=cmap, extend='both')
    plt.colorbar()
    plt.contour(time, z, variable.T, linewidths = 1, levels=[0], colors = 'black')
    # plt.plot(time, PBL1, color='red', label='PBL (th = surface_th + 0.5K)', linewidth=2)
    # plt.plot(time, PBL2, color='blue', label='PBL (dth/dz max)', linewidth=2)
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel('Height (m)')
    
    # plt.legend(loc='upper right')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=1))
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(exp_name + '_' +title + '.png', dpi=300)
    plt.show()
    # plt.close()


# -----------------------------------------------------------------------------
#                                main                                         |
# -----------------------------------------------------------------------------
exp_name = 'HW6a'
VVMObjects = VVMTools('/data/mlcloud/r13229010/VVM/DATA/case_HW5a')
newtools = newVVMTools('/data/mlcloud/r13229010/VVM/DATA/case_HW5a')
# (k1, k2, j1, j2, i1, i2)
domain_range=(0, 50, 0, 128, 0, 128)
z = np.arange(0,2000,40)
time = np.arange(0,721)
time_steps = 721
start_time = datetime(2024, 9, 18, 0, 0)  # Start at 00:00
time_datetime = [start_time + timedelta(minutes=2*i) for i in range(time_steps)]  # 2-minute intervals



u,w = newtools.UW(time)
np.save(exp_name + '_' +'u.npy', u)
np.save(exp_name + '_' +'w.npy', w)

# w_all ,th_all = newtools.Turbulance(time = time)
# PBL1, PBL2 = calc_PBL_heights(np.nanmean(th_all, axis=(2,3)), z)

# np.save(exp_name + '_' +'wth.npy', wth(w_all, th_all))
# np.save(exp_name + '_' +'PBL 1.npy',calc_PBL_heights(np.nanmean(th_all[:,:,:,:64], axis=(2,3)), z))
# np.save(exp_name + '_' +'PBL 2.npy',calc_PBL_heights(np.nanmean(th_all[:,:,:,64:], axis=(2,3)), z))
# np.save(exp_name + '_' +'PBL 1+2.npy',calc_PBL_heights(np.nanmean(th_all, axis=(2,3)), z))

# Turbulance_1 = wthbar(w_all[:,:,:,:64],th_all[:,:,:,:64])
# Turbulance_2 = wthbar(w_all[:,:,:,64:],th_all[:,:,:,64:])
# Turbulance_12 = wthbar(w_all,th_all)
# np.save(exp_name + '_' +'Turbulance 1.npy', Turbulance_1)
# np.save(exp_name + '_' +'Turbulance 2.npy', Turbulance_2)
# np.save(exp_name + '_' +'Turbulance 1+2.npy', Turbulance_12)
print(exp_name + ' done!')