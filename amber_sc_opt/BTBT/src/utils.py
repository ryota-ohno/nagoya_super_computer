import numpy as np

def get_E(file):
    i=0;E_list=[]
    with open(file)as f:
        for line in f:
            if i==1:
                s=line.split()
                E=float(s[1])
                E_list.append(E)
                break
            if (line.find(' RMS ')>-1):
                i=1
    return E_list

def Rod(n,theta_in):
    nx,ny,nz=n
    theta_t=np.radians(theta_in)
    Rod=np.array([[np.cos(theta_t)+(nx**2)*(1-np.cos(theta_t)),nx*ny*(1-np.cos(theta_t))-nz*np.sin(theta_t),nx*nz*(1-np.cos(theta_t))+ny*np.sin(theta_t)],
                [nx*ny*(1-np.cos(theta_t))+nz*np.sin(theta_t),np.cos(theta_t)+(ny**2)*(1-np.cos(theta_t)),ny*nz*(1-np.cos(theta_t))-nx*np.sin(theta_t)],
                [nx*nz*(1-np.cos(theta_t))-ny*np.sin(theta_t),ny*nz*(1-np.cos(theta_t))+nx*np.sin(theta_t),np.cos(theta_t)+(nz**2)*(1-np.cos(theta_t))]])
    return Rod

def R2atom(R):
    if R==1.8:
        return 'S'
    elif R==1.7:
        return 'C'
    elif R==1.2:
        return 'H'
    else:
        return 'X'

def check_calc_status(df_cur,A1,A2,A3,a,b):
    try:        
        return df_cur.loc[
                        (df_cur['A1']==A1)&
                        (df_cur['A2']==A2)&
                        (df_cur['A3']==A3)&
                        (df_cur['a']==a)&
                        (df_cur['b']==b), 'status'].values[0] == 'Done'
    except IndexError:
        return False

