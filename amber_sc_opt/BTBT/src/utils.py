import math

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
    theta_t=math.radians(theta_in)
    Rod=[[math.cos(theta_t)+(nx**2)*(1-math.cos(theta_t)),nx*ny*(1-math.cos(theta_t))+nz*math.sin(theta_t),nx*nz*(1-math.cos(theta_t))-ny*math.sin(theta_t)],
        [nx*ny*(1-math.cos(theta_t))-nz*math.sin(theta_t),math.cos(theta_t)+(ny**2)*(1-math.cos(theta_t)),ny*nz*(1-math.cos(theta_t))+nx*math.sin(theta_t)],
        [nx*nz*(1-math.cos(theta_t))+ny*math.sin(theta_t),ny*nz*(1-math.cos(theta_t))-nx*math.sin(theta_t),math.cos(theta_t)+(nz**2)*(1-math.cos(theta_t))]]
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

