import os
import numpy as np
import pandas as pd
import subprocess

def get_E1(path_file):
    with open(path_file,'r') as f:
        lines=f.readlines()
    lines_E=[]
    for line in lines:
        if line.find('E(R')>-1 and len(line.split())>5:
            lines_E.append(float(line.split()[4])*627.510)
    E_list=[lines_E[5*i]-lines_E[5*i+1]-lines_E[5*i+2] for i in range(int(len(lines_E)/5))]
    return E_list[0]

auto_dir='/data/group1/z40145w/Working/nagoya_super_computer/perylene/sandwich_step2_para/2/gaussian'
z_list1=[np.round(z,1) for z in np.linspace(np.round(-2.5,1),np.round(2.5,1),int(np.round(np.round(5,1)/0.1))+1)]
file_base_name='perylene_step2_a=11.3_b=10.6_theta=35_'
with open('/data/group1/z40145w/Working/nagoya_super_computer/perylene/sandwich_step2_para/2/result1.txt','w')as f:
    for z in z_list1:
        file_basename1 = file_base_name;file_basename2 = file_base_name;file_basename3 = file_base_name
        file_basename1 +='z={}_1.log'.format(z);file_basename2 +='z={}_2.log'.format(z);file_basename3 +='z={}_3.log'.format(z)
        result1=os.path.join(auto_dir,file_basename1);result2=os.path.join(auto_dir,file_basename2);result3=os.path.join(auto_dir,file_basename3)
        E1=get_E1(result1);E2=get_E1(result2);E3=get_E1(result3)
        f.write('{} {} {} {}\n'.format(z,E1,E2,E3))


############################################################################################