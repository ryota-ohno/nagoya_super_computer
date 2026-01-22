import os
import subprocess
import csv
import time

def concatenate(array_list):
    total_array=[]
    for arr in array_list:
        total_array.extend(arr)
    return total_array

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc):
    csv_path=f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/PDI/monomer/{monomer_name}.csv';cols=['atom','X','Y','Z'];atoms_array_axyz = []
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 指定列の値をfloatに変換してリスト化
            values = [float(row[col]) for col in cols]
            atoms_array_axyz.append(values)
    xyz_array = [];a_array=[]
    for atom,x,y,z, in atoms_array_axyz:
        xyz_array.append([x,y,z]);a_array.append(atom)
    xyz_array_ = []
    for x,y,z in xyz_array:
        xyz_array_.append([x+Ta,y+Tb,z+Tc])
    axyz_array=[]
    for i in range(len(xyz_array_)):
        axyz_array.append([a_array[i],xyz_array_[i][0],xyz_array_[i][1],xyz_array_[i][2]])
    return axyz_array

## PDI num_atom 40 num_bond 46 

line1='@<TRIPOS>MOLECULE\nPDI\n   80    92     2     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n'
line2='@<TRIPOS>BOND\n'


###変更
bond_lines=[[1, 1, 2, 'ar'], [2, 1, 3, 'ar'], [3, 1, 7, 'ar'], [4, 2, 4, 'ar'], [5, 2, 9, 'ar'], [6, 3, 5, 'ar'], [7, 3, 22, '1'], [8, 4, 6, 'ar'],
            [9, 4, 12, 'ar'], [10, 5, 6, 'ar'], [11, 5, 38, '1'], [12, 6, 37, '1'], [13, 7, 8, 'ar'], [14, 7, 21, '1'], [15, 8, 10, 'ar'], [16, 8, 31, '1'], 
            [17, 9, 10, 'ar'], [18, 9, 11, 'ar'], [19, 10, 32, '1'], [20, 11, 13, 'ar'], [21, 11, 18, 'ar'], [22, 12, 13, 'ar'], [23, 12, 14, 'ar'], [24, 13, 16, 'ar'], 
            [25, 14, 15, 'ar'], [26, 14, 36, '1'], [27, 15, 17, 'ar'], [28, 15, 35, '1'], [29, 16, 17, 'ar'], [30, 16, 20, 'ar'], [31, 17, 23, '1'], [32, 18, 19, 'ar'], 
            [33, 18, 33, '1'], [34, 19, 20, 'ar'], [35, 19, 34, '1'], [36, 20, 24, '1'], [37, 21, 25, '1'], [38, 21, 26, '2'], [39, 22, 25, '1'], [40, 22, 27, '2'], 
            [41, 23, 28, '1'], [42, 23, 30, '2'], [43, 24, 28, '1'], [44, 24, 29, '2'], [45, 25, 40, '1'], [46, 28, 39, '1']]

line3='@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n     2 RES2       41 GROUP             0 ****  ****    0 \n\n'

para_list=[]
with open(r'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/PDI/monomer/PDI_mono.mol2')as f:
    for line in f:
        #print(line)
        s=line.split()
        if len(s)==9:
            para_list.append([s[5],float(s[8])])
        if (line.find('BOND')>-1):
            break

def make_dimer_mol2(xyzr_list_1,xyzr_list_2):
    num_mol1=int(len(xyzr_list_1)/40);num_mol2=int(len(xyzr_list_2)/40)
    lines=[f'@<TRIPOS>MOLECULE\nPDI\n   {int(40*(num_mol1+num_mol2))}    {int(46*(num_mol1+num_mol2))}     1     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n']
    for i in range(num_mol1):
        for j in range(40):
            a,x,y,z=xyzr_list_1[40*i+j]
            atom_type,charge=para_list[j]
            lines.append(f'  {40*i+j+1} {a} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    for i in range(num_mol2):
        for j in range(40):
            a,x,y,z=xyzr_list_2[40*i+j]
            atom_type,charge=para_list[j]
            lines.append(f'  {40*(i+num_mol1)+j+1} {a} {x} {y} {z} {atom_type} 2 RES2 {charge}\n')
    lines.append('@<TRIPOS>BOND\n')
    for i in range(num_mol1+num_mol2):
        for j in range(46):
            a,b,c,type=bond_lines[j]
            a+=46*i;b+=40*i;c+=40*i
            lines.append(f'{a} {b} {c} {type} \n')
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name):
    file_basename = file_name
    lines_job=['#!/bin/bash\n','\n',
    'module load amber\n','\n',]
    lines_job.append(f'tleap -f {file_basename}_tleap.in\n')
    for i in range(1,4):
        file_basename_=file_basename+f'_{i}'
        lines_job.append(f'sander -O -i FF_calc.in -o {file_basename_}.out -p {file_basename_}.prmtop -c {file_basename_}.inpcrd -r min.rst -ref {file_basename_}.inpcrd\n')
        lines_job.append(f'rm {file_basename_}.inpcrd\n')
        lines_job.append(f'rm {file_basename_}.prmtop\n')
    
    lines_tleap=['source /home/center/opt/aarch64/apps/amber/19.0/dat/leap/cmd/leaprc.gaff\n']
    for i in range(1,4):
        file_basename_=file_basename+f'_{i}'
        lines_tleap.append(f'MOL = loadmol2 {file_basename_}.mol2\n')
        lines_tleap.append(f'loadamberparams PDI_mono.frcmod\n')
        lines_tleap.append(f'saveamberparm MOL {file_basename_}.prmtop {file_basename_}.inpcrd\n')
    lines_tleap.append('quit\n\n')
        
    file_tleap = os.path.join(auto_dir,f'amber/{file_basename}_tleap.in')
    with open(file_tleap,'w')as f:
        f.writelines(lines_tleap)
    file_job = os.path.join(auto_dir,f'amber/job_{file_basename}.sh')
    with open(file_job,'w')as f:
        f.writelines(lines_job)
    file_job_base = os.path.join(auto_dir,f'amber/job_{file_basename}')
    return file_job_base,f'{file_basename}'

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict):
    x = float(params_dict.get('x',0.0));y = float(params_dict.get('y',0.0)); z = float(params_dict.get('z',0.0))
    
    monomer_array_0 = get_monomer_xyzR(monomer_name,0,0,0)
    monomer_array_1 = get_monomer_xyzR(monomer_name,x,y,z)
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    monomers_array_4 = concatenate([monomer_array_0,monomer_array_1])
    
    for atom,x,y,z in monomers_array_4:
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['x','y','z']:
            val = round(val,1)
        xyzfile_name += '_{}_{}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):
    x = float(params_dict.get('x',0.0));y = float(params_dict.get('y',0.0)); z = float(params_dict.get('z',0.0))
    
    monomer_array_0 = get_monomer_xyzR(monomer_name,0,0,0)
    monomer_array_1 = get_monomer_xyzR(monomer_name,x,y,z)
    
    dimer_1=make_dimer_mol2(monomer_array_0,monomer_array_1)
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    file_mol2 = file_name+ '.mol2'
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'amber',file_mol2)
    line=dimer_1
    with open(gij_xyz_path,'w') as f:
        f.writelines(line)
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['x','y','z']:
            val = round(val,1)
        file_name += '_{}_{}'.format(key,val)
    return file_name
    
def exec_gjf(auto_dir, monomer_name, params_dict,isTest=True):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
    file_job_base,log_file_name = get_one_exe(auto_dir,file_name)
    file_job=file_job_base+'.sh'
    if not(isTest):
        subprocess.run(['chmod','+x',file_job])
        subprocess.run([file_job])
    time.sleep(0.1)
    return log_file_name

############################################################################################