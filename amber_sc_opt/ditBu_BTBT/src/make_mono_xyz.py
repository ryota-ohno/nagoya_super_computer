import os
import subprocess
from utils import Rod, R2atom
import csv
import math

def concatenate(array_list):
    total_array=[]
    for arr in array_list:
        total_array.extend(arr)
    return total_array

def matmul(mat1, mat2):
    n = len(mat1)  # 行数
    result = [[0.0 for _ in range(3)] for _ in range(n)]
    for i in range(n):        # mat1 の各行
        for j in range(3):    # mat2 の各列
            for k in range(3):
                result[i][j] += mat1[i][k] * mat2[k][j]
    return result

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3,phi):
    csv_path=f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{monomer_name}/monomer/{monomer_name}.csv';cols=['X','Y','Z','R'];atoms_array_xyzR = []
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 指定列の値をfloatに変換してリスト化
            values = [float(row[col]) for col in cols]
            atoms_array_xyzR.append(values)
    xyz_array = [];R_array=[]
    for x,y,z,r in atoms_array_xyzR:
        xyz_array.append([x,y,z]);R_array.append(r)
    xyz_array = matmul(xyz_array,Rod([-1,0,0],A2))
    xyz_array = matmul(xyz_array,Rod([0,0,1],A3))
    
    C0_index = 5;C1_index = 35;C2_index = 13;C3_index = 22
    C0=xyz_array[C0_index];C1=xyz_array[C1_index]
    nx=C1[0]-C0[0];ny=C1[1]-C0[1];nz=C1[2]-C0[2]
    n1_=math.sqrt(nx**2+ny**2+nz**2)
    n1=[nx/n1_,ny/n1_,nz/n1_]
    
    C2=xyz_array[C2_index];C3=xyz_array[C3_index]
    nx=C3[0]-C2[0];ny=C3[1]-C2[1];nz=C3[2]-C2[2]
    n2_=math.sqrt(nx**2+ny**2+nz**2)
    n2=[nx/n2_,ny/n2_,nz/n2_]
    
    xyz_array_=[]; i=0
    for x,y,z in xyz_array:
        if i < 22: #0~21 BTBT
            xyz_array_.append([x,y,z])
        elif i > 21 and i < 35: #22~34 tBu 1
            xyz_array_1=[[x-C2[0],y-C2[1],z-C2[2]]]
            xyz_array_2=matmul(xyz_array_1,Rod(n2,phi))
            xyz_array_.append([xyz_array_2[0][0]+C2[0],xyz_array_2[0][1]+C2[1],xyz_array_2[0][2]+C2[2]])
        elif i > 34 : #35~47 tBu 2
            xyz_array_1=[[x-C0[0],y-C0[1],z-C0[2]]]
            xyz_array_2=matmul(xyz_array_1,Rod(n1,phi))
            xyz_array_.append([xyz_array_2[0][0]+C0[0],xyz_array_2[0][1]+C0[1],xyz_array_2[0][2]+C0[2]])
        i+=1

    xyz_array_f = []
    for x,y,z in xyz_array_:
        xyz_array_f.append([x+Ta,y+Tb,z+Tc])
    xyzR_array=[]
    for i in range(len(xyz_array_f)):
        xyzR_array.append([xyz_array_f[i][0],xyz_array_f[i][1],xyz_array_f[i][2],R_array[i]])
    return xyzR_array
        
line1='@<TRIPOS>MOLECULE\nmono-C9-BTBT\n   48    51     1     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n'
line2='@<TRIPOS>BOND\n'
bond_lines=[[1, 1, 2, 'ar'], [2, 1, 3, 'ar'], [3, 1, 7, '1'], [4, 2, 4, 'ar'], [5, 2, 9, '1'], [6, 3, 5, 'ar'], [7, 3, 20, '1'], [8, 4, 6, 'ar'], 
            [9, 4, 22, '1'], [10, 5, 6, 'ar'], [11, 5, 21, '1'], [12, 6, 36, '1'], [13, 7, 8, '2'], [14, 7, 10, '1'], [15, 8, 9, '1'], [16, 8, 12, '1'], 
            [17, 10, 11, '1'], [18, 11, 12, 'ar'], [19, 11, 13, 'ar'], [20, 12, 15, 'ar'], [21, 13, 14, 'ar'], [22, 13, 17, '1'], [23, 14, 16, 'ar'], [24, 14, 23, '1'], 
            [25, 15, 16, 'ar'], [26, 15, 19, '1'], [27, 16, 18, '1'], [28, 23, 24, '1'], [29, 23, 28, '1'], [30, 23, 32, '1'], [31, 24, 25, '1'], [32, 24, 26, '1'], 
            [33, 24, 27, '1'], [34, 28, 29, '1'], [35, 28, 30, '1'], [36, 28, 31, '1'], [37, 32, 33, '1'], [38, 32, 34, '1'], [39, 32, 35, '1'], [40, 36, 37, '1'], 
            [41, 36, 41, '1'], [42, 36, 45, '1'], [43, 37, 38, '1'], [44, 37, 39, '1'], [45, 37, 40, '1'], [46, 41, 42, '1'], [47, 41, 43, '1'], [48, 41, 44, '1'], 
            [49, 45, 46, '1'], [50, 45, 47, '1'], [51, 45, 48, '1']]
line3='@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n\n'

para_list=[]
with open(r'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/mono-C9-BTBT/monomer/mono_esp.mol2')as f:
    for line in f:
        #print(line)
        s=line.split()
        if len(s)==9:
            para_list.append([s[5],float(s[8])])
        if (line.find('BOND')>-1):
            break

def get_xyzR_lines(xyzr_array):
    lines=[]
    lines.append(line1)
    mol=int(len(xyzr_array))
    for i in range(mol):
        x,y,z,r=xyzr_array[i]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1} {R2atom(r)} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    lines.append(line2)
    for bond,atom1,atom2,type in bond_lines:
        line=f'{bond} {atom1} {atom2} {type}\n'
        lines.append(line)
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name):
    file_basename = os.path.splitext(file_name)[0]
    lines_job=[
'#!/bin/bash\n','\n',
'module load amber\n','\n',
f'antechamber -i {file_basename}.mol2 -fi mol2 -o {file_basename}_.mol2 -fo mol2 -s 2\n',
#f'parmchk2 -i {file_basename}_.mol2 -f mol2 -o {file_basename}.frcmod\n',
f'tleap -f {file_basename}_tleap.in\n',
f'sander -O -i FF_calc.in -o {file_basename}.out -p {file_basename}.prmtop -c {file_basename}.inpcrd -r min.rst -ref {file_basename}.inpcrd\n',
f'rm {file_basename}.inpcrd\n',
f'rm {file_basename}.prmtop\n',
]    
    lines_tleap=['source /home/center/opt/aarch64/apps/amber/19.0/dat/leap/cmd/leaprc.gaff\n',
f'MOL = loadmol2 {file_basename}.mol2\n',
f'loadamberparams ditBu_BTBT_dimer.frcmod\n',
f'saveamberparm MOL {file_basename}.prmtop {file_basename}.inpcrd\n',
'quit\n']
    file_job = os.path.join(auto_dir,f'amber/job_{file_basename}.sh')
    file_tleap = os.path.join(auto_dir,f'amber/{file_basename}_tleap.in')
    
    with open(file_job,'w')as f:
        f.writelines(lines_job)
    with open(file_tleap,'w')as f:
        f.writelines(lines_tleap)

    return file_job,f'{file_basename}.out'

def make_xyzfile(monomer_name,params_dict):
    phi = float(params_dict.get('phi',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,0,0,phi)
    xyz_list=[]
    for x,y,z,R in monomer_array_i:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['phi']:
            val = int(val)
        xyzfile_name += '_mono_{}'.format(val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):
    phi = float(params_dict.get('phi',0.0))

    monomer_array = get_monomer_xyzR(monomer_name,0,0,0,0,0,phi)
    line_list_monomer = get_xyzR_lines(monomer_array)
    gij_xyz_lines = line_list_monomer 
    
    file_name = get_file_name_from_dict(monomer_name,params_dict)
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'amber',file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['phi']:
            val = int(val)
        xyzfile_name += '_mono_{}'.format(val)
    return xyzfile_name + '.mol2'
    
def exec_gjf_mono(auto_dir, monomer_name, params_dict,isTest=True):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict)
    file_job,log_file_name = get_one_exe(auto_dir,file_name)
    if not(isTest):
        subprocess.run(['chmod','+x',file_job])
        subprocess.run([file_job])
    return log_file_name
    
############################################################################################