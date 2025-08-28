import os
import subprocess
from utils import Rod, R2atom
import csv

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

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3):
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
    xyz_array_ = []
    for x,y,z in xyz_array:
        xyz_array_.append([x+Ta,y+Tb,z+Tc])
    xyzR_array=[]
    for i in range(len(xyz_array_)):
        xyzR_array.append([xyz_array_[i][0],xyz_array_[i][1],xyz_array_[i][2],R_array[i]])
    return xyzR_array
        
line1='@<TRIPOS>MOLECULE\npentacene\n   48    54     2     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n'
line2='@<TRIPOS>BOND\n'
bond_lines=[' 1 1 2 ar\n', ' 2 1 3 ar\n', ' 3 1 7 1\n', ' 4 2 4 ar\n', ' 5 2 9 1\n', ' 6 3 5 ar\n', ' 7 3 21 1\n', ' 8 4 6 ar\n', ' 9 4 24 1\n',
            ' 10 5 6 ar\n', ' 11 5 22 1\n', ' 12 6 23 1\n', ' 13 7 8 2\n', ' 14 7 10 1\n', ' 15 8 9 1\n', ' 16 8 12 1\n', ' 17 10 11 1\n', ' 18 11 12 ar\n',
            ' 19 11 13 ar\n', ' 20 12 15 ar\n', ' 21 13 14 ar\n', ' 22 13 17 1\n', ' 23 14 16 ar\n', ' 24 14 18 1\n', ' 25 15 16 ar\n', ' 26 15 20 1\n', ' 27 16 19 1\n',
            ' 28 25 26 ar\n', ' 29 25 27 ar\n', ' 30 25 31 1\n', ' 31 26 28 ar\n', ' 32 26 33 1\n', ' 33 27 29 ar\n', ' 34 27 45 1\n', ' 35 28 30 ar\n', ' 36 28 48 1\n',
            ' 37 29 30 ar\n', ' 38 29 46 1\n', ' 39 30 47 1\n', ' 40 31 32 2\n', ' 41 31 34 1\n', ' 42 32 33 1\n', ' 43 32 36 1\n', ' 44 34 35 1\n', ' 45 35 36 ar\n',
            ' 46 35 37 ar\n', ' 47 36 39 ar\n', ' 48 37 38 ar\n', ' 49 37 41 1\n', ' 50 38 40 ar\n', ' 51 38 42 1\n', ' 52 39 40 ar\n', ' 53 39 44 1\n', ' 54 40 43 1\n']
line3='@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n     2 RES2       25 GROUP             0 ****  ****    0 \n\n'

para_list=[]
with open(r'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT/monomer/btbt_mono.mol2')as f:
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
    mol=int(len(xyzr_array)/2)
    for i in range(mol):
        x,y,z,r=xyzr_array[i]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1} {R2atom(r)} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    for i in range(mol):
        x,y,z,r=xyzr_array[i+mol]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1+mol} {R2atom(r)} {x} {y} {z} {atom_type} 2 RES2 {charge}\n')   
    lines.append(line2)
    for line in bond_lines:
        lines.append(line)
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name):
    file_basename = file_name
    for i in range(1,10):
        file_basename_=file_basename+f'_{i}'
        lines_job=[
    '#!/bin/bash\n','\n',
    'module load amber\n','\n',
    f'parmchk2 -i {file_basename_}.mol2 -f mol2 -o {file_basename_}.frcmod\n',
    f'tleap -f {file_basename_}_tleap.in\n',
    f'sander -O -i FF_calc.in -o {file_basename_}.out -p {file_basename_}.prmtop -c {file_basename_}.inpcrd -r min.rst -ref {file_basename_}.inpcrd\n',
    f'rm {file_basename_}.frcmod\n',
    f'rm {file_basename_}.inpcrd\n',
    f'rm {file_basename_}.prmtop\n',
    ]
        
        lines_tleap=['source /home/center/opt/aarch64/apps/amber/19.0/dat/leap/cmd/leaprc.gaff\n',
    f'MOL = loadmol2 {file_basename_}.mol2\n',
    f'loadamberparams {file_basename_}.frcmod\n',
    f'saveamberparm MOL {file_basename_}.prmtop {file_basename_}.inpcrd\n',
    'quit\n']
        file_job = os.path.join(auto_dir,f'amber/job_{file_basename_}.sh')
        file_tleap = os.path.join(auto_dir,f'amber/{file_basename_}_tleap.in')
        
        with open(file_job,'w')as f:
            f.writelines(lines_job)
        with open(file_tleap,'w')as f:
            f.writelines(lines_tleap)
    file_job_base = os.path.join(auto_dir,f'amber/job_{file_basename}')
    return file_job_base,f'{file_basename}'

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    cx = float(params_dict.get('cx',0.0));cy = float(params_dict.get('cy',0.0)); cz = float(params_dict.get('cz',0.0))
    
    monomer_array_c = get_monomer_xyzR(monomer_name,cx,cy,cz,A2,A3)
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,-a,0,0,A2,A3)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3)##1,2がb方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,0,-b,-2*z,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z,A2,-A3)##1,2がb方向
    monomer_array_t3 = get_monomer_xyzR(monomer_name,a/2,-b/2,-z,A2,-A3)##1,2がb方向
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    monomers_array_4 = concatenate([monomer_array_c,monomer_array_i,monomer_array_p1,monomer_array_p2,monomer_array_p3,monomer_array_p4,
                                    monomer_array_t1,monomer_array_t2,monomer_array_t3,monomer_array_t4])
    
    for x,y,z,R in monomers_array_4:
        atom = R2atom(R)
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z','cx','cy','cz']:
            val = round(val,1)
        elif key in ['A1','A2','theta']:
            val = int(val)
        xyzfile_name += '_{}_{}'.format(key,val)
    return xyzfile_name + '.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    cx = float(params_dict.get('cx',0.0));cy = float(params_dict.get('cy',0.0)); cz = float(params_dict.get('cz',0.0))
    
    monomer_array_c = get_monomer_xyzR(monomer_name,cx,cy,cz,A2,A3)
    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3)
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,-a,0,0,A2,A3)##1,2がb方向
    monomer_array_p3 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3)##1,2がb方向
    monomer_array_p4 = get_monomer_xyzR(monomer_name,0,-b,-2*z,A2,A3)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,b/2,z,A2,-A3)##1,2がb方向
    monomer_array_t3 = get_monomer_xyzR(monomer_name,a/2,-b/2,-z,A2,-A3)##1,2がb方向
    monomer_array_t4 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3)##1,2がb方向
    
    dimer_array_i = concatenate([monomer_array_c,monomer_array_i])
    dimer_array_p1 = concatenate([monomer_array_c,monomer_array_p1]);dimer_array_p2 = concatenate([monomer_array_c,monomer_array_p2]);dimer_array_p3 = concatenate([monomer_array_c,monomer_array_p3]);dimer_array_p4 = concatenate([monomer_array_c,monomer_array_p4])
    dimer_array_t1 = concatenate([monomer_array_c,monomer_array_t1]);dimer_array_t2 = concatenate([monomer_array_c,monomer_array_t2]);dimer_array_t3 = concatenate([monomer_array_c,monomer_array_t3]);dimer_array_t4 = concatenate([monomer_array_c,monomer_array_t4])
    
    line_i= get_xyzR_lines(dimer_array_i);line_p1= get_xyzR_lines(dimer_array_p1);line_p2= get_xyzR_lines(dimer_array_p2);line_p3= get_xyzR_lines(dimer_array_p3);line_p4= get_xyzR_lines(dimer_array_p4)
    line_t1= get_xyzR_lines(dimer_array_t1);line_t2= get_xyzR_lines(dimer_array_t2);line_t3= get_xyzR_lines(dimer_array_t3);line_t4= get_xyzR_lines(dimer_array_t4)
    
    lines=[line_i,line_p1,line_p2,line_p3,line_p4,line_t1,line_t2,line_t3,line_t4]
    i=1
    for line in lines:
        line = line + ['\n\n\n']
        file_name = get_file_name_from_dict(monomer_name,params_dict)
        file_mol2 = file_name+ f'_{i}.mol2'
        os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
        gij_xyz_path = os.path.join(auto_dir,'amber',file_mol2)
        with open(gij_xyz_path,'w') as f:
            f.writelines(line)
        i+=1
    return file_name

def get_file_name_from_dict(monomer_name,params_dict):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        if key in ['a','b','z','cx','cy','cz']:
            val = round(val,1)
        elif key in ['A2','theta']:
            val = int(val)
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
    for i in range(1,10):
        file_job=file_job_base+f'_{i}.sh'
        if not(isTest):
            subprocess.run(['chmod','+x',file_job])
            subprocess.run([file_job])
    return log_file_name
    
############################################################################################