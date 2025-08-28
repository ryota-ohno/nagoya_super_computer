##tetracene層内計算
import os
os.environ['HOME'] ='/home/HasegawaLab'
import time
from make_xyz import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse
import shutil
import csv
## new params ('z' and 'A2') are used

def main_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT_step3/{args.auto_dir}'
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'amber'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    amber_path=os.path.join(auto_dir,'amber')
    shutil.copy(f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{args.monomer_name}_step3/src/FF_calc.in',amber_path)
    
    auto_csv_path = os.path.join(auto_dir,'step3.csv')
    if not os.path.exists(auto_csv_path): 
        header = ['cx','cy','cz','theta','A2','a','b','z','E','E1','E2','E3','E4','E5','E6','E7','E8','E9','status','file_name']
        with open(auto_csv_path, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=header)
            writer.writeheader()
        
    os.chdir(os.path.join(auto_dir,'amber'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(auto_dir,args.monomer_name,args.num_nodes,args.isTest)##argsの中身を取る
        #time.sleep(0.1)

def listen(auto_dir,monomer_name,num_nodes,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    fixed_param_keys = ['theta','A2','a','b','z'];opt_param_keys_1 = ['cx','cy','cz']
    
    mono_file=f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{monomer_name}/monomer/{monomer_name}_mono.out'
    E_mono=get_E(mono_file)[0]
    
    auto_csv_1 = os.path.join(auto_dir,'step3.csv');rows_1 = []
    with open(auto_csv_1, mode='r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows_1.append(row)
    for idx, row in enumerate(rows_1):
        if row['status'] != 'InProgress':
            continue
        params_dict1_ = {key: row[key] for key in fixed_param_keys + opt_param_keys_1 + ['file_name']}
        file_name1 = params_dict1_['file_name'];E_list=[]
        for i in range(1,10):
            log_filepath = os.path.join(auto_dir, 'amber', file_name1+f'_{i}.out')
            if not os.path.exists(log_filepath):
                continue
            E_list_ = get_E(log_filepath)
            if len(E_list_) != 1:
                continue
            E = round(float(E_list_[0]) - 2 * E_mono, 4)
            E_list.append(E)
        if len(E_list)!=9:
            continue
        E_columns=['E1','E2','E3','E4','E5','E6','E7','E8','E9']
        for key, value in zip(E_columns, E_list):
            rows_1[idx][key] = f"{value:.4f}"
        rows_1[idx]['E'] = f"{sum(E_list):.4f}";rows_1[idx]['status'] = 'Done'
        with open(auto_csv_1, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows_1[0].keys())
            writer.writeheader()
            writer.writerows(rows_1)
        break  # 1件で処理終了

    dict_matrix = get_params_dict(auto_dir,num_nodes)##更新分を流す a1/HOME/HASEGAWALABz2まで取得
    if len(dict_matrix)!=0:#終わりがまだ見えないなら
        for i in range(len(dict_matrix)):
            params_dict=dict_matrix[i]#print(params_dict)
            params_dict1 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys_1)}
            alreadyCalculated = check_calc_status(auto_dir,params_dict)
            if not(alreadyCalculated):
                dictlist_E1,fieldnames1 = read_csv_to_dictlist(os.path.join(auto_dir,'step3.csv'))
                dict_list_E1_f=filter_dictlist(dictlist_E1, params_dict1)
                if len(dict_list_E1_f) == 0:
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict1},isTest=isTest)
                    new_dict_1 = {**params_dict1,'E':'0.0','E1':'0.0','E2':'0.0','E3':'0.0','E4':'0.0','E5':'0.0','E6':'0.0','E7':'0.0','E8':'0.0','E9':'0.0', 'status': 'InProgress','file_name':file_name}
                    dictlist_E1.append(new_dict_1);write_dictlist_to_csv(os.path.join(auto_dir,'step3.csv'), dictlist_E1, fieldnames1)
    
    dictlist_init_params,field_name_i=read_csv_to_dictlist(os.path.join(auto_dir, 'step3_init_params.csv'))
    dictlist_init_params_done = filter_dictlist(dictlist_init_params,{'status':'Done'})
    isOver = True if len(dictlist_init_params_done)==len(dictlist_init_params) else False
    return isOver

def check_calc_status(auto_dir,params_dict):
    dict_list, fieldnames = read_csv_to_dictlist(os.path.join(auto_dir,'step3.csv'))
    if len(dict_list)==0:
        return False
    dict_list_f=filter_dictlist(dict_list, params_dict)
    if len(dict_list_f)==0:
        return False
    try:
        status = dict_list_f[0]['status']
        return status=='Done'
    except KeyError:
        return False

def get_params_dict(auto_dir, num_nodes):
    init_params_csv = os.path.join(auto_dir, 'step3_init_params.csv')
    dictlist_init_params,fieldnames = read_csv_to_dictlist(init_params_csv)
    dictlist_init_params_inprogress = filter_dictlist(dictlist_init_params, {'status':'InProgress'})
    fixed_param_keys = ['theta','A2','a','b','z'];opt_param_keys_1 = ['cx','cy','cz']

    #最初の立ち上がり時
    dictlist_init_params_notyet = filter_dictlist(dictlist_init_params, {'status':'NotYet'})
    if len(dictlist_init_params_notyet) != 0:
        if len(dictlist_init_params_inprogress) < num_nodes:
            notyet_row = dictlist_init_params_notyet[0]
            notyet_index = dictlist_init_params.index(notyet_row)
            dictlist_init_params = update_row_value(dictlist_init_params, notyet_index, 'status', 'InProgress')
            write_dictlist_to_csv(init_params_csv, dictlist_init_params, fieldnames)
            selected_keys = fixed_param_keys + opt_param_keys_1
            params_dict = {k: notyet_row[k] for k in selected_keys if k in notyet_row}
            return [params_dict]
    init_params_csv = os.path.join(auto_dir, 'step3_init_params.csv');dictlist_init_params,_ = read_csv_to_dictlist(init_params_csv)
    
    dict_matrix = []
    for index, row in enumerate(dictlist_init_params):
        if row['status'] != 'InProgress':
            continue
        print(index)
        dictlist_init_params,_ = read_csv_to_dictlist(init_params_csv)
        init_params_dict = get_values_from_dictlist(dictlist_init_params, index, fixed_param_keys + opt_param_keys_1)
        fixed_params_dict = get_values_from_dictlist(dictlist_init_params, index, fixed_param_keys)
        isDone, opt_params_matrix = get_opt_params_dict(auto_dir, init_params_dict, fixed_params_dict)
        with open(os.path.join(auto_dir, 'debug4.txt'),'w') as f:
            f.write(f'debug4 {isDone} {len(opt_params_matrix)}')
        if isDone:
            opt_params_dict = {'cx': round(opt_params_matrix[0][0], 1),'cy': round(opt_params_matrix[0][1], 1),'cz': round(opt_params_matrix[0][2], 1)}
            dictlist_init_params = update_row_value(dictlist_init_params, index, 'status', 'Done')
            if index + 1 >= len(dictlist_init_params):
                status = 'Done'
            else:
                status = dictlist_init_params[index + 1]['status']
            fieldnames = list(dictlist_init_params[0].keys())
            write_dictlist_to_csv(init_params_csv, dictlist_init_params, fieldnames)
            if status == 'NotYet':
                opt_params_dict = get_values_from_dictlist(dictlist_init_params, index + 1, opt_param_keys_1)
                dictlist_init_params = update_row_value(dictlist_init_params, index + 1, 'status', 'InProgress')
                write_dictlist_to_csv(init_params_csv, dictlist_init_params, fieldnames)
                dict_matrix.append({**fixed_params_dict, **opt_params_dict})
            else:
                continue
        else:
            for i in range(len(opt_params_matrix)):
                opt_params_dict = {'cx': round(opt_params_matrix[i][0], 1),'cy': round(opt_params_matrix[i][1], 1),'cz': round(opt_params_matrix[i][2], 1)}
                d = {**fixed_params_dict, **opt_params_dict};dict_matrix.append(d)
    return dict_matrix
        
def get_opt_params_dict(auto_dir, init_params_dict, fixed_params_dict):
    cur_csv = os.path.join(auto_dir, 'step1.csv');dictlist_cur,_ = read_csv_to_dictlist(cur_csv)
    filtered = filter_dictlist(dictlist_cur, fixed_params_dict)
    cx_init_prev = round(float(init_params_dict['cx']), 1);cy_init_prev = round(float(init_params_dict['cy']), 1);cz_init_prev = round(float(init_params_dict['cz']), 1)
    while True:
        E_list = [];xyz_list = [];para_list = []
        for dcx in [0.0]:
            for dcy in [-0.1, 0.0, 0.1]:
                for dcz in [-0.1, 0.0, 0.1]:
                    cx = round(cx_init_prev + dcx, 1);cy = round(cy_init_prev + dcy, 1);cz = round(cz_init_prev + dcz, 1);E = find_entry(filtered, cx, cy, cz)
                    if E is None:
                        para_list.append([cx, cy, cz])
                    else:
                        xyz_list.append([cx, cy, cz]);E_list.append(E)
        if len(para_list) != 0:
            return False, para_list
        min_idx = int(argmin(E_list));cx_init, cy_init, cz_init = xyz_list[min_idx]
        if cx_init == cx_init_prev and cy_init == cy_init_prev and cz_init == cz_init_prev:
            return True, [[cx_init, cy_init, cz_init]]
        else:
            cx_init_prev = cx_init;cy_init_prev = cy_init;cz_init_prev = cz_init

def filter_df(df, dict_filter):
    for k, v in dict_filter.items():
        if type(v)==str:
            df=df[df[k]==v]
        else:
            df=df[df[k]==v]
    df_filtered=df
    return df_filtered

def read_csv_to_dictlist(path):
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        return list(reader), reader.fieldnames

def write_dictlist_to_csv(path, dict_list, fieldnames):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(dict_list)

def filter_dictlist(dict_list, filter_dict):
    result = []
    for row in dict_list:
        if all(row.get(k) == str(v) for k, v in filter_dict.items()):
            result.append(row)
    return result

def find_entry(dict_list, cx, cy, cz):
    for row in dict_list:
        if (float(row['cx']) == cx and float(row['cy']) == cy and float(row['cz']) == cz and row['status'] == 'Done'):
            return float(row['E'])
    return None

def filter_data(data, dict_filter):
    filtered_data = []
    for row in data:
        match = True
        for key, value in dict_filter.items():
            if row.get(key) != value:
                match = False
                break
        if match:
            filtered_data.append(float(row))
    return filtered_data

def update_row_value(dictlist, row_index, key, new_value):
    dictlist[row_index][key] = new_value
    return dictlist

def get_values_from_dictlist(dictlist, index, keys):
    return {k: float(dictlist[index][k]) for k in keys}

def argmin(lst):
    if not lst:
        raise ValueError("Empty list has no minimum.")
    min_index = 0
    min_value = lst[0]
    for i, value in enumerate(lst):
        if value < min_value:
            min_index = i
            min_value = value
    return min_index

def final_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{args.monomer_name}_step3/{args.auto_dir}'
    with open(os.path.join(auto_dir,'done.txt'),'w')as f:
        f.write('Done')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes amber, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    args = parser.parse_args()

    print("----main process----")
    main_process(args)
    final_process(args)
    print("----finish process----")
    