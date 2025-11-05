##tetracene層内計算
import os
os.environ['HOME'] ='/home/HasegawaLab'
import time
from make_8_xyz import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse
import shutil
import csv
## new params ('z' and 'A2') are used

def main_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/pentacene/{args.auto_dir}'
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'amber'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    amber_path=os.path.join(auto_dir,'amber')
    shutil.copy(f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{args.monomer_name}/src/FF_calc.in',amber_path)
    
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path): 
        header = ['theta','a','b','z1','z2','E','E1','E2','E3','status']
        with open(auto_csv_path, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=header)
            writer.writeheader()
        
    auto_csv_path1 = os.path.join(auto_dir,'step1_1.csv')
    if not os.path.exists(auto_csv_path1):
        header = ['theta','a','z1','z2','E1','status','file_name']
        with open(auto_csv_path1, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=header)
            writer.writeheader()
        
    auto_csv_path2 = os.path.join(auto_dir,'step1_2.csv')
    if not os.path.exists(auto_csv_path2):
        header = ['theta','b','z1','z2','E2','status','file_name']
        with open(auto_csv_path2, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=header)
            writer.writeheader()
        
    auto_csv_path3 = os.path.join(auto_dir,'step1_3.csv')
    if not os.path.exists(auto_csv_path3):        
        header = ['theta','a','b','z1','z2','E3','status','file_name']
        with open(auto_csv_path3, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=header)
            writer.writeheader()

    auto_csv_path4 = os.path.join(auto_dir,'step1_4.csv')
    if not os.path.exists(auto_csv_path4):        
        header = ['theta','a','b','z1','z2','E4','status','file_name']
        with open(auto_csv_path4, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=header)
            writer.writeheader()

    os.chdir(os.path.join(auto_dir,'amber'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(auto_dir,args.monomer_name,args.num_nodes,args.isTest)##argsの中身を取る
        #time.sleep(0.1)

def listen(auto_dir,monomer_name,num_nodes,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    fixed_param_keys = ['theta','z1','z2'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b']
    
    mono_file=f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{monomer_name}/monomer/{monomer_name}_mono.out'
    E_mono=get_E(mono_file)[0]
    
    auto_csv_1 = os.path.join(auto_dir,'step1_1.csv');rows_1 = []
    with open(auto_csv_1, mode='r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows_1.append(row)
    for idx, row in enumerate(rows_1):
        if row['status'] != 'InProgress':
            continue
        params_dict1_ = {key: row[key] for key in fixed_param_keys + opt_param_keys_1 + ['file_name']}
        file_name1 = params_dict1_['file_name'];log_filepath1 = os.path.join(auto_dir, 'amber', file_name1)
        if not os.path.exists(log_filepath1):
            continue
        E_list1 = get_E(log_filepath1)
        if len(E_list1) != 1:
            continue
        E1 = round(float(E_list1[0]) - 2 * E_mono, 4)
        rows_1[idx]['E1'] = E1;rows_1[idx]['status'] = 'Done'
        with open(auto_csv_1, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows_1[0].keys())
            writer.writeheader()
            writer.writerows(rows_1)
        break  # 1件で処理終了
    
    auto_csv_2 = os.path.join(auto_dir,'step1_2.csv');rows_2 = []
    with open(auto_csv_2, mode='r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows_2.append(row)
    for idx, row in enumerate(rows_2):
        if row['status'] != 'InProgress':
            continue
        params_dict2_ = {key: row[key] for key in fixed_param_keys + opt_param_keys_2 + ['file_name']}
        file_name2 = params_dict2_['file_name'];log_filepath2 = os.path.join(auto_dir, 'amber', file_name2)
        if not os.path.exists(log_filepath2):
            continue
        E_list2 = get_E(log_filepath2)
        if len(E_list2) != 1:
            continue
        E2 = round(float(E_list2[0]) - 2 * E_mono, 4)
        rows_2[idx]['E2'] = E2;rows_2[idx]['status'] = 'Done'
        with open(auto_csv_2, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows_2[0].keys())
            writer.writeheader()
            writer.writerows(rows_2)
        break  # 1件で処理終了
    
    auto_csv_3 = os.path.join(auto_dir,'step1_3.csv');rows_3 = []
    with open(auto_csv_3, mode='r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows_3.append(row)
    for idx, row in enumerate(rows_3):
        if row['status'] != 'InProgress':
            continue
        params_dict3_ = {key: row[key] for key in fixed_param_keys + opt_param_keys_1 + opt_param_keys_2 + ['file_name']}
        file_name3 = params_dict3_['file_name'];log_filepath3 = os.path.join(auto_dir, 'amber', file_name3)
        if not os.path.exists(log_filepath3):
            continue
        E_list3 = get_E(log_filepath3)
        if len(E_list3) != 1:
            continue
        E3 = round(float(E_list3[0]) - 2 * E_mono, 4)
        rows_3[idx]['E3'] = E3;rows_3[idx]['status'] = 'Done'
        with open(auto_csv_3, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows_3[0].keys())
            writer.writeheader()
            writer.writerows(rows_3)
        break  # 1件で処理終了

    auto_csv_4 = os.path.join(auto_dir,'step1_4.csv');rows_4 = []
    with open(auto_csv_4, mode='r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows_4.append(row)
    for idx, row in enumerate(rows_4):
        if row['status'] != 'InProgress':
            continue
        params_dict4_ = {key: row[key] for key in fixed_param_keys + opt_param_keys_1 + opt_param_keys_2 + ['file_name']}
        file_name4 = params_dict4_['file_name'];log_filepath4 = os.path.join(auto_dir, 'amber', file_name4)
        if not os.path.exists(log_filepath4):
            continue
        E_list4 = get_E(log_filepath4)
        if len(E_list4) != 1:
            continue
        E4 = round(float(E_list4[0]) - 2 * E_mono, 4)
        rows_4[idx]['E4'] = E4;rows_4[idx]['status'] = 'Done'
        with open(auto_csv_4, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows_4[0].keys())
            writer.writeheader()
            writer.writerows(rows_4)
        break  # 1件で処理終了

    auto_csv = os.path.join(auto_dir, 'step1.csv')
    df, fieldnames = read_csv_to_dictlist(auto_csv)
    df_1, _ = read_csv_to_dictlist(os.path.join(auto_dir, 'step1_1.csv'))
    df_2, _ = read_csv_to_dictlist(os.path.join(auto_dir, 'step1_2.csv'))
    df_3, _ = read_csv_to_dictlist(os.path.join(auto_dir, 'step1_3.csv'))
    df_4, _ = read_csv_to_dictlist(os.path.join(auto_dir, 'step1_4.csv'))
    updated = False
    for row in df:
        if row['status'] != 'InProgress':
            continue
        params_dict1_ = {k: row[k] for k in fixed_param_keys + opt_param_keys_1}
        params_dict2_ = {k: row[k] for k in fixed_param_keys + opt_param_keys_2}
        params_dict3_ = {k: row[k] for k in fixed_param_keys + opt_param_keys_1 + opt_param_keys_2}
        params_dict4_ = {k: row[k] for k in fixed_param_keys + opt_param_keys_1 + opt_param_keys_2}
        s1 = [r for r in filter_dictlist(df_1, params_dict1_) if r['status'] == 'Done']
        s2 = [r for r in filter_dictlist(df_2, params_dict2_) if r['status'] == 'Done']
        s3 = [r for r in filter_dictlist(df_3, params_dict3_) if r['status'] == 'Done']
        s4 = [r for r in filter_dictlist(df_4, params_dict4_) if r['status'] == 'Done']
        if not s1 or not s2 or not s3 or not s4:
            continue
        E1 = float(s1[0]['E1']);E2 = float(s2[0]['E2']);E3 = float(s3[0]['E3']);E4 = float(s4[0]['E4'])
        E = 2 * E1 + 2 * E2 + 2 * E3 + 2 * E4
        row['E1'] = f"{E1:.4f}";row['E2'] = f"{E2:.4f}";row['E3'] = f"{E3:.4f}";row['E4'] = f"{E4:.4f}"
        row['E'] = f"{E:.4f}";row['status'] = 'Done'
        updated = True
        break
    if updated:
        write_dictlist_to_csv(auto_csv, df, fieldnames)

    dict_matrix = get_params_dict(auto_dir,num_nodes)##更新分を流す a1/HOME/HASEGAWALABz2まで取得
    if len(dict_matrix)!=0:#終わりがまだ見えないなら
        for i in range(len(dict_matrix)):
            params_dict=dict_matrix[i]#print(params_dict)
            params_dict1 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys_1)}
            params_dict2 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys_2)}
            params_dict3 = params_dict;params_dict4 = params_dict
            alreadyCalculated = check_calc_status(auto_dir,params_dict)
            if not(alreadyCalculated):
                dictlist_E,fieldnames = read_csv_to_dictlist(os.path.join(auto_dir,'step1.csv'))
                dictlist_E_f=filter_dictlist(dictlist_E, params_dict)
                if len(dictlist_E_f) == 0:
                    new_dict = {**params_dict, 'E': '0.0', 'E1': '0.0', 'E2': '0.0', 'E3': '0.0', 'E4': '0.0', 'status': 'InProgress'}
                    dictlist_E.append(new_dict);write_dictlist_to_csv(os.path.join(auto_dir,'step1.csv'), dictlist_E, fieldnames)

                ## 1の実行　##
                dictlist_E1,fieldnames1 = read_csv_to_dictlist(os.path.join(auto_dir,'step1_1.csv'))
                dict_list_E1_f=filter_dictlist(dictlist_E1, params_dict1)
                if len(dict_list_E1_f) == 0:
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict1}, structure_type=1,isTest=isTest)
                    new_dict_1 = {**params_dict1, 'E1': '0.0', 'status': 'InProgress','file_name':file_name}
                    dictlist_E1.append(new_dict_1);write_dictlist_to_csv(os.path.join(auto_dir,'step1_1.csv'), dictlist_E1, fieldnames1)
                    
                ## 2の実行　##
                dictlist_E2,fieldnames2 = read_csv_to_dictlist(os.path.join(auto_dir,'step1_2.csv'))
                dict_list_E2_f=filter_dictlist(dictlist_E2, params_dict2)
                if len(dict_list_E2_f) == 0:
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict2}, structure_type=2,isTest=isTest)
                    new_dict_2 = {**params_dict2, 'E2': '0.0', 'status': 'InProgress','file_name':file_name}
                    dictlist_E2.append(new_dict_2);write_dictlist_to_csv(os.path.join(auto_dir,'step1_2.csv'), dictlist_E2, fieldnames2)
                    
                ## 3の実行　##
                dictlist_E3,fieldnames3 = read_csv_to_dictlist(os.path.join(auto_dir,'step1_3.csv'))
                dict_list_E3_f=filter_dictlist(dictlist_E3, params_dict3)
                if len(dict_list_E3_f) == 0:
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict3}, structure_type=3,isTest=isTest)
                    new_dict_3 = {**params_dict3, 'E3': '0.0', 'status': 'InProgress','file_name':file_name}
                    dictlist_E3.append(new_dict_3);write_dictlist_to_csv(os.path.join(auto_dir,'step1_3.csv'), dictlist_E3, fieldnames3)
                
                dictlist_E4,fieldnames4 = read_csv_to_dictlist(os.path.join(auto_dir,'step1_4.csv'))
                dict_list_E4_f=filter_dictlist(dictlist_E4, params_dict4)
                if len(dict_list_E4_f) == 0:
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict4}, structure_type=4,isTest=isTest)
                    new_dict_4 = {**params_dict4, 'E4': '0.0', 'status': 'InProgress','file_name':file_name}
                    dictlist_E4.append(new_dict_4);write_dictlist_to_csv(os.path.join(auto_dir,'step1_4.csv'), dictlist_E4, fieldnames4)
    
    dictlist_init_params,field_name_i=read_csv_to_dictlist(os.path.join(auto_dir, 'step1_init_params.csv'))
    dictlist_init_params_done = filter_dictlist(dictlist_init_params,{'status':'Done'})
    isOver = True if len(dictlist_init_params_done)==len(dictlist_init_params) else False
    return isOver

def check_calc_status(auto_dir,params_dict):
    dict_list, fieldnames = read_csv_to_dictlist(os.path.join(auto_dir,'step1.csv'))
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
    init_params_csv = os.path.join(auto_dir, 'step1_init_params.csv')
    dictlist_init_params,fieldnames = read_csv_to_dictlist(init_params_csv)
    dictlist_init_params_inprogress = filter_dictlist(dictlist_init_params, {'status':'InProgress'})
    fixed_param_keys = ['theta','z1','z2'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b']

    #最初の立ち上がり時
    dictlist_init_params_notyet = filter_dictlist(dictlist_init_params, {'status':'NotYet'})
    if len(dictlist_init_params_notyet) != 0:
        if len(dictlist_init_params_inprogress) < num_nodes:
            notyet_row = dictlist_init_params_notyet[0]
            notyet_index = dictlist_init_params.index(notyet_row)
            dictlist_init_params = update_row_value(dictlist_init_params, notyet_index, 'status', 'InProgress')
            write_dictlist_to_csv(init_params_csv, dictlist_init_params, fieldnames)
            selected_keys = fixed_param_keys + opt_param_keys_1 + opt_param_keys_2
            params_dict = {k: notyet_row[k] for k in selected_keys if k in notyet_row}
            return [params_dict]
    init_params_csv = os.path.join(auto_dir, 'step1_init_params.csv');dictlist_init_params,_ = read_csv_to_dictlist(init_params_csv)
    
    dict_matrix = []
    for index, row in enumerate(dictlist_init_params):
        if row['status'] != 'InProgress':
            continue
        print(index)
        dictlist_init_params,_ = read_csv_to_dictlist(init_params_csv)
        init_params_dict = get_values_from_dictlist(dictlist_init_params, index, fixed_param_keys + opt_param_keys_1 + opt_param_keys_2)
        fixed_params_dict = get_values_from_dictlist(dictlist_init_params, index, fixed_param_keys)
        isDone, opt_params_matrix = get_opt_params_dict(auto_dir, init_params_dict, fixed_params_dict)
        with open(os.path.join(auto_dir, 'debug4.txt'),'w') as f:
            f.write(f'debug4 {isDone} {len(opt_params_matrix)}')
        if isDone:
            opt_params_dict = {'a': round(opt_params_matrix[0][0], 1),'b': round(opt_params_matrix[0][1], 1)}
            dictlist_init_params = update_row_value(dictlist_init_params, index, 'status', 'Done')
            if index + 1 >= len(dictlist_init_params):
                status = 'Done'
            else:
                status = dictlist_init_params[index + 1]['status']
            fieldnames = list(dictlist_init_params[0].keys())
            write_dictlist_to_csv(init_params_csv, dictlist_init_params, fieldnames)
            if status == 'NotYet':
                opt_params_dict = get_values_from_dictlist(dictlist_init_params, index + 1, opt_param_keys_1 + opt_param_keys_2)
                dictlist_init_params = update_row_value(dictlist_init_params, index + 1, 'status', 'InProgress')
                write_dictlist_to_csv(init_params_csv, dictlist_init_params, fieldnames)
                dict_matrix.append({**fixed_params_dict, **opt_params_dict})
            else:
                continue
        else:
            for i in range(len(opt_params_matrix)):
                opt_params_dict = {'a': round(opt_params_matrix[i][0], 1),'b': round(opt_params_matrix[i][1], 1)}
                d = {**fixed_params_dict, **opt_params_dict};dict_matrix.append(d)
    return dict_matrix
        
def get_opt_params_dict(auto_dir, init_params_dict, fixed_params_dict):
    cur_csv = os.path.join(auto_dir, 'step1.csv');dictlist_cur,_ = read_csv_to_dictlist(cur_csv)
    filtered = filter_dictlist(dictlist_cur, fixed_params_dict)
    a_init_prev = round(float(init_params_dict['a']), 1);b_init_prev = round(float(init_params_dict['b']), 1)
    while True:
        E_list = [];xyz_list = [];para_list = []
        for da in [-0.1, 0.0, 0.1]:
            for db in [-0.1, 0.0, 0.1]:
                a = round(a_init_prev + da, 1);b = round(b_init_prev + db, 1);E = find_entry(filtered, a, b)
                if E is None:
                    para_list.append([a, b])
                else:
                    xyz_list.append([a, b]);E_list.append(E)
        if len(para_list) != 0:
            return False, para_list
        min_idx = int(argmin(E_list));a_init, b_init = xyz_list[min_idx]
        if a_init == a_init_prev and b_init == b_init_prev:
            return True, [[a_init, b_init]]
        else:
            a_init_prev = a_init;b_init_prev = b_init

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

def find_entry(dict_list, a, b):
    for row in dict_list:
        if (float(row['a']) == a and float(row['b']) == b and row['status'] == 'Done'):
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
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/{args.monomer_name}/{args.auto_dir}'
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
    