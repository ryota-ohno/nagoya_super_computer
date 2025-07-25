##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import time
from make_8 import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse
import numpy as np

def main_process(args):
    auto_dir = args.auto_dir
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path):        
        df_E = pd.DataFrame(columns = ['theta','A2','a','b','z','E','E1','E2','E3','E4','status'])##いじる
        df_E.to_csv(auto_csv_path,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    auto_csv_path1 = os.path.join(auto_dir,'step1_1.csv')
    if not os.path.exists(auto_csv_path1):        
        df_E_1 = pd.DataFrame(columns = ['theta','A2','a','E1','machine_type','status','file_name'])##いじる
        df_E_1.to_csv(auto_csv_path1,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    auto_csv_path2 = os.path.join(auto_dir,'step1_2.csv')
    if not os.path.exists(auto_csv_path2):        
        df_E_2 = pd.DataFrame(columns = ['theta','A2','b','z','E2','machine_type','status','file_name'])##いじる
        df_E_2.to_csv(auto_csv_path2,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    auto_csv_path3 = os.path.join(auto_dir,'step1_3.csv')
    if not os.path.exists(auto_csv_path3):        
        df_E_3 = pd.DataFrame(columns = ['theta','A2','a','b','z','E3','E4','machine_type','status','file_name'])##いじる
        df_E_3.to_csv(auto_csv_path3,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(args.auto_dir,args.monomer_name,args.num_nodes,args.max_nodes,args.isTest)##argsの中身を取る
        #time.sleep(1)

def listen(auto_dir,monomer_name,num_nodes,max_nodes,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    fixed_param_keys = ['theta','A2'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b','z']
    
    auto_csv_1 = os.path.join(auto_dir,'step1_1.csv');df_E_1 = pd.read_csv(auto_csv_1)
    df_prg_1 = df_E_1.loc[df_E_1['status']=='InProgress',fixed_param_keys+opt_param_keys_1+['machine_type','file_name']]
    len_prg_1=len(df_prg_1)
    for idx, row in df_prg_1.iterrows():
        params_dict1_ = row[fixed_param_keys + opt_param_keys_1 + ['file_name']].to_dict()
        file_name1=params_dict1_['file_name']##辞書をつくってそこにopt_1とopt_2でファイル名作成
        log_filepath1 = os.path.join(*[auto_dir,'gaussian',file_name1])
        if not(os.path.exists(log_filepath1)):#logファイルが生成される直前だとまずいので
            continue
        E_list1=get_E(log_filepath1)
        if len(E_list1)!=1 :##get Eの長さは計算した分子の数
            continue
        else:
            len_prg_1-=1
            E1=float(E_list1[0])##8分子に向けてep1,ep2作成　ep1:b ep2:a
            df_E_1.loc[idx, ['E1','status']] = [E1,'Done']
            df_E_1.to_csv(auto_csv_1,index=False)
            #time.sleep(1)
            break#2つ同時に計算終わったりしたらまずいので一個で切る
    
    
    auto_csv_2 = os.path.join(auto_dir,'step1_2.csv')
    df_E_2 = pd.read_csv(auto_csv_2)
    df_prg_2 = df_E_2.loc[df_E_2['status']=='InProgress', fixed_param_keys+opt_param_keys_2+['machine_type','file_name']]
    len_prg_2 = len(df_prg_2)

    for idx, row in df_prg_2.iterrows():
        params_dict2_ = row[fixed_param_keys + opt_param_keys_2 + ['file_name']].to_dict()
        file_name2=params_dict2_['file_name']##辞書をつくってそこにopt_1とopt_2でファイル名作成
        log_filepath2 = os.path.join(*[auto_dir, 'gaussian', file_name2])
        if not(os.path.exists(log_filepath2)):
            continue
        E_list2 = get_E(log_filepath2)
        if len(E_list2) != 1:
            continue
        else:
            len_prg_2 -= 1
            E2 = float(E_list2[0])  # Updated to E2
            df_E_2.loc[idx, ['E2', 'status']] = [E2, 'Done']
            df_E_2.to_csv(auto_csv_2, index=False)  # Updated to auto_csv_2
            #time.sleep(1)
            break  # Break after one iteration

    
    auto_csv_3 = os.path.join(auto_dir, 'step1_3.csv')
    df_E_3 = pd.read_csv(auto_csv_3)
    df_prg_3 = df_E_3.loc[df_E_3['status'] == 'InProgress', fixed_param_keys + opt_param_keys_1 + opt_param_keys_2 + ['machine_type','file_name']]
    len_prg_3 = len(df_prg_3)

    for idx, row in df_prg_3.iterrows():
        params_dict3_ = row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2 + ['file_name']].to_dict()
        file_name3=params_dict3_['file_name']##辞書をつくってそこにopt_1とopt_2でファイル名作成
        log_filepath3 = os.path.join(*[auto_dir, 'gaussian', file_name3])
        if not (os.path.exists(log_filepath3)):
            continue
        E_list3 = get_E(log_filepath3)
        if len(E_list3) != 1:
            continue
        else:
            len_prg_3 -= 1
            E3 = float(E_list3[0]);E4 = float(E_list3[0])  # Updated to E3
            df_E_3.loc[idx, ['E3', 'E4', 'status']] = [E3, E4, 'Done']
            df_E_3.to_csv(auto_csv_3, index=False)  # Updated to auto_csv_3
            break  # Break after one iteration


    auto_csv = os.path.join(auto_dir,'step1.csv')
    df_E = pd.read_csv(auto_csv)
    df_prg = df_E.loc[df_E['status']=='InProgress',fixed_param_keys+opt_param_keys_1+opt_param_keys_2]
    
    for idx,row in df_prg.iterrows():
        params_dict1_ = row[fixed_param_keys + opt_param_keys_1].to_dict()
        params_dict2_ = row[fixed_param_keys + opt_param_keys_2].to_dict()
        params_dict3_ = row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2].to_dict()
        s1=filter_df(df_E_1, params_dict1_)['file_name'];s2=filter_df(df_E_2, params_dict2_)['file_name'];s3=filter_df(df_E_3, params_dict3_)['file_name']
        if (len(s1) == 0) or (len(s2) == 0) or (len(s3) == 0):
            continue
        file_name1 = str(s1.values[0]);file_name2 = str(s2.values[0]);file_name3 = str(s3.values[0])
        #print(file_name1)
        log_filepath1 = os.path.join(*[auto_dir,'gaussian',file_name1])
        #print(log_filepath1)
        if not(os.path.exists(log_filepath1)):#logファイルが生成される直前だとまずいので
            continue
        E_list1=get_E(log_filepath1)
        
        log_filepath2 = os.path.join(*[auto_dir,'gaussian',file_name2])
        if not(os.path.exists(log_filepath2)):#logファイルが生成される直前だとまずいので
            continue
        E_list2=get_E(log_filepath2)
        
        log_filepath3 = os.path.join(*[auto_dir,'gaussian',file_name3])
        if not(os.path.exists(log_filepath3)):#logファイルが生成される直前だとまずいので
            continue
        E_list3=get_E(log_filepath3)
        
        if (len(E_list1)!=1) or (len(E_list2)!=1) or (len(E_list3)!=1):##get Eの長さは計算した分子の数
            continue
        else:
            E1=float(E_list1[0]);E2=float(E_list2[0]);E3=float(E_list3[0]);E4=float(E_list3[0])##8分子に向けてep1,ep2作成　ep1:b ep2:a
            E=E1+E2+E3+E4
            df_E.loc[idx, ['E','E1','E2','E3','E4','status']] = [E,E1,E2,E3,E4,'Done']
            df_E.to_csv(auto_csv,index=False)
            break#2つ同時に計算終わったりしたらまずいので一個で切る
    
    df_qw_1 = df_E_1[df_E_1['status'] == 'qw'];df_qw_2 = df_E_2[df_E_2['status'] == 'qw'];df_qw_3 = df_E_3[df_E_3['status'] == 'qw']
    len_queue = len_prg_1 + len_prg_2 + len_prg_3
    len_qw_1 = len(df_qw_1);len_qw_2 = len(df_qw_2);len_qw_3 = len(df_qw_3)
    margin = max_nodes - len_queue;machine_type=3

    if len_qw_1 > 0 and margin > 0:# 進行中ジョブのマシンタイプをカウント
        for index, row in df_qw_1.iterrows():
            if margin == 0:
                break
            params_dict = row[fixed_param_keys + opt_param_keys_1].to_dict()# パラメータの辞書を作成
            file_name = exec_gjf(auto_dir, monomer_name, {**params_dict}, machine_type, structure_type=1, isTest=isTest)# ジョブの実行 structure type
            len_queue += 1;margin -= 1
            df_E_1.at[index, 'machine_type'] = machine_type
            df_E_1.at[index, 'status'] = 'InProgress'
            df_E_1.at[index, 'file_name'] = file_name
        df_E_1.to_csv(auto_csv_1, index=False)
        #time.sleep(2)#
    
    if len_qw_2 > 0 and margin > 0:# 進行中ジョブのマシンタイプをカウント
        for index, row in df_qw_2.iterrows():
            if margin == 0:
                break
            params_dict = row[fixed_param_keys + opt_param_keys_2].to_dict()# パラメータの辞書を作成
            
            file_name = exec_gjf(auto_dir, monomer_name, {**params_dict}, machine_type, structure_type=2,isTest=isTest)# ジョブの実行
            len_queue += 1;margin -= 1
            df_E_2.at[index, 'machine_type'] = machine_type
            df_E_2.at[index, 'status'] = 'InProgress'
            df_E_2.at[index, 'file_name'] = file_name
        df_E_2.to_csv(auto_csv_2, index=False)
        #time.sleep(2)# データフレームをCSVに保存
    
    if len_qw_3 > 0 and margin > 0:# 進行中ジョブのマシンタイプをカウント
        for index, row in df_qw_3.iterrows():
            if margin == 0:
                break
            params_dict = row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2].to_dict()# パラメータの辞書を作成
            
            file_name = exec_gjf(auto_dir, monomer_name, {**params_dict}, machine_type, structure_type=3,isTest=isTest)# ジョブの実行
            len_queue += 1;margin -= 1
            df_E_3.at[index, 'machine_type'] = machine_type
            df_E_3.at[index, 'status'] = 'InProgress'
            df_E_3.at[index, 'file_name'] = file_name
        df_E_3.to_csv(auto_csv_3, index=False)
        #time.sleep(2)# データフレームをCSVに保存
    
    dict_matrix = get_params_dict(auto_dir,num_nodes)##更新分を流す a1~z2まで取得
    if len(dict_matrix)!=0:#終わりがまだ見えないなら
        for i in range(len(dict_matrix)):
            params_dict=dict_matrix[i]#print(params_dict)
            params_dict1 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys_1)}
            params_dict2 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys_2)}
            params_dict3 = params_dict
            alreadyCalculated = check_calc_status(auto_dir,params_dict)
            if not(alreadyCalculated):##ここで各点について1~3を実行しつつ余ったものをqwにぶち込む
                df_E= pd.read_csv(os.path.join(auto_dir,'step1.csv'))
                df_E_filtered = filter_df(df_E, params_dict)
                if len(df_E_filtered) == 0:
                    df_newline = pd.Series({**params_dict,'E':0.,'E1':0.,'E2':0.,'E3':0.,'E4':0.,'status':'InProgress'})
                    df_E_new=pd.concat([df_E,df_newline.to_frame().T],axis=0,ignore_index=True);df_E_new.to_csv(auto_csv,index=False)
                
                ## 1の実行　##
                auto_csv_1 = os.path.join(auto_dir,'step1_1.csv');df_E_1 = pd.read_csv(auto_csv_1)
                df_sub_1 = filter_df(df_E_1, params_dict1)
                #df_done_1_ = filter_df(df_E_1, {**params_dict1,'status':'Done'});df_inpr_1_ = filter_df(df_E_1, {**params_dict1,'status':'InProgress'});df_qw_1_ = filter_df(df_E_1, {**params_dict1,'status':'qw'})
                if len(df_sub_1) == 0:#(len(df_done_1_)>=1) or(len(df_inpr_1_)>=1) or (len(df_qw_1_)>=1):
                    print('debug1')
                    isAvailable = len_queue < max_nodes
                    if isAvailable:
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict1}, machine_type, structure_type=1,isTest=isTest);len_queue +=1
                        df_newline_1 = pd.Series({**params_dict1,'E1':0.,'machine_type':machine_type,'status':'InProgress','file_name':file_name})
                        df_E_new_1=pd.concat([df_E_1,df_newline_1.to_frame().T],axis=0,ignore_index=True);df_E_new_1.to_csv(auto_csv_1,index=False)
                        #time.sleep(1)
                    else:
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict1}, machine_type, structure_type=1,isTest=True)
                        df_newline_1 = pd.Series({**params_dict1,'E1':0.,'machine_type':3,'status':'qw','file_name':file_name})
                        df_E_new_1=pd.concat([df_E_1,df_newline_1.to_frame().T],axis=0,ignore_index=True);df_E_new_1.to_csv(auto_csv_1,index=False)
                        #time.sleep(1)

                ## 2の実行　##
                auto_csv_2 = os.path.join(auto_dir,'step1_2.csv');df_E_2 = pd.read_csv(auto_csv_2)
                df_sub_2 = filter_df(df_E_2, params_dict2)
                #df_done_1_ = filter_df(df_E_1, {**params_dict1,'status':'Done'});df_inpr_1_ = filter_df(df_E_1, {**params_dict1,'status':'InProgress'});df_qw_1_ = filter_df(df_E_1, {**params_dict1,'status':'qw'})
                if len(df_sub_2) == 0:#(len(df_done_1_)>=1) or(len(df_inpr_1_)>=1) or (len(df_qw_1_)>=1):
                    print('debug2')
                    isAvailable = len_queue < max_nodes
                    if isAvailable:
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict2}, machine_type, structure_type=2,isTest=isTest);len_queue += 1
                        df_newline_2 = pd.Series({**params_dict2,'E2':0.,'machine_type':machine_type,'status':'InProgress','file_name':file_name})
                        df_E_new_2=pd.concat([df_E_2,df_newline_2.to_frame().T],axis=0,ignore_index=True);df_E_new_2.to_csv(auto_csv_2,index=False)
                    else:
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict2}, machine_type, structure_type=2,isTest=True)
                        df_newline_2 = pd.Series({**params_dict2,'E2':0.,'machine_type':3,'status':'qw','file_name':file_name})
                        df_E_new_2=pd.concat([df_E_2,df_newline_2.to_frame().T],axis=0,ignore_index=True);df_E_new_2.to_csv(auto_csv_2,index=False)

                ## 3の実行　##
                auto_csv_3 = os.path.join(auto_dir,'step1_3.csv');df_E_3 = pd.read_csv(auto_csv_3)
                df_sub_3 = filter_df(df_E_3, params_dict3)
                #df_done_1_ = filter_df(df_E_1, {**params_dict1,'status':'Done'});df_inpr_1_ = filter_df(df_E_1, {**params_dict1,'status':'InProgress'});df_qw_1_ = filter_df(df_E_1, {**params_dict1,'status':'qw'})
                if len(df_sub_3) == 0:#(len(df_done_1_)>=1) or(len(df_inpr_1_)>=1) or (len(df_qw_1_)>=1):
                    print('debug3')
                    isAvailable = len_queue < max_nodes
                    if isAvailable:
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict3}, machine_type, structure_type=3,isTest=isTest);len_queue +=1
                        df_newline_3 = pd.Series({**params_dict3,'E3':0.,'E4':0.,'machine_type':machine_type,'status':'InProgress','file_name':file_name})
                        df_E_new_3=pd.concat([df_E_3,df_newline_3.to_frame().T],axis=0,ignore_index=True);df_E_new_3.to_csv(auto_csv_3,index=False)
                    else:
                        file_name = exec_gjf(auto_dir, monomer_name, {**params_dict3}, machine_type, structure_type=3,isTest=True)
                        df_newline_3 = pd.Series({**params_dict3,'E3':0.,'E4':0.,'machine_type':3,'status':'qw','file_name':file_name})
                        df_E_new_3=pd.concat([df_E_3,df_newline_3.to_frame().T],axis=0,ignore_index=True);df_E_new_3.to_csv(auto_csv_3,index=False)
    
    init_params_csv=os.path.join(auto_dir, 'step1_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params_done = filter_df(df_init_params,{'status':'Done'})
    isOver = True if len(df_init_params_done)==len(df_init_params) else False
    return isOver

def check_calc_status(auto_dir,params_dict):
    df_E= pd.read_csv(os.path.join(auto_dir,'step1.csv'))
    if len(df_E)==0:
        return False
    df_E_filtered = filter_df(df_E, params_dict)
    df_E_filtered = df_E_filtered.reset_index(drop=True)
    try:
        status = get_values_from_df(df_E_filtered,0,'status')
        return status=='Done'
    except KeyError:
        return False

def get_params_dict(auto_dir, num_nodes):
    """
    前提:
        step1_init_params.csvとstep1.csvがauto_dirの下にある
    """
    init_params_csv=os.path.join(auto_dir, 'step1_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_cur = pd.read_csv(os.path.join(auto_dir, 'step1.csv'))
    df_init_params_inprogress = df_init_params[df_init_params['status']=='InProgress']
    fixed_param_keys = ['theta','A2'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b','z']
    
    #最初の立ち上がり時
    if len(df_init_params_inprogress) < num_nodes:
        #print(1)
        df_init_params_notyet = df_init_params[df_init_params['status']=='NotYet']
        for index in df_init_params_notyet.index:
            df_init_params = update_value_in_df(df_init_params,index,'status','InProgress')
            df_init_params.to_csv(init_params_csv,index=False)
            params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys_1+opt_param_keys_2].to_dict()
            return [params_dict]
    dict_matrix=[]
    for index in df_init_params_inprogress.index:##こちら側はinit_params内のある業に関する探索が終わった際の新しい行での探索を開始するもの ###ここを改良すればよさそう
        df_init_params = pd.read_csv(init_params_csv)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys_1+opt_param_keys_2].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_matrix = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            opt_params_dict={'a':np.round(opt_params_matrix[0][0],1),
                            'b':np.round(opt_params_matrix[0][1],1),'z':np.round(opt_params_matrix[0][2],1)
                            }
            # df_init_paramsのstatusをupdate
            df_init_params = update_value_in_df(df_init_params,index,'status','Done')
            if np.max(df_init_params.index) < index+1:##もうこれ以上は新しい計算は進まない
                status = 'Done'
            else:
                status = get_values_from_df(df_init_params,index+1,'status')
            df_init_params.to_csv(init_params_csv,index=False)
            
            if status=='NotYet':##計算が始まっていないものがあったらこの時点で開始する　ここでダメでもまた直にlistenでgrt_params_dictまでいけば新しいのが始まる            
                opt_params_dict = get_values_from_df(df_init_params,index+1,fixed_param_keys+opt_param_keys_1+opt_param_keys_2)
                df_init_params = update_value_in_df(df_init_params,index+1,'status','InProgress')
                df_init_params.to_csv(init_params_csv,index=False)
                dict_matrix.append({**fixed_params_dict,**opt_params_dict})
            else:
                continue

        else:
            for i in range(len(opt_params_matrix)):
                opt_params_dict={'a':np.round(opt_params_matrix[i][0],1),
                            'b':np.round(opt_params_matrix[i][1],1),'z':np.round(opt_params_matrix[i][2],1)
                            }
                d={**fixed_params_dict,**opt_params_dict}
                dict_matrix.append(d)
                    #print(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    a_init_prev = init_params_dict['a']
    b_init_prev = init_params_dict['b']; z_init_prev = init_params_dict['z']
    while True:
        E_list=[];xyz_list=[]
        para_list=[]
        for a in [a_init_prev-0.1,a_init_prev,a_init_prev+0.1]:#]
                    for b in [b_init_prev-0.1,b_init_prev,b_init_prev+0.1]:
                        for z in [z_init_prev]:
                            a = np.round(a,1)
                            b = np.round(b,1);z = np.round(z,1)
                            df_val_xyz = df_val[(df_val['a']==a)&
                                                (df_val['b']==b)&(df_val['z']==z)&
                                                (df_val['status']=='Done')]
                            if len(df_val_xyz)==0:
                                para_list.append([a,b,z])
                                continue
                            xyz_list.append([a,b,z]);E_list.append(df_val_xyz['E'].values[0])
        if len(para_list) != 0:
            print(para_list)
            return False,para_list
        a_init,b_init,z_init = xyz_list[np.argmin(np.array(E_list))]
        if a_init==a_init_prev and b_init==b_init_prev and z_init==z_init_prev:
            return True,[[a_init,b_init,z_init]]
        else:
            a_init_prev=a_init;b_init_prev=b_init;z_init_prev=z_init

def get_values_from_df(df,index,key):
    return df.loc[index,key]

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

def filter_df(df, dict_filter):
    for k, v in dict_filter.items():
        if type(v)==str:
            df=df[df[k]==v]
        else:
            df=df[df[k]==v]
    df_filtered=df
    return df_filtered

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    parser.add_argument('--max-nodes',type=int,help='max nodes')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    main_process(args)
    print("----finish process----")
    