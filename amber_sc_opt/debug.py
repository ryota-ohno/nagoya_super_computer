##tetracene層内計算
import os
os.environ['HOME'] ='/home/HasegawaLab'
import pandas as pd
import argparse
import numpy as np

def main_process(args):
    auto_dir = f'/data/group1/z40145w/Working/nagoya_super_computer/amber_sc_opt/BTBT/{args.auto_dir}'
    os.chdir(os.path.join(auto_dir,'amber'))
    isOver = False
    auto_csv_1 = os.path.join(auto_dir,'step1_1.csv');df_E_1 = pd.read_csv(auto_csv_1)
    auto_csv_2 = os.path.join(auto_dir,'step1_2.csv');df_E_2 = pd.read_csv(auto_csv_2)
    auto_csv_3 = os.path.join(auto_dir, 'step1_3.csv');df_E_3 = pd.read_csv(auto_csv_3)
    
    while not(isOver):
        #check
        isOver = listen(auto_dir,df_E_1,df_E_2,df_E_3)##argsの中身を取る
        #time.sleep(0.1)

def listen(auto_dir,df_E_1,df_E_2,df_E_3):##args自体を引数に取るか中身をばらして取るかの違い
    fixed_param_keys = ['theta'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b']
    
    auto_csv = os.path.join(auto_dir,'step1.csv')
    df_E = pd.read_csv(auto_csv)
    df_prg = df_E.loc[df_E['status']=='InProgress',fixed_param_keys+opt_param_keys_1+opt_param_keys_2]
    
    for idx,row in df_prg.iterrows():
        params_dict1_ = row[fixed_param_keys + opt_param_keys_1].to_dict()
        params_dict2_ = row[fixed_param_keys + opt_param_keys_2].to_dict()
        params_dict3_ = row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2].to_dict()
        s1=filter_df(df_E_1, params_dict1_);s2=filter_df(df_E_2, params_dict2_);s3=filter_df(df_E_3, params_dict3_)#['file_name']
        s1=s1[s1['status']=='Done'];s2=s2[s2['status']=='Done'];s3=s3[s3['status']=='Done']
    
        if (len(s1) == 0) or (len(s2) == 0) or (len(s3) == 0):
            continue
        else:
            E1 = s1['E1'].values.tolist()[0]
            E2 = s2['E2'].values.tolist()[0]
            E3 = s3['E3'].values.tolist()[0]
            
            E=2*E1+2*E2+4*E3
            df_E.loc[idx, ['E','E1','E2','E3','status']] = [round(E,4),round(E1,4),round(E2,4),round(E3,4),'Done']
        df_E.to_csv(auto_csv,index=False)
        
    dict_matrix = get_params_dict(auto_dir)##更新分を流す a1/HOME/HASEGAWALABz2まで取得
    if len(dict_matrix)!=0:#終わりがまだ見えないなら
        for i in range(len(dict_matrix)):
            params_dict=dict_matrix[i]#print(params_dict)
            alreadyCalculated = check_calc_status(auto_dir,params_dict)
            if not(alreadyCalculated):
                df_E= pd.read_csv(os.path.join(auto_dir,'step1.csv'))
                df_E_filtered = filter_df(df_E, params_dict)
                if len(df_E_filtered) == 0:
                    df_newline = pd.Series({**params_dict,'E':0.,'E1':0.,'E2':0.,'E3':0.,'status':'InProgress'})
                    df_E_new=pd.concat([df_E,df_newline.to_frame().T],axis=0,ignore_index=True);df_E_new.to_csv(auto_csv,index=False)
    
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

def get_params_dict(auto_dir):
    """
    前提:
        step1_init_params.csvとstep1.csvがauto_dirの下にある
    """
    init_params_csv=os.path.join(auto_dir, 'step1_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_cur = pd.read_csv(os.path.join(auto_dir, 'step1.csv'))
    df_init_params_inprogress = df_init_params[df_init_params['status']=='Done']
    fixed_param_keys = ['theta','A2'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b','z']
    print(df_init_params_inprogress)
    dict_matrix=[]
    for index in df_init_params_inprogress.index:##こちら側はinit_params内のある業に関する探索が終わった際の新しい行での探索を開始するもの ###ここを改良すればよさそう
        df_init_params = pd.read_csv(init_params_csv)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys_1+opt_param_keys_2].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_matrix = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            opt_params_dict={'a':np.round(opt_params_matrix[0][0],1),'b':np.round(opt_params_matrix[0][1],1),'z':np.round(opt_params_matrix[0][2],1)}
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
                opt_params_dict={'a':np.round(opt_params_matrix[i][0],1),'b':np.round(opt_params_matrix[i][1],1),'z':np.round(opt_params_matrix[i][2],1)}
                d={**fixed_params_dict,**opt_params_dict}
                dict_matrix.append(d)
                    #print(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    a_init_prev = init_params_dict['a'];b_init_prev = init_params_dict['b'];z_init_prev = init_params_dict['z']
    while True:
        E_list=[];xyz_list=[]
        para_list=[]
        for a in [a_init_prev-0.1,a_init_prev,a_init_prev+0.1]:
                    for b in [b_init_prev-0.1,b_init_prev,b_init_prev+0.1]:
                        for z in [z_init_prev]:
                            a = np.round(a,1);b = np.round(b,1);z=np.round(z,1)
                            df_val_xyz = df_val[(df_val['a']==a)&(df_val['b']==b)&(df_val['z']==z)&(df_val['status']=='Done')]
                            if len(df_val_xyz)==0:
                                para_list.append([a,b,z])
                                continue
                            xyz_list.append([a,b,z]);E_list.append(df_val_xyz['E'].values[0])
        if len(para_list) != 0:
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
    
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes amber, gaussview and csv')
    args = parser.parse_args()

    print("----main process----")
    main_process(args)
    print("----finish process----")
    