##tetracene層内計算
import os
os.environ['HOME'] ='/home/ohno'
import pandas as pd
import time
import sys
sys.path.append(os.path.join(os.environ['HOME'],'Working/interaction/'))
from make_xyz_new import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse
import numpy as np

def main_process(args):
    auto_dir = args.auto_dir
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    auto_csv_path = os.path.join(auto_dir,'step3.csv')
    if not os.path.exists(auto_csv_path):        
        df_E = pd.DataFrame(columns = ['cx','cy','cz','a','b','theta1','theta2','z','E','Ei0','Eip1','Eip2','Eip3','Eip4','Eit1','Eit2','Eit3','Eit4','Ei0_','Eip1_','Eip2_','Eip3_','Eip4_','Eit1_','Eit2_','Eit3_','Eit4_','status','file_name'])##いじる
        df_E.to_csv(auto_csv_path,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(args.auto_dir,args.monomer_name,args.num_nodes,args.num_init,args.isTest)##argsの中身を取る
        #time.sleep(1)

def listen(auto_dir,monomer_name,num_nodes,num_init,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    auto_csv = os.path.join(auto_dir,'step3.csv')
    df_E = pd.read_csv(auto_csv)
    df_queue = df_E.loc[df_E['status']=='InProgress',['file_name']]
    len_queue = len(df_queue)
    
    for idx,row in zip(df_queue.index,df_queue.values):
        file_name = row[0]
        log_filepath = os.path.join(*[auto_dir,'gaussian',file_name])
        if not(os.path.exists(log_filepath)):#logファイルが生成される直前だとまずいので
            continue
        E_list=get_E(log_filepath)
        if len(E_list)!=18:##get Eの長さは計算した分子の数
            continue
        else:
            len_queue-=1
            Ei0=float(E_list[0]);Eip1=float(E_list[1]);Eip2=float(E_list[2]);Eip3=float(E_list[3]);Eip4=float(E_list[4])
            Ei0_=float(E_list[5]);Eip1_=float(E_list[6]);Eip2_=float(E_list[7]);Eip3_=float(E_list[8]);Eip4_=float(E_list[9])
            Eit1=float(E_list[10]);Eit2=float(E_list[11]);Eit3=float(E_list[12]);Eit4=float(E_list[13])
            Eit1_=float(E_list[14]);Eit2_=float(E_list[15]);Eit3_=float(E_list[16]);Eit4_=float(E_list[17])
            E = ((Ei0 + Eip1 + Eip2 + Eip3 + Eip4 + Eit1 + Eit2 + Eit3 + Eit4) + (Ei0_ + Eip1_ + Eip2_ + Eip3_ + Eip4_ + Eit1_ + Eit2_ + Eit3_ + Eit4_))/2##エネルギーの値も変える
            df_E.loc[idx, ['Ei0','Eip1','Eip2','Eip3','Eip4','Ei0_','Eip1_','Eip2_','Eip3_','Eip4_','Eit1','Eit2','Eit3','Eit4','Eit1_','Eit2_','Eit3_','Eit4_','E','status']] = [Ei0,Eip1,Eip2,Eip3,Eip4,Ei0_,Eip1_,Eip2_,Eip3_,Eip4_,Eit1,Eit2,Eit3,Eit4,Eit1_,Eit2_,Eit3_,Eit4_,E,'Done']
            df_E.to_csv(auto_csv,index=False)
            #time.sleep(2)
            #break#2つ同時に計算終わったりしたらまずいので一個で切る
    isAvailable = len_queue < num_nodes 
    if isAvailable:
        dict_matrix = get_params_dict(auto_dir,num_init)
        if len(dict_matrix)!=0:#終わりがまだ見えないなら
            for i in range(len(dict_matrix)):
                params_dict=dict_matrix[i]
                alreadyCalculated = check_calc_status(auto_dir,params_dict)
                if not(alreadyCalculated):
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict},isInterlayer=False,isTest=isTest)
                    df_newline = pd.Series({**params_dict,'E':0.,'Ei0':0.,'Eip1':0.,'Eip2':0.,'Eip3':0.,'Eip4':0.,'Ei0_':0.,'Eip1_':0.,'Eip2_':0.,'Eip3_':0.,'Eip4_':0.,'Eit1':0.,'Eit2':0.,'Eit3':0.,'Eit4':0.,'machine_type':'3','status':'InProgress','file_name':file_name})
                    df_E=df_E.append(df_newline,ignore_index=True)
                    df_E.to_csv(auto_csv,index=False)
                    #time.sleep(2)
    
    init_params_csv=os.path.join(auto_dir, 'step3_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params_done = filter_df(df_init_params,{'status':'Done'})
    isOver = True if len(df_init_params_done)==len(df_init_params) else False
    return isOver

def check_calc_status(auto_dir,params_dict):
    df_E= pd.read_csv(os.path.join(auto_dir,'step3.csv'))
    if len(df_E)==0:
        return False
    df_E_filtered = filter_df(df_E, params_dict)
    df_E_filtered = df_E_filtered.reset_index(drop=True)
    try:
        status = get_values_from_df(df_E_filtered,0,'status')
        return status=='Done'
    except KeyError:
        return False

def get_params_dict(auto_dir, num_init):
    """
    前提:
        step3_init_params.csvとstep3.csvがauto_dirの下にある
    """
    init_params_csv=os.path.join(auto_dir, 'step3_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_cur = pd.read_csv(os.path.join(auto_dir, 'step3.csv'))
    df_init_params_inprogress = df_init_params[df_init_params['status']=='InProgress']
    fixed_param_keys = ['a','b','theta1','theta2','z']
    opt_param_keys = ['cx','cy','cz']

    #最初の立ち上がり時
    if len(df_init_params_inprogress) < num_init:
        df_init_params_notyet = df_init_params[df_init_params['status']=='NotYet']
        for index in df_init_params_notyet.index:
            df_init_params = update_value_in_df(df_init_params,index,'status','InProgress')
            df_init_params.to_csv(init_params_csv,index=False)
            params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
            return [params_dict]
    dict_matrix=[]
    for index in df_init_params_inprogress.index:##こちら側はinit_params内のある業に関する探索が終わった際の新しい行での探索を開始するもの ###ここを改良すればよさそう
        df_init_params = pd.read_csv(init_params_csv)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_matrix = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            opt_params_dict={'cx':opt_params_matrix[0][0],'cy':opt_params_matrix[0][1],'cz':opt_params_matrix[0][2]}
            # df_init_paramsのstatusをupdate
            df_init_params = update_value_in_df(df_init_params,index,'status','Done')
            if np.max(df_init_params.index) < index+1:##もうこれ以上は新しい計算は進まない
                status = 'Done'
            else:
                status = get_values_from_df(df_init_params,index+1,'status')
            df_init_params.to_csv(init_params_csv,index=False)
            
            if status=='NotYet':##計算が始まっていないものがあったらこの時点で開始する　ここでダメでもまた直にlistenでgrt_params_dictまでいけば新しいのが始まる            
                opt_params_dict = get_values_from_df(df_init_params,index+1,opt_param_keys)
                df_init_params = update_value_in_df(df_init_params,index+1,'status','InProgress')
                df_init_params.to_csv(init_params_csv,index=False)
                dict_matrix.append({**fixed_params_dict,**opt_params_dict})
            else:
                continue

        else:
            for i in range(len(opt_params_matrix)):
                opt_params_dict={'cx':opt_params_matrix[i][0],'cy':opt_params_matrix[i][1],'cz':opt_params_matrix[i][2]}
                df_inprogress = filter_df(df_cur, {**fixed_params_dict,**opt_params_dict,'status':'InProgress'})
                if len(df_inprogress)>=1:
                    continue
                else:
                    d={**fixed_params_dict,**opt_params_dict}
                    dict_matrix.append(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    cx_init_prev = init_params_dict['cx']; cy_init_prev = init_params_dict['cy']; cz_init_prev = init_params_dict['cz']
    while True:
        E_list=[];cxyz_list=[]
        para_list=[]
        for cx in [cx_init_prev]:#-0.1,cx_init_prev,cx_init_prev+0.1
            for cy in [cy_init_prev]:#-0.1,cy_init_prev,cy_init_prev+0.1
                for cz in [cz_init_prev]:#-0.1,cz_init_prev,cz_init_prev+0.1
                    cx = np.round(cx,1);cy = np.round(cy,1);cz = np.round(cz,1)
                    df_val_cxyz = df_val[
                        (df_val['cx']==cx)&(df_val['cy']==cy)&(df_val['cz']==cz)&(df_val['status']=='Done')]
                    if len(df_val_cxyz)==0:
                        para_list.append([cx,cy,cz])
                        continue
                    cxyz_list.append([cx,cy,cz]);E_list.append(df_val_cxyz['E'].values[0])
        if len(para_list) != 0:
            return False,para_list
        cx_init,cy_init,cz_init = cxyz_list[np.argmin(np.array(E_list))]
        if cx_init==cx_init_prev and cy_init==cy_init_prev and cz_init==cz_init_prev:
            return True,[[cx_init,cy_init,cz_init]]
        else:
            cx_init_prev=cx_init;cy_init_prev=cy_init;cz_init_prev=cz_init

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
    parser.add_argument('--num-init',type=int,help='num nodes')
    args = parser.parse_args()
    print("----main process----")
    main_process(args)
    print("----finish process----")
    