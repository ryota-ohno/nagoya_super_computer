##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import time
import sys
from tqdm import tqdm
sys.path.append(os.path.join(os.environ['HOME'],'Working/interaction/'))
from make_new_2 import exec_gjf##計算した点のxyzfileを出す
from vdw import vdw_R##同様
from utils import get_E
from utils import get_E_mono_1
from utils import get_E_mono_2
import argparse
import numpy as np
from scipy import signal
import scipy.spatial.distance as distance
import random

def init_process(args):
    # 数理モデル的に自然な定義の元のparams initリスト: not yet
    # 結晶学的に自然なパラメータへ変換: not yet
    auto_dir = args.auto_dir
    monomer_name = args.monomer_name
    order = 5
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)

    def get_init_para_csv(auto_dir,monomer_name):
        init_params_csv = os.path.join(auto_dir, 'step1_init_params.csv')
        
        init_para_list = []
        A1 = 0; A2 = 0
        theta_list = [5,10,15,20,25,30,35,40,45]#list(range(0,95,5)) ##+ [22.5,23.0,23.5,24.0,24.5] + [25.5,26.0,26.5,27.0,27.5]##全体を掃く＋細かいところ
        for theta in tqdm(theta_list):
            a_list = []; b_list = []; S_list = []
            a_clps=vdw_R(A1,A2,theta,0.0,'a',monomer_name)
            b_clps=vdw_R(A1,A2,theta,90.0,'b',monomer_name)
            for theta_ab in range(0,91):
                R_clps=vdw_R(A1,A2,theta,theta_ab,'t',monomer_name)
                a=2*R_clps*np.cos(np.radians(theta_ab))
                b=2*R_clps*np.sin(np.radians(theta_ab))
                if (a_clps > a) or (b_clps > b):
                    continue
                else:
                    a1 = np.round(a,1);b1 = np.round(b,1)
                    a_list.append(a1);b_list.append(b1);S_list.append(a*b)##Sを丸めずに出力
            local_minidx_list = signal.argrelmin(np.array(S_list), order=order)
            if len(local_minidx_list[0])>0:
                for local_minidx in local_minidx_list[0]:
                    init_para_list.append([a_list[local_minidx],b_list[local_minidx],theta,'NotYet'])
            init_para_list.append([a_list[0],b_list[0],theta,'NotYet'])
            init_para_list.append([a_list[-1],b_list[-1],theta,'NotYet'])
            
        df_init_params = pd.DataFrame(np.array(init_para_list),columns = ['a','b','theta','status'])
        df_init_params.to_csv(init_params_csv,index=False)
    
    get_init_para_csv(auto_dir,monomer_name)
    
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path):        
        df_E_init = pd.DataFrame(columns = ['a','b','theta','E','E_p1','E_p2','E_t','machine_type','status','file_name'])##隣接8分子
    else:
        df_E_init = pd.read_csv(auto_csv_path)
        df_E_init = df_E_init[df_E_init['status']!='InProgress']
    df_E_init.to_csv(auto_csv_path,index=False)

    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    df_init['status']='NotYet'
    df_init.to_csv(os.path.join(auto_dir,'step1_init_params.csv'),index=False)

def main_process(args):
    auto_dir = args.auto_dir
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussian'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    auto_csv_path = os.path.join(auto_dir,'step2_twist.csv')
    if not os.path.exists(auto_csv_path):        
        df_E = pd.DataFrame(columns = ['a','b','theta','A1','A2','phi_r','phi_b_r','phi_b_p','E','E_p1','E_p2','E_t1','E_t3','E_m','machine_type','status','file_name'])##いじる
        df_E.to_csv(auto_csv_path,index=False)##step3を二段階でやる場合二段階目ではinitをやらないので念のためmainにも組み込んでおく

    os.chdir(os.path.join(args.auto_dir,'gaussian'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(args.auto_dir,args.monomer_name,args.num_nodes,args.isTest)##argsの中身を取る
        time.sleep(1)

def listen(auto_dir,monomer_name,num_nodes,isTest):##args自体を引数に取るか中身をばらして取るかの違い
    num_init = args.num_init
    fixed_param_keys = ['theta','A1','A2','phi_b_r']
    opt_param_keys = ['a','b','phi_r','phi_b_p']
    df_mono_e=pd.read_csv('~/Working/nagoya_super_computer/BTBT/monomer/{}_e.csv'.format(monomer_name))
    auto_csv = os.path.join(auto_dir,'step2_twist.csv')
    df_E = pd.read_csv(auto_csv)
    df_queue = df_E.loc[df_E['status']=='InProgress',['file_name','phi_r']]
    len_queue = len(df_queue)
    
    for idx,row in zip(df_queue.index,df_queue.values):
        file_name = row[0]
        phi_r=row[1]
        log_filepath = os.path.join(*[auto_dir,'gaussian',file_name])
        if not(os.path.exists(log_filepath)):#logファイルが生成される直前だとまずいので
            continue
        E_list=get_E(log_filepath)
        if len(E_list)!=4:
            continue
        else:
            len_queue-=1
            Ep1=float(E_list[0]);Ep2=float(E_list[1]);Et1=float(E_list[2]);Et3=float(E_list[3])
            if monomer_name == 'mono-C4-BTBT':
                Em=get_E_mono_1(log_filepath)
            elif monomer_name == 'mono-C9-BTBT':
                Em=get_E_mono_2(log_filepath)
            else:
                Em=0
            E = (Ep1+Ep2+Et1+Et3)+Em
            df_E.loc[idx, ['E_t1','E_t3','E_p1','E_p2','E_m','E','status']] = [Et1,Et3,Ep1,Ep2,Em,E,'Done']
            df_E.to_csv(auto_csv,index=False)
            break#2つ同時に計算終わったりしたらまずいので一個で切る
    isAvailable = len_queue < num_nodes 
    if isAvailable:
        dict_matrix = get_params_dict(auto_dir,num_init, fixed_param_keys, opt_param_keys)
        if len(dict_matrix)!=0:#終わりがまだ見えないなら
            for i in range(len(dict_matrix)):
                params_dict=dict_matrix[i]
                alreadyCalculated = check_calc_status(auto_dir,params_dict)
                if not(alreadyCalculated):
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict,'cx':0,'cy':0,'cz':0},isInterlayer=False,isTest=isTest)
                    df_newline = pd.Series({**params_dict,'E':0.,'E_p1':0.,'E_p2':0.,'E_t1':0.,'E_t3':0.,'E_m':0.,'machine_type':'3','status':'InProgress','file_name':file_name})
                    df_E=df_E.append(df_newline,ignore_index=True)
                    df_E.to_csv(auto_csv,index=False)
    
    init_params_csv=os.path.join(auto_dir, 'step2_twist_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params_done = filter_df(df_init_params,{'status':'Done'})
    isOver = True if len(df_init_params_done)==len(df_init_params) else False
    return isOver

def check_calc_status(auto_dir,params_dict):
    df_E= pd.read_csv(os.path.join(auto_dir,'step2_twist.csv'))
    if len(df_E)==0:
        return False
    df_E_filtered = filter_df(df_E, params_dict)
    df_E_filtered = df_E_filtered.reset_index(drop=True)
    try:
        status = get_values_from_df(df_E_filtered,0,'status')
        return status=='Done'
    except KeyError:
        return False

def get_params_dict(auto_dir, num_init,fixed_param_keys,opt_param_keys):
    """
    前提:
        step1_init_params.csvとstep1.csvがauto_dirの下にある
    """
    init_params_csv=os.path.join(auto_dir, 'step2_twist_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_cur = pd.read_csv(os.path.join(auto_dir, 'step2_twist.csv'))
    df_init_params_inprogress = df_init_params[df_init_params['status']=='InProgress']

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
            opt_params_dict={'a':opt_params_matrix[0][0],'b':opt_params_matrix[0][1],'phi_r':opt_params_matrix[0][2],'phi_b_p':opt_params_matrix[0][3]}
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
                opt_params_dict={'a':opt_params_matrix[i][0],'b':opt_params_matrix[i][1],'phi_r':opt_params_matrix[i][2],'phi_b_p':opt_params_matrix[i][3]}
                df_inprogress = filter_df(df_cur, {**fixed_params_dict,**opt_params_dict,'status':'InProgress'})
                if len(df_inprogress)>=1:
                    continue
                else:
                    d={**fixed_params_dict,**opt_params_dict}
                    dict_matrix.append(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    a_init_prev = init_params_dict['a']; b_init_prev = init_params_dict['b']; phi_r_init_prev = init_params_dict['phi_r']; phi_b_init_prev = init_params_dict['phi_b_p']
    A1 = init_params_dict['A1']; A2 = init_params_dict['A2']; theta = init_params_dict['theta']; phi_b_r = init_params_dict['phi_b_r']
    
    while True:
        E_list=[];heri_list=[]
        para_list=[]
        for a in [a_init_prev-0.1,a_init_prev,a_init_prev+0.1]:
            for b in [b_init_prev-0.1,b_init_prev,b_init_prev+0.1]:
                for phi_r in [phi_r_init_prev-1,phi_r_init_prev,phi_r_init_prev+1]:
                    for phi_b in [phi_b_init_prev]:
                        a = np.round(a,1);b = np.round(b,1);phi_r = np.round(phi_r,1)
                        df_val_ab = df_val[
                            (df_val['a']==a)&(df_val['b']==b)&(df_val['theta']==theta)&
                            (df_val['A1']==A1)&(df_val['A2']==A2)&(df_val['phi_r']==phi_r)&(df_val['phi_b_p']==phi_b)&(df_val['phi_b_r']==phi_b_r)&
                            (df_val['status']=='Done')
                                    ]
                        if len(df_val_ab)==0:
                            para_list.append([a,b,phi_r,phi_b])
                            continue
                        heri_list.append([a,b,phi_r,phi_b]);E_list.append(df_val_ab['E'].values[0])
        if len(para_list) != 0:
            return False,para_list
        a_init,b_init,phi_r_init,phi_b_init = heri_list[np.argmin(np.array(E_list))]
        if a_init==a_init_prev and b_init==b_init_prev and phi_r_init==phi_r_init_prev and phi_b_init==phi_b_init_prev:
            return True,[[a_init,b_init,phi_r_init,phi_b_init]]
        else:
            a_init_prev=a_init;b_init_prev=b_init;phi_r_init_prev=phi_r_init;phi_b_init_prev=phi_b_init

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
    
    parser.add_argument('--init',action='store_true')
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    parser.add_argument('--num-init',type=int,help='number of parameters in progress at init_params.csv')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    if args.init:
        print("----initial process----")
        init_process(args)
    
    print("----main process----")
    main_process(args)
    print("----finish process----")
    