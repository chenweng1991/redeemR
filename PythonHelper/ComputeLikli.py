"""
Utilities to compute the liklihood for assigning variants to nodes
Work together with R function Add_AssignVariant in scMitoTracing

By Chen, Created 2022-7-8
"""
import numpy as np
import pandas as pd
import time
import sys
from multiprocessing import Pool

path=sys.argv[1]  # The path for intermediate csv files created by R
cores=int(sys.argv[2]) # Number of cores to be used

#path="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Old1Old2/Shell_Old1Old2.basics/MakeTree/tmp"

## Load in the maties from R
mtr_bi_pd=pd.read_csv(path+"/mtr.bi.csv",index_col=0)
mtr_bi=pd.DataFrame.to_numpy(pd.read_csv(path+"/mtr.bi.csv",index_col=0))
df_profile_mtx_t=pd.DataFrame.to_numpy(pd.read_csv(path+"/df_profile_mtx.t.csv",index_col=0))
NonZero_logV=pd.DataFrame.to_numpy(pd.read_csv(path+"/NonZero.logV.csv",index_col=0))
Zero_logV=pd.DataFrame.to_numpy(pd.read_csv(path+"/Zero.logV.csv",index_col=0))
NonZero_logP=pd.DataFrame.to_numpy(pd.read_csv(path+"/NonZero.logP.csv",index_col=0))
Zero_logP=pd.DataFrame.to_numpy(pd.read_csv(path+"/Zero.logP.csv",index_col=0))
print("Intermediate files in ComputeLikli.py")
# mtr_bi_pd.index[1]  ## mtr_bi_pd is variant(each row) VS cell(each column)

## Get the number of variants, or the number of rows in mtr_bi
NumVariant=mtr_bi.shape[0] 

def ComputeLoglik(i):
    x=df_profile_mtx_t+mtr_bi[i,:][:,None]
    y=df_profile_mtx_t-mtr_bi[i,:][:,None]
    z=mtr_bi[i,:][:,None]-df_profile_mtx_t
    ## Make indicator matries for 1-1 (Inside the clade, and is mutation) and 0-0 (outside of clade, and is not mutation) 
    x_11=np.where(x==2,1,0)
    x_00=np.where(x==0,1,0)
    ## Make indicator matries for 1-0 (Inside the clade, and is NOT mutation) and 0-1 (outside of clade, and is mutation) 
    x_10=np.where(y==-1,0,y)
    x_01=np.where(z==-1,0,z)
    ## Compute the Loglikihood across all nodes
    Loglik_11=x_11*NonZero_logV[i,][:,None]
    Loglik_10=x_10*Zero_logV[i,][:,None]
    Loglik_01=x_01*NonZero_logP
    Loglik_00=x_00*Zero_logP
    Loglik_i=(Loglik_11+Loglik_10+Loglik_01+Loglik_00).sum(axis=0)
    #print(str(i)+"\n")
    return(Loglik_i.tolist())

##
start=time.time()
print("In ComputeLikli.py, start pooling and parallele computing")
p = Pool(cores)
Res=p.map(ComputeLoglik, range(0,NumVariant))
p.close()
end=time.time()
print("ComputeLikli.py complete, running time is "+str(end-start))

Res_pd=pd.DataFrame(Res)
Res_pd.index=mtr_bi_pd.index
Res_pd.to_csv(path+"/likliRes.csv")
