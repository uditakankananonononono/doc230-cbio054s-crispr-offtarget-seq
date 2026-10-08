import numpy as np,pandas as pd,json,torch,torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score
torch.set_num_threads(2); rng=np.random.default_rng(0); torch.manual_seed(0)
D='/tmp/crispr/'; SYM={c:i for i,c in enumerate('ACGTN')}
def enc(s,L=24):
    s=s.upper(); s=('-'*(L-len(s)))+s
    return np.array([SYM.get(c,5) for c in s],dtype=np.int8)
def load(f,sg,off,lab,grp=None):
    d=pd.read_csv(D+f); S=np.stack([enc(x) for x in d[sg]]); O=np.stack([enc(x) for x in d[off]])
    return S,O,(d[lab].values>0).astype(int),(d[grp].values if grp else None)
S,O,y,g=load('I_1_CIRCLE_seq_10gRNA_wholeDataset.csv','sgRNA_seq','off_seq','label','sgRNA_type')
S[:,0]=5;O[:,0]=5  # first char is a placeholder
Sk,Ok,yk,_=load('II_5_Kleinstiver_5gRNA_wholeDataset.csv','sgRNA_seq','off_seq','label')
Sl,Ol,yl,_=load('II_6_Listgarten_22gRNA_wholeDataset.csv','sgRNA_seq','off_seq','label')
def mm(S,O): return ((S!=O)&(S!=4)&(O!=4)).astype(np.float32)  # mismatch or gap-vs-base
def onehot(S,O):
    X=np.zeros((len(S),12,24),dtype=np.float32); ar=np.arange(24)
    for i,A in enumerate((S,O)):
        for j in range(6): X[:,i*6+j,:]=(A==j)
    return X
class Net(nn.Module):
    def __init__(s):
        super().__init__(); s.a=nn.Conv1d(12,64,5,padding=2); s.b=nn.Conv1d(64,64,3,padding=1); s.o=nn.Linear(128,1)
    def forward(s,x):
        h=torch.relu(s.b(torch.relu(s.a(x)))); return s.o(torch.cat([h.max(2)[0],h.mean(2)],1)).squeeze(-1)
def train_idx(idx,yy,ratio,r):
    p=idx[yy[idx]==1]; n=idx[yy[idx]==0]; n=r.choice(n,min(len(n),ratio*len(p)),replace=False); return np.concatenate([p,n])
def fitM(S_,O_,y_,seed=0,shuffle=False):
    r=np.random.default_rng(seed); torch.manual_seed(seed); X=torch.tensor(onehot(S_,O_)); yt=y_.copy()
    if shuffle: yt=r.permutation(yt)
    yt=torch.tensor(yt,dtype=torch.float32); net=Net(); opt=torch.optim.Adam(net.parameters(),2e-3); lf=nn.BCEWithLogitsLoss()
    for ep in range(3):
        pm=torch.randperm(len(X))
        for i in range(0,len(X),128):
            b=pm[i:i+128]; opt.zero_grad(); lf(net(X[b]),yt[b]).backward(); opt.step()
    net.eval()
    def f(S2,O2):
        out=[]
        with torch.no_grad():
            for i in range(0,len(S2),20000): out.append(net(torch.tensor(onehot(S2[i:i+20000],O2[i:i+20000]))).numpy())
        return np.concatenate(out)
    return f
def met(yt,s): return dict(auprc=round(float(average_precision_score(yt,s)),4),auroc=round(float(roc_auc_score(yt,s)),4),base=round(float(yt.mean()),5))
def evaluate(tr,te,S_,O_,y_,ratio=10,seed=0,posmask=None,shuffle=False,only=('B1','B2','M')):
    r=np.random.default_rng(seed); ti=train_idx(tr,y_,ratio,r); res={}
    if 'B1' in only: res['B1']=met(y_[te],-mm(S_[te],O_[te]).sum(1))
    if 'B2' in only:
        F=lambda a,b:(mm(a,b) if posmask is None else mm(a,b)*posmask)
        lr=LogisticRegression(C=1,max_iter=300).fit(F(S_[ti],O_[ti]),y_[ti]); res['B2']=met(y_[te],lr.decision_function(F(S_[te],O_[te])))
    if 'M' in only:
        f=fitM(S_[ti],O_[ti],y_[ti],seed,shuffle); res['M']=met(y_[te],f(S_[te],O_[te]))
    return res
guides=sorted(set(g)); r0=np.random.default_rng(1); r0.shuffle(guides); folds=[guides[i::5] for i in range(5)]
def agg(rs): return {k:{m:round(float(np.mean([r[k][m] for r in rs])),4) for m in rs[0][k]} for k in rs[0]}
R={'d_stats':dict(n=len(y),pos=int(y.sum()),guides=len(guides),kleinstiver=[len(yk),int(yk.sum())],listgarten=[len(yl),int(yl.sum())])}
cold=[]
for k,fg in enumerate(folds):
    te=np.where(np.isin(g,fg))[0]; tr=np.where(~np.isin(g,fg))[0]; cold.append(evaluate(tr,te,S,O,y,seed=k)); print('cold fold',k,cold[-1],flush=True)
R['S1_cold_guide']=agg(cold); R['S1_folds']=cold
bb=max(R['S1_cold_guide'][k]['auprc'] for k in ('B1','B2')); R['S1_gate_M_beats_best_baseline_by_0.02']=bool(R['S1_cold_guide']['M']['auprc']>=bb+.02)
perm=np.random.default_rng(2).permutation(len(y)); rnd=[]
for k in range(5):
    te=perm[k::5]; tr=np.setdiff1d(perm,te); rnd.append(evaluate(tr,te,S,O,y,seed=k)); print('rand fold',k,rnd[-1],flush=True)
R['S2_random_split']=agg(rnd); R['S2_gate_gap_le_0.05']=bool(R['S2_random_split']['M']['auprc']-R['S1_cold_guide']['M']['auprc']<=.05)
pm_seed=np.zeros(24,dtype=np.float32); pm_seed[14:23]=1; pm_dist=np.zeros(24,dtype=np.float32); pm_dist[1:14]=1
R['S3']={}
for nm,pmk in (('seed_PAM_proximal_only',pm_seed),('distal_only',pm_dist)):
    R['S3'][nm]=agg([evaluate(np.where(~np.isin(g,fg))[0],np.where(np.isin(g,fg))[0],S,O,y,seed=k,posmask=pmk,only=('B2',)) for k,fg in enumerate(folds)])['B2']
allidx=np.arange(len(y)); ti_r=np.random.default_rng(5); R['S4']={}
for nm,(S2,O2,y2) in (('kleinstiver',(Sk,Ok,yk)),('listgarten',(Sl,Ol,yl))):
    R['S4'][nm]=evaluate(allidx,np.arange(len(y2)),np.concatenate([S,S2]),np.concatenate([O,O2]),np.concatenate([y,y2]),seed=0) if False else None
    ti=train_idx(allidx,y,10,ti_r); lr=LogisticRegression(max_iter=300).fit(mm(S[ti],O[ti]),y[ti]); f=fitM(S[ti],O[ti],y[ti],0)
    R['S4'][nm]=dict(B1=met(y2,-mm(S2,O2).sum(1)),B2=met(y2,lr.decision_function(mm(S2,O2))),M=met(y2,f(S2,O2)),n_pos=int(y2.sum())); print(nm,R['S4'][nm],flush=True)
te=np.where(np.isin(g,folds[0]))[0]; tr=np.where(~np.isin(g,folds[0]))[0]
R['S5_shuffle_control_M_fold0']=evaluate(tr,te,S,O,y,seed=0,shuffle=True,only=('M',))['M']
R['S5_ratio_sensitivity_B2_fold0']={str(rt):evaluate(tr,te,S,O,y,ratio=rt,seed=0,only=('B2',))['B2'] for rt in (1,10,50)}
json.dump(R,open('results/results.json','w'),indent=1); print(json.dumps({k:v for k,v in R.items() if k!='S1_folds'},indent=1))
