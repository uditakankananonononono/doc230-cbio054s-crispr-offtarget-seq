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

import sys
SEEDX=int(sys.argv[1]); res=[]
for k,fg in enumerate(folds):
    te=np.where(np.isin(g,fg))[0]; tr=np.where(~np.isin(g,fg))[0]
    r=evaluate(tr,te,S,O,y,seed=k+10*SEEDX,shuffle=True,only=('M',))['M']; res.append(r); print('shuffle fold',k,r,flush=True)
out=dict(per_fold=res,mean_auroc=round(float(np.mean([r['auroc'] for r in res])),4),mean_auprc=round(float(np.mean([r['auprc'] for r in res])),4),mean_base=round(float(np.mean([r['base'] for r in res])),5))
json.dump(out,open('results/shuffle_control_5fold_seed%d.json'%SEEDX,'w'),indent=1); print(out)
