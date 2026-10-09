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

import time
t0=time.time()
def ece(p,yv,nb=10):
    o=np.argsort(p); bins=np.array_split(o,nb); return float(sum(len(b)/len(p)*abs(p[b].mean()-yv[b].mean()) for b in bins))
sig=lambda z:1/(1+np.exp(-z))
sc={k:[] for k in ('M','B2','M25','M50','Msh','y','mm','g','fold')}
for k,fg in enumerate(folds):
    te=np.where(np.isin(g,fg))[0]; tr=np.where(~np.isin(g,fg))[0]; r=np.random.default_rng(k)
    ti=train_idx(tr,y,10,r)
    f=fitM(S[ti],O[ti],y[ti],0); sc['M'].append(f(S[te],O[te]))
    lr=LogisticRegression(C=1,max_iter=300).fit(mm(S[ti],O[ti]),y[ti]); sc['B2'].append(lr.decision_function(mm(S[te],O[te])))
    for frac,key in ((0.25,'M25'),(0.5,'M50')):
        gs_=sorted(set(g[tr])); rr=np.random.default_rng(0); keep=rr.choice(gs_,max(1,int(round(frac*len(gs_)))),replace=False)
        sub=tr[np.isin(g[tr],keep)]; ti2=train_idx(sub,y,10,np.random.default_rng(k)); f2=fitM(S[ti2],O[ti2],y[ti2],0); sc[key].append(f2(S[te],O[te]))
    Sp=S[ti].copy(); Sp=Sp[np.random.default_rng(0).permutation(len(Sp))]; f3=fitM(Sp,O[ti],y[ti],0); sc['Msh'].append(f3(S[te],O[te]))
    sc['y'].append(y[te]); sc['mm'].append(mm(S[te],O[te]).sum(1)); sc['g'].append(g[te]); sc['fold'].append(np.full(len(te),k))
    print('fold',k,round(time.time()-t0),flush=True)
sc={k:np.concatenate(v) for k,v in sc.items()}; yy=sc['y']; R={}
ap=lambda yv,s:float(average_precision_score(yv,s))
# E6
R['E6']={}
for nm,mask in (('<=3mm',sc['mm']<=3),('>3mm',sc['mm']>3)):
    R['E6'][nm]=dict(n=int(mask.sum()),pos=int(yy[mask].sum()),M=ap(yy[mask],sc['M'][mask]) if yy[mask].sum()>0 else None,B2=ap(yy[mask],sc['B2'][mask]) if yy[mask].sum()>0 else None)
a=R['E6']['<=3mm']; R['E6']['gate_met']=bool(a['M'] is not None and a['M']>=a['B2']+0.02); print('E6',R['E6'],flush=True)
# E7
gl=sorted(set(sc['g'])); wins=0; per={}
for gg in gl:
    m=sc['g']==gg
    if yy[m].sum()==0: continue
    per[str(gg)]=dict(M=ap(yy[m],sc['M'][m]),B2=ap(yy[m],sc['B2'][m]),pos=int(yy[m].sum())); wins+=per[str(gg)]['M']>per[str(gg)]['B2']
R['E7']=dict(per_guide=per,M_wins=int(wins),n_guides=len(per),gate_met=bool(wins>=7)); print('E7',wins,len(per),flush=True)
# E8
R['E8']=dict(auprc_25=ap(yy,sc['M25']),auprc_50=ap(yy,sc['M50']),auprc_100=ap(yy,sc['M'])); R['E8']['gate_met']=bool(R['E8']['auprc_100']>=R['E8']['auprc_25']+0.02); print('E8',R['E8'],flush=True)
# E9
p=sig(sc['M']); R['E9']=dict(ece=ece(p,yy),mean_pred=float(p.mean()),true_rate=float(yy.mean())); print('E9',R['E9'],flush=True)
# E10
R['E10']=dict(auprc_original=ap(yy,sc['M']),auprc_shuffled_guide=ap(yy,sc['Msh'])); R['E10']['drop']=R['E10']['auprc_original']-R['E10']['auprc_shuffled_guide']; R['E10']['gate_met']=bool(R['E10']['drop']>=0.05); print('E10',R['E10'],flush=True)
R['pooled_note']='AUPRC pooled over all cold-guide OOF rows (original S1 reported mean of per-fold AUPRC); pooled numbers are not directly comparable to S1'
json.dump(R,open('results/ext_results.json','w'),indent=1,default=float); print('done',round(time.time()-t0))
