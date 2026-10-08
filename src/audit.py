import json
src=open('src/run.py').read().split("guides=sorted(set(g))")[0]
exec(src)
import pandas as pd
d=pd.read_csv(D+'I_1_CIRCLE_seq_10gRNA_wholeDataset.csv')
guides=sorted(set(g)); r0=np.random.default_rng(1); r0.shuffle(guides); folds=[guides[i::5] for i in range(5)]
gseq={gg:d.loc[d['sgRNA_type']==gg,'sgRNA_seq'].iloc[0] for gg in guides}
print('guides per fold',folds, flush=True)
def ham(a,b): return sum(x!=y for x,y in zip(a,b))
pairkey=(d['sgRNA_seq']+'|'+d['off_seq']).values; offk=d['off_seq'].values
R={'A1':{},'A2':{},'A3':{},'A4':{}}
for k,fg in enumerate(folds):
    te=np.isin(g,fg); trg=[x for x in guides if x not in fg]
    dists=[ham(gseq[a],gseq[b]) for a in fg for b in trg]
    R['A1'][k]=dict(test_guides=fg,min_hamming=int(min(dists)),identical=int(sum(gseq[a]==gseq[b] for a in fg for b in trg)))
    tp=te&(y==1); trpairs=set(pairkey[~te]); trpos_pairs=set(pairkey[(~te)&(y==1)]); troff=set(offk[~te]); troffpos=set(offk[(~te)&(y==1)])
    R['A2'][k]=dict(n_test_pos=int(tp.sum()),pair_in_train=float(np.mean([p in trpairs for p in pairkey[tp]])),pair_in_train_pos=float(np.mean([p in trpos_pairs for p in pairkey[tp]])),
        off_in_train=float(np.mean([p in troff for p in offk[tp]])),off_in_train_pos=float(np.mean([p in troffpos for p in offk[tp]])))
    print(k,R['A1'][k],R['A2'][k],flush=True)
for k,fg in enumerate(folds):
    te=np.where(np.isin(g,fg))[0]; b1=met(y[te],-mm(S[te],O[te]).sum(1))['auroc']; un=[]
    for sd in range(5):
        torch.manual_seed(100+sd); net=Net().eval()
        with torch.no_grad(): sc=np.concatenate([net(torch.tensor(onehot(S[te][i:i+20000],O[te][i:i+20000]))).numpy() for i in range(0,len(te),20000)])
        un.append(met(y[te],sc))
    R['A3'][k]=dict(B1_auroc=b1,untrained_auroc=[u['auroc'] for u in un],untrained_auroc_mean=float(np.mean([u['auroc'] for u in un])),untrained_auprc_mean=float(np.mean([u['auprc'] for u in un])))
    print('A3',k,R['A3'][k],flush=True)
a=[R['A3'][k]['B1_auroc'] for k in range(5)]; b=[R['A3'][k]['untrained_auroc_mean'] for k in range(5)]
R['A3_corr_untrained_vs_B1']=float(np.corrcoef(a,b)[0,1]); print('corr',R['A3_corr_untrained_vs_B1'],flush=True)
def fit_within(S_,O_,y_,gg,seed):
    r=np.random.default_rng(seed); yt=y_.copy()
    for x in np.unique(gg): m=gg==x; yt[m]=r.permutation(yt[m])
    return fitM(S_,O_,yt,seed,shuffle=False)
for k in (4,0):
    fg=folds[k]; te=np.where(np.isin(g,fg))[0]; tr=np.where(~np.isin(g,fg))[0]; out=[]
    for sd in range(3):
        ti=train_idx(tr,y,10,np.random.default_rng(sd)); f=fit_within(S[ti],O[ti],y[ti],g[ti],sd); out.append(met(y[te],f(S[te],O[te])))
        print('A4',k,sd,out[-1],flush=True)
    R['A4'][k]=out
json.dump(R,open('results/audit_fold4.json','w'),indent=1,default=float); print('done')
