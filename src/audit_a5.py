import json
src=open('src/run.py').read().split("guides=sorted(set(g))")[0]
exec(src)
guides=sorted(set(g)); r0=np.random.default_rng(1); r0.shuffle(guides); folds=[guides[i::5] for i in range(5)]
R={}
for k in (4,0):
    fg=folds[k]; te=np.where(np.isin(g,fg))[0]; tr=np.where(~np.isin(g,fg))[0]; out=[]
    for sd in range(10,18):
        ti=train_idx(tr,y,10,np.random.default_rng(sd)); f=fitM(S[ti],O[ti],y[ti],sd,shuffle=True); out.append(met(y[te],f(S[te],O[te]))); print('A5',k,sd,out[-1],flush=True)
    a=[o['auroc'] for o in out]; R[k]=dict(runs=out,auroc_mean=float(np.mean(a)),auroc_min=float(min(a)),auroc_max=float(max(a)),auroc_sd=float(np.std(a)))
    print(k,R[k]['auroc_mean'],R[k]['auroc_min'],R[k]['auroc_max'],flush=True)
json.dump(R,open('results/audit_a5.json','w'),indent=1); print('done')
