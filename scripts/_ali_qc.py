import scanpy as sc, pandas as pd, numpy as np, glob, os, re
sc.settings.verbosity=0
rows=[]
for f in sorted(glob.glob("data/raw/ali_gse121600/*.h5")):
    s=os.path.basename(f).replace(".h5","")
    a=sc.read_10x_h5(f); a.var_names_make_unique()
    n_raw=a.n_obs
    sc.pp.filter_cells(a,min_counts=1000); sc.pp.filter_cells(a,min_genes=500)
    mt=[g for g in a.var_names if g.startswith("MT-")]
    if mt: a.obs["pmt"]=np.asarray(a[:,mt].X.sum(1)).ravel()/np.asarray(a.X.sum(1)).ravel()*100
    else: a.obs["pmt"]=0.0
    a=a[a.obs.pmt<25].copy()
    g=lambda x: float(np.asarray(a[:,x].X.sum(0)).ravel()[0]) if x in a.var_names else np.nan
    frac=lambda x: float((np.asarray(a[:,x].X.todense()).ravel()>0).mean()) if x in a.var_names else np.nan
    rows.append(dict(sample=s, barcodes_raw=n_raw, cells=a.n_obs,
        median_umi=float(np.median(a.obs.n_counts)), median_genes=float(np.median(a.obs.n_genes)),
        KRT5=frac("KRT5"), SCGB1A1=frac("SCGB1A1"), SCGB3A1=frac("SCGB3A1"),
        MUC5AC=frac("MUC5AC"), MUC5B=frac("MUC5B"), FOXJ1=frac("FOXJ1"), TP63=frac("TP63"),
        n_genes_detected=int((np.asarray(a.X.sum(0)).ravel()>0).sum())))
    print("done",s,a.n_obs,flush=True)
df=pd.DataFrame(rows); df.to_csv("results/ali_qc.csv",index=False)
print(df.round(3).to_string(index=False))
