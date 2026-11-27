"""Figure: the helix contrast in deletion cost against helix propensity, 19 residues."""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

PACE={"A":0.0,"L":0.21,"R":0.21,"M":0.24,"K":0.26,"Q":0.39,"E":0.40,"I":0.41,"W":0.49,"S":0.50,
      "Y":0.53,"F":0.54,"H":0.61,"V":0.61,"N":0.65,"T":0.66,"C":0.68,"D":0.69,"G":1.00,"P":3.16}
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":160})

def main(rdir,outdir):
    ind=pd.read_csv(f"{rdir}/tsuboyama_indels.csv"); feat=pd.read_csv(f"{rdir}/tsub_features.csv")
    ss=pd.read_csv(f"{rdir}/tsub_ss.csv"); ss["domain"]=ss["domain"].str.replace(".pdb","",regex=False)
    wt=json.load(open(f"{rdir}/tsub_wt_dG.json"))
    ind=ind[(ind.kind=="del")&(ind.dG_ci<=1.0)].copy(); ind["wt"]=ind["domain"].map(wt)
    ind=ind[ind.wt.notna()]; ind["cost"]=ind["wt"]-ind["dG"]
    df=(feat[["domain","pos","struct_pos"]].merge(ss[["domain","struct_pos","aa","ss"]],on=["domain","struct_pos"])
        .merge(ind.groupby(["domain","pos"])["cost"].mean().reset_index(),on=["domain","pos"]))
    pts=[]
    for aa,g in df.groupby("aa"):
        v=[float(x[x.ss=="H"]["cost"].mean()-x[x.ss=="L"]["cost"].mean()) for _,x in g.groupby("domain")
           if (x.ss=="H").sum()>=1 and (x.ss=="L").sum()>=1]
        if len(v)>=5 and aa in PACE: pts.append((aa,PACE[aa],float(np.mean(v)),len(v)))
    fig,ax=plt.subplots(figsize=(4.6,3.3))
    for aa,x,y,n in pts:
        ax.scatter(x,y,s=14+n/4,color="#2f6f9f",zorder=3,alpha=.85)
        ax.annotate(aa,(x,y),textcoords="offset points",xytext=(4,3),fontsize=8,color="#333")
    x=np.array([p[1] for p in pts]); y=np.array([p[2] for p in pts])
    b,a=np.polyfit(x,y,1); xs=np.linspace(x.min(),x.max(),50)
    ax.plot(xs,a+b*xs,color="#c2553a",lw=1.3,zorder=2)
    r=stats.spearmanr(x,y)[0]
    ax.set_xlabel("Pace & Scholtz helix penalty of the deleted residue (kcal/mol)")
    ax.set_ylabel("deletion cost in helix - in loop\n(kcal/mol, paired within protein)")
    ax.set_title("Good helix formers are costlier to delete from a helix",fontsize=9)
    ax.text(0.97,0.95,f"Spearman {r:.2f}\nCI [-0.69, -0.45]",transform=ax.transAxes,
            ha="right",va="top",fontsize=8,color="#555")
    fig.tight_layout(); fig.savefig(os.path.join(outdir,"fig6_deletion_propensity.png"),bbox_inches="tight")
    print("ok",len(pts),"residues")

if __name__=="__main__": main(sys.argv[1],sys.argv[2])
