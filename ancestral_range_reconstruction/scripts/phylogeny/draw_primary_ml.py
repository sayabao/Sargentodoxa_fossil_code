"""Draw the existing primary ML topology; do not infer a new tree."""
from pathlib import Path
import csv, hashlib, json, copy
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2]
SOURCE = BASE / "data/phylogeny/primary50.treefile"
OUT = BASE / "output/phylogeny/recomputed_figure"
OUT.mkdir(parents=True, exist_ok=True)
tree = Phylo.read(SOURCE, "newick")
taxa = frozenset(t.name for t in tree.get_terminals())
assert len(taxa) == 42

def key(tips):
    a = frozenset(tips); b = taxa - a
    return tuple(sorted(a if (len(a), tuple(sorted(a))) <= (len(b), tuple(sorted(b))) else b))

def splits(t):
    return {key(x.name for x in c.get_terminals()): c.name for c in t.get_nonterminals()
            if c is not t.root and 1 < len(c.get_terminals()) < len(taxa)-1}

original_splits = splits(tree)
support = dict(original_splits)
tree.root_with_outgroup("Euptelea_polyandra")
tree.ladderize(reverse=True)
assert set(splits(tree)) == set(original_splits)
for c in tree.get_nonterminals():
    c.name = support.get(key(t.name for t in c.get_terminals())) if c is not tree.root else None
    c.confidence = None
Phylo.write(tree, OUT / "primary_ml_rooted.newick", "newick", format_branch_length="%1.10f")
with (BASE / "data/phylogeny/wang2020_taxon_inventory.csv").open(encoding="utf-8-sig", newline="") as f:
    families = {r["taxon"].rstrip("."): r["family"] for r in csv.DictReader(f)}

def label(t):
    return "Sargentodoxa cuneata" if t.name.startswith("Sargentodoxa_cuneata") else t.name.replace("_", " ") + ("." if t.name.endswith("_sp") else "")

leaves=tree.get_terminals()
y={t:i for i,t in enumerate(leaves)}
x={}
for c in tree.find_clades(order="postorder"):
    if c.is_terminal(): x[c]=0
    else:
        x[c]=min(x[k] for k in c.clades)-1
        y[c]=(min(y[k] for k in c.clades)+max(y[k] for k in c.clades))/2
shift=-x[tree.root]
x={c:v+shift for c,v in x.items()}
maxx=max(x.values()); textx=maxx+.35
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9, "pdf.fonttype":42, "svg.fonttype":"none"})
fig,ax=plt.subplots(figsize=(11.7,15))
for c in tree.get_nonterminals():
    ys=[y[k] for k in c.clades]
    ax.plot([x[c],x[c]],[min(ys),max(ys)],color="black",lw=.8)
    for k in c.clades:
        ax.plot([x[c],x[k]],[y[k],y[k]],color="black",lw=.8)
        if not k.is_terminal() and k.name:
            ax.text((x[c]+x[k])/2,y[k]-.18,k.name,ha="center",va="bottom",fontsize=7,
                    bbox=dict(facecolor="white",edgecolor="none",pad=.25))
for t in leaves:
    s=label(t); highlight=s.startswith("Sargentodoxa")
    ax.text(textx,y[t],s,va="center",fontstyle="italic",fontsize=9,
            color="#9b2727" if highlight else "black",fontweight="bold" if highlight else "normal")
    if x[t]<maxx: ax.plot([x[t],maxx],[y[t],y[t]],color=".8",lw=.4,ls=":")
family_groups={}
for t in leaves: family_groups.setdefault(families[label(t).rstrip(".")],[]).append(y[t])
brace=maxx+8.3
for fam,ys in family_groups.items():
    lo,hi=min(ys)-.27,max(ys)+.27
    ax.plot([brace-.15,brace,brace,brace-.15],[lo,lo,hi,hi],color=".25",lw=.7)
    ax.text(brace+.25,(lo+hi)/2,fam,va="center",fontsize=9)
ax.text(0,-1.9,"Maximum-likelihood phylogeny of the five plastid regions",fontsize=12,weight="bold")
ax.text(0,-1.05,"Node labels: SH-aLRT / ultrafast bootstrap (%)",fontsize=9)
ax.text(0,len(leaves)+.6,"Branch lengths are not drawn to scale.",fontsize=9)
ax.set_xlim(-.3,brace+5.9);ax.set_ylim(len(leaves)+1.1,-2.4);ax.axis("off")
fig.subplots_adjust(left=.025,right=.995,top=.99,bottom=.025)
for suffix in ["pdf","svg","png"]:
    fig.savefig(OUT / ("primary_ml_tree_support."+suffix),dpi=400,facecolor="white")
plt.close(fig)
rows=[]
for c in tree.get_nonterminals():
    if c.name:
        a,b=c.name.split("/")
        rows.append(dict(descendants="|".join(sorted(t.name for t in c.get_terminals())),SH_aLRT=float(a),UFBoot=float(b)))
with (OUT / "node_support.csv").open("w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["descendants","SH_aLRT","UFBoot"]);w.writeheader();w.writerows(rows)
assert len(rows)==len(original_splits)==39
manifest={"source":str(SOURCE.relative_to(BASE)),"source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          "tips":42,"ingroup":23,"outgroups":19,"supported_bipartitions":39,
          "root":"Euptelea polyandra","display":"cladogram; lengths not proportional; no branches collapsed",
          "support_order":"SH-aLRT / UFBoot (%)","topology_verified_against_source":True,
          "new_inference_performed":False,"alignment_sites":6399,"retained_fragments":199}
(OUT / "provenance.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
print(json.dumps(manifest,indent=2))
