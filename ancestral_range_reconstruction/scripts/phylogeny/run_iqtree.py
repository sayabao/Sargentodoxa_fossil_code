"""Opt-in fresh ML inference; never overwrites frozen outputs."""
from pathlib import Path
import argparse, shutil, subprocess
p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--iqtree',default='iqtree3');a=p.parse_args()
root=Path(__file__).resolve().parents[2]
out=root/'output/phylogeny/fresh_run'
cmd=[a.iqtree,'-s',str(root/'data/phylogeny/concatenated.fasta'),'-m','GTR+G4','-B','1000','--alrt','1000','--bnni','-T','2','--seed','20261005','-o','Euptelea_polyandra','-p',str(root/'data/phylogeny/partitions.nex'),'--prefix',str(out/'primary50')]
print(subprocess.list2cmdline(cmd))
if a.execute:
    if out.exists():raise SystemExit('fresh_run already exists; choose a separate directory before rerunning')
    out.mkdir(parents=True)
    subprocess.run(cmd,check=True)
