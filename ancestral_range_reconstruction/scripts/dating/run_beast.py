"""Print or explicitly execute a fresh 18-million-state run, not a historical resume."""
from pathlib import Path
import argparse, subprocess
p=argparse.ArgumentParser();p.add_argument('--repeat',type=int,choices=[1,2],required=True);p.add_argument('--beast',default='beast');p.add_argument('--execute',action='store_true');a=p.parse_args()
root=Path(__file__).resolve().parents[2];out=root/f'output/dating/fresh_repeat{a.repeat}'
xml=root/f'data/dating/fresh_run_xml/Walbeck_repeat{a.repeat}_18m.xml'
cmd=[a.beast,'-seed',str(20262899+a.repeat),'-threads','4',str(xml)]
print('Working directory:',out);print(subprocess.list2cmdline(cmd))
if a.execute:
    if out.exists():raise SystemExit('Run directory already exists; refusing to overwrite')
    out.mkdir(parents=True);subprocess.run(cmd,cwd=out,check=True)
