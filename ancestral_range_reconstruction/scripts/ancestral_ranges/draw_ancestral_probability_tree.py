"""Traceable range-probability tree: fixed representative tree, six-model AIC weighting."""
from pathlib import Path
import json, hashlib, argparse
import numpy as np
import pandas as pd
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Patch, Rectangle
from matplotlib.offsetbox import DrawingArea, AnnotationBbox

BASE=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser()
parser.add_argument('--with-outgroups',action='store_true')
args=parser.parse_args()
DATA=BASE/'data/ancestral_ranges/frozen_fits'
OUT=BASE/'output/ancestral_ranges/recomputed_figure';OUT.mkdir(parents=True,exist_ok=True)
import shutil
for support_name in ['node_age_intervals.csv','display_tree_with_outgroups.nwk','display_outgroup_selection.csv']:
    shutil.copy2(BASE/'output/ancestral_ranges/frozen'/support_name, OUT/support_name)
assert (DATA/'completed.txt').exists()
models=pd.read_csv(DATA/'model_summary.csv')
nodes=pd.read_csv(DATA/'node_probabilities.csv')
treepath=BASE/'data/ancestral_ranges/representative_21_tip_tree.nwk'
tree=Phylo.read(treepath,'newick');tree.ladderize()
height=max(tree.depths().values())
keys={c:'|'.join(sorted(t.name for t in c.get_terminals())) for c in tree.get_nonterminals()}
weights=models[['input','model','AIC_weight']]
merged=nodes.merge(weights,on=['input','model'],validate='many_to_one')
merged['weighted_probability']=merged.probability*merged.AIC_weight
averaged=merged.groupby(['input','clade_key','area_state'],as_index=False).weighted_probability.sum().rename(columns={'weighted_probability':'probability'})
check=averaged.groupby(['input','clade_key']).probability.sum()
assert np.allclose(check,1,atol=1e-7)
averaged.to_csv(OUT/'model_weighted_probabilities.csv',index=False)
source='Walbeck_minimum_all_dated'
primary=averaged[averaged.input==source].copy()
assert set(primary.clade_key)==set(keys.values()), 'Clade mapping must match the 21-tip tree exactly'
primary.to_csv(OUT/'primary_node_probabilities.csv',index=False)
models[models.input==source].to_csv(OUT/'primary_model_weights.csv',index=False)
stats=[]
for i,c in enumerate(tree.get_nonterminals(),1):
    p=primary[primary.clade_key==keys[c]].set_index('area_state').probability
    stats.append(dict(node_label=f'N{i:02}',clade_key=keys[c],age_Ma=height-tree.depths()[c],top_range=p.idxmax(),top_probability=p.max()))
pd.DataFrame(stats).to_csv(OUT/'node_labels_and_ages.csv',index=False)
Phylo.write(tree,OUT/'display_tree.nwk','newick',format_branch_length='%1.12f')
states=['EA','SA','NAM','EU','EA+SA','EA+NAM','EA+EU','SA+NAM','SA+EU','NAM+EU','EA+SA+NAM','EA+SA+EU','EA+NAM+EU','SA+NAM+EU','EA+SA+NAM+EU','empty']
palette=['#4477AA','#EE7733','#228833','#CC3311','#66CCEE','#009988','#EE99AA','#BBBB33','#994455','#AA4499','#88CCAA','#DDCC77','#332288','#882255','#777777','#FFFFFF']
colors=dict(zip(states,palette))
nodeids={c:f'N{i:02}' for i,c in enumerate(tree.get_nonterminals(),1)}
age_intervals=pd.read_csv(OUT/'node_age_intervals.csv').set_index('node_label')
assert set(age_intervals.index)==set(nodeids.values())
for c, label in nodeids.items():
    assert age_intervals.loc[label,'clade_key']==keys[c]
def annotated(c):
    if c.is_terminal():return c.name+f':{c.branch_length:.12f}'
    p=primary[primary.clade_key==keys[c]].set_index('area_state').probability
    quoted_states=','.join('"'+s.replace('NAM','NA')+'"' for s in states)
    values=','.join(f'{p[s]:.12g}' for s in states)
    interval=age_intervals.loc[nodeids[c]]
    age_note=(f',height_95%_HPD={{{interval.HPD95_lower_Ma:.12g},{interval.HPD95_upper_Ma:.12g}}}'
              if interval.pass_check else '')
    comment=(f'[&node_id="{nodeids[c]}",height={height-tree.depths()[c]:.12g},'
             f'range="{p.idxmax().replace("NAM","NA")}",range.prob={p.max():.12g},'
             f'range.set={{{quoted_states}}},range.set.prob={{{values}}}'+age_note+']')
    branch=':0.000000000000' if c is tree.root else f':{c.branch_length:.12f}'
    return '('+','.join(annotated(child) for child in c.clades)+')'+nodeids[c]+comment+branch
newick=annotated(tree.root)+';'
(OUT/'ancestral_probability_tree_annotated.nwk').write_text(newick+'\n')
(OUT/'ancestral_probability_tree_annotated.nex').write_text('#NEXUS\n[Range probabilities are AIC-weighted marginal values conditional on the fixed tree.]\nBegin trees;\nTree ancestral_ranges = [&R] '+newick+'\nEnd;\n')
parsed=Phylo.read(OUT/'ancestral_probability_tree_annotated.nex','nexus')
assert {t.name for t in parsed.get_terminals()}=={t.name for t in tree.get_terminals()}
assert len(parsed.get_nonterminals())==20
assert np.isclose(max(parsed.depths().values()),height,atol=1e-8)
assert all(c.comment and 'range.set.prob=' in c.comment for c in parsed.get_nonterminals())
shown=[s for s in states if primary[primary.area_state==s].probability.max()>=.015]
others=[s for s in states if s not in shown]

def draw(filename, probabilities, subtitle):
    plot_tree=Phylo.read(OUT/'display_tree_with_outgroups.nwk','newick') if args.with_outgroups else tree
    depths=plot_tree.depths(); tips=plot_tree.get_terminals()
    plot_height=max(depths[t] for t in tips)
    matched={'|'.join(sorted(t.name for t in c.get_terminals())):c for c in plot_tree.get_nonterminals()}
    outgroups=(pd.read_csv(OUT/'display_outgroup_selection.csv').set_index('tree_tip').family.to_dict()
               if args.with_outgroups else {})
    fig,ax=plt.subplots(figsize=(16,14.5) if args.with_outgroups else (12.5,11.8))
    ypos={t:len(tips)-1-i for i,t in enumerate(tips)}
    def assign(c):
        if c not in ypos:ypos[c]=(assign(c.clades[0])+assign(c.clades[-1]))/2
        return ypos[c]
    assign(plot_tree.root)
    for c in tree.get_nonterminals():
        pc=matched[keys[c]]
        interval=age_intervals.loc[nodeids[c]]
        if interval.pass_check:
            ax.plot([-interval.HPD95_upper_Ma,-interval.HPD95_lower_Ma],
                    [ypos[pc],ypos[pc]],color='#8EAFC4',lw=5,alpha=.65,
                    solid_capstyle='butt',zorder=2)
    for c in plot_tree.get_nonterminals():
        x=depths[c]-plot_height
        ax.plot([x,x],[ypos[c.clades[0]],ypos[c.clades[-1]]],color='#444444',lw=1,zorder=1)
        for child in c.clades:ax.plot([x,depths[child]-plot_height],[ypos[child],ypos[child]],color='#444444',lw=1,zorder=1)
    for t in tips:
        label=t.name.replace('_Cui_MT898426','').replace('_',' ')
        if t.name in outgroups:label+=f'  ({outgroups[t.name]})'
        ax.text(2.5 if args.with_outgroups else 1.7,ypos[t],label,va='center',fontsize=9 if args.with_outgroups else 10,fontstyle='italic',color='#666666' if t.name in outgroups else '#111111')
        area='SA' if t.name in ['Boquila_trifoliata','Lardizabala_biternata'] else 'EA'
        ax.scatter([0],[ypos[t]],s=30,c='#888888' if t.name in outgroups else colors[area],edgecolors='white',linewidths=.4,zorder=3)
    for i,c in enumerate(tree.get_nonterminals(),1):
        pc=matched[keys[c]]
        p=probabilities[probabilities.clade_key==keys[c]].set_index('area_state').probability
        assert abs(p.sum()-1)<1e-7
        diameter=19 if i<=6 else 14
        center=diameter/2;radius=center-.5
        da=DrawingArea(diameter,diameter,0,0);angle=90.
        for state in shown+['other']:
            prob=float(p.get(state,0)) if state!='other' else float(p.reindex(others,fill_value=0).sum())
            if prob<1e-10:continue
            da.add_artist(Wedge((center,center),radius,angle,angle+360*prob,facecolor=colors.get(state,'#DDDDDD'),edgecolor='none'))
            angle+=360*prob
        da.add_artist(Circle((center,center),radius,fill=False,edgecolor='white',linewidth=.5))
        xy=(depths[pc]-plot_height,ypos[pc]);ax.add_artist(AnnotationBbox(da,xy,frameon=False,zorder=4))
    # Historical observations are marks on genus branches, not inferred new species.
    obs=BASE/'data/ancestral_ranges/Walbeck_minimum_all_dated.nwk'
    ot=Phylo.read(obs,'newick');oh=max(ot.depths().values());od=ot.depths()
    for tip in ot.get_terminals():
        if not tip.name.startswith('OBS_'):continue
        genus=tip.name.split('_')[1]
        modern=next(t for t in tips if t.name.startswith(genus+'_'))
        age=oh-od[tip]
        ax.scatter([-age],[ypos[modern]],marker='|',s=95,color='#111111',linewidths=1.2,zorder=5)
    root=probabilities[probabilities.clade_key==keys[tree.root]].set_index('area_state').probability
    top=root.idxmax().replace('NAM','NA')
    ax.text(-plot_height,len(tips)+1.0,f'Lardizabalaceae crown (N01): {top} ({root.max():.3f})',fontsize=11,weight='bold')
    oldest=max(plot_height,age_intervals.loc[age_intervals.pass_check,'HPD95_upper_Ma'].max())
    ax.set_xlim(-oldest-3,85 if args.with_outgroups else 32);ax.set_ylim(-1.2,len(tips)+2.2)
    ticks=list(range(0,int(oldest)+1,20 if args.with_outgroups else 10))
    ax.set_xticks([-x for x in ticks],[str(x) for x in ticks]);ax.set_xlabel('Time before present (Ma)',fontsize=10,labelpad=4)
    # ICS 2024/12 boundaries (Ma); terminal units abbreviated to avoid overlap.
    time_units=[('Cretaceous',143.1,66.,'#B6D99B'),('Paleocene',66.,56.,'#F3B77C'),
                ('Eocene',56.,33.9,'#F7C790'),('Oligocene',33.9,23.04,'#F9D9AF'),
                ('Miocene',23.04,5.333,'#F3E995'),('Pl.',5.333,2.58,'#F8EDBD'),
                ('Q.',2.58,0.,'#EEE8CF')]
    for label,old,young,color in time_units:
        left=min(old,oldest);right=max(young,0.)
        if left<=right:continue
        ax.add_patch(Rectangle((-left,-.115),left-right,.040,
                               transform=ax.get_xaxis_transform(),facecolor=color,
                               edgecolor='#555555',linewidth=.6,clip_on=False))
        ax.text(-(left+right)/2,-.095,label,transform=ax.get_xaxis_transform(),
                ha='center',va='center',fontsize=7 if args.with_outgroups else 6.5,clip_on=False)
    ax.set_yticks([]);ax.spines[['left','top','right']].set_visible(False)
    ax.spines['bottom'].set_bounds(-oldest,0)
    ax.set_title('Lardizabalaceae ancestral range reconstruction\n'+subtitle,fontsize=13,pad=12)
    handles=[Patch(facecolor=colors[s],label=s.replace('NAM','NA')) for s in shown]
    if others:handles.append(Patch(facecolor='#DDDDDD',label='Other ranges'))
    if args.with_outgroups:handles.append(Patch(facecolor='#888888',label='Outgroup (display only)'))
    fig.legend(handles=handles,loc='lower center',ncol=5,fontsize=9,frameon=False,bbox_to_anchor=(.5,.052))
    fig.text(.07,.033,'EA = East Asia; SA = South America; NA = North America; EU = Europe.',fontsize=8.5)
    fig.text(.07,.117,'Geological time: ICS 2024/12. Pl. = Pliocene; Q. = Quaternary (Pleistocene and Holocene).',fontsize=8.0)
    fig.text(.07,.017,'Blue bars: 95% HPD intervals for node MRCA ages (diagnostically passing nodes). Pies: conditional range probabilities.',fontsize=8.0)
    fig.text(.07,.003,'Black ticks: assumed genus observations. Outgroups: dating context only; no range pies or age bars. Full-tree diagnostic limits remain.' if args.with_outgroups else 'Black ticks: dated genus observations under assumed lineage placements. Full-tree clock diagnostics retain limitations.',fontsize=8.0)
    fig.subplots_adjust(left=.07,right=.985,top=.92,bottom=.23)
    # Label the ages of this selected tree, not posterior median ages. Position
    # labels in display coordinates to keep the recent radiation readable.
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    occupied=[]
    from matplotlib.transforms import Bbox
    for c in tree.get_nonterminals():
        pc=matched[keys[c]]
        px,py=ax.transData.transform((depths[pc]-plot_height,ypos[pc]))
        radius=(10 if int(nodeids[c][1:])<=6 else 7.5)*fig.dpi/72
        occupied.append(Bbox.from_bounds(px-radius,py-radius,2*radius,2*radius))
    age_rows=[]
    label_by_key={keys[c]:nodeids[c] for c in tree.get_nonterminals()}
    for pc in plot_tree.get_nonterminals():
        key='|'.join(sorted(t.name for t in pc.get_terminals()))
        node_id=label_by_key.get(key,'')
        age=plot_height-depths[pc]
        label=(f'{node_id}: ' if node_id and int(node_id[1:])<=6 else '')+f'{age:.2f}'
        xy=(depths[pc]-plot_height,ypos[pc])
        annotation=ax.annotate(label,xy,xytext=(-12,11),textcoords='offset points',
            fontsize=7.2,color='#333333',ha='right',va='center',zorder=7,
            bbox=dict(boxstyle='round,pad=0.12',fc='white',ec='none',alpha=.9),
            arrowprops=dict(arrowstyle='-',color='#777777',lw=.45,shrinkA=2,shrinkB=7))
        # Prefer a nearby label above the node; move left for dense tip nodes.
        placed=False
        for dx in (-12,-38,-64,-90,-116):
            for dy in (11,-12,23,-24,35,-36):
                annotation.set_position((dx,dy))
                annotation.update_positions(renderer)
                annotation.update_bbox_position_size(renderer)
                box=annotation.get_bbox_patch().get_window_extent(renderer).expanded(1.08,1.18)
                if not any(box.overlaps(other) for other in occupied):
                    occupied.append(box)
                    placed=True
                    break
            if placed:break
        assert placed, f'No nonoverlapping age-label position for {key}'
        age_rows.append(dict(node_label=node_id,clade_key=key,age_Ma=age,
                             displayed_age=f'{age:.2f}',label=label))
    pd.DataFrame(age_rows).to_csv(OUT/f'{filename}_node_ages.csv',index=False)
    fig.text(.07,.132,'Node labels: ages (Ma) on the representative tree, rounded to two decimals; N01-N06 identify key nodes.',fontsize=8.0)
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{filename}.{ext}',dpi=240)
    plt.close(fig)

suffix='_with_outgroups_HPD' if args.with_outgroups else '_with_HPD'
draw('ancestral_probability_tree'+suffix,primary,'Six-model AIC-weighted probabilities; Walbeck-calibrated representative tree')
root=primary[primary.clade_key==keys[tree.root]].sort_values('probability',ascending=False)
metadata={'tree_sha256':hashlib.sha256(treepath.read_bytes()).hexdigest(),'input':source,'tip_count':len(tree.get_terminals()),'internal_nodes':len(keys),'tree_role':'sample MCC with original sampled heights, not posterior median heights','probability_role':'AIC-weighted marginal ancestral range probabilities on a fixed tree; not joint Bayesian posterior probabilities','root_probabilities':dict(zip(root.area_state,root.probability))}
metadata['age_intervals']='95% HPD of MRCA ages from both complete late Walbeck chains; show only passing node-age diagnostics; see node_age_intervals.csv'
metadata['geological_time_scale']={'source':'https://stratigraphy.org/ICSchart/ChronostratChart2024-12.pdf','boundaries_Ma':[66.,56.,33.9,23.04,5.333,2.58,0.],'abbreviations':{'Pl.':'Pliocene','Q.':'Quaternary (Pleistocene and Holocene)'},'oldest_visible_unit':'Cretaceous truncated at displayed-tree age'}
if args.with_outgroups:
    displaypath=OUT/'display_tree_with_outgroups.nwk'
    metadata['display_tree']={'path':displaypath.name,'sha256':hashlib.sha256(displaypath.read_bytes()).hexdigest(),'tip_count':23,'outgroups':['Circaeaster_agrestis','Kingdonia_uniflora'],'role':'same posterior sample, pruned display only; no outgroup range inference','root_role':'MRCA of sampled Lardizabalaceae and Circaeasteraceae; distinct from ingroup crown N01'}
metadata['primary_figure']='ancestral_probability_tree'+suffix
metadata['node_age_labels']={'unit':'Ma','decimal_places':2,'scope':'all displayed internal nodes, including outgroup nodes','source':'original node heights of the displayed representative tree; not posterior median ages','table':'ancestral_probability_tree'+suffix+'_node_ages.csv'}
(OUT/'figure_provenance.json').write_text(json.dumps(metadata,indent=2))
print(root.head(7).to_string(index=False))
