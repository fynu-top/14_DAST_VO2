
#!/usr/bin/env python3
"""Two-stage NO-Au Janus RCWA screen; pure Python, original material folders, no exports.
Air/PMMA cap/VO2 squares/SiO2 (aligned squares OR continuous)/ITO/DAST/ITO/SiO2 half-space.
Same-height spacer comparison, NOT equal-volume comparison. Trial effective scalar materials.
Python 3.9+; numpy scipy pandas openpyxl matplotlib torch torcwa (API 0.1.4.2).
Run beside optical_constants, Atmospheric data, VO2_temperature_data. No other .py required.
"""
from pathlib import Path
import argparse,time
import numpy as np,pandas as pd,torch,torcwa
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
# ===== All lengths nm; only change these for a different exploratory run =====
flag=0; STEP_NM=40.; ORDER=2; NX=NY=32; SHOW_FIGURE=True
P=400.; WIDTH=180.; H_VO2=12.; H_PMMA_TOP=15.
H_ITO1=5.; H_DAST=95.; H_ITO2=5.; VEND=8.5; E_LIMIT=100.
SPACERS=[20.,40.,60.,80.]; WIDTHS=[140.,180.,220.]; HEIGHTS=[8.,12.,20.]
RUN_STAGE2=True; DETAIL_POINTS=7; TLUM_FLOOR=.50; MARGIN=.005; FONT_SIZE=11
EO_SIGN=1.; DAST_SOURCE='hybrid'; ITO_FILE='4TIO.xlsx'; LOSSLESS_PMMA=True
PROJECT=Path(__file__).resolve().parent; DTYPE=torch.complex128
TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
def datafile(path,w,cols=(0,1,2),unit='nm',zero_outside=False):
    if not path.is_file():raise FileNotFoundError(f'Missing data: {path}')
    f=pd.read_csv(path,header=None,comment='#',on_bad_lines='skip') if path.suffix.lower()=='.csv' else pd.read_excel(path,header=None,engine='openpyxl')
    if f.shape[1]<=max(cols):raise ValueError(f'{path.name}: missing expected columns {cols}')
    a=f.iloc[:,list(cols)].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    a=a[np.isfinite(a).all(axis=1)];a=a[np.argsort(a[:,0])];a=a[np.unique(a[:,0],return_index=True)[1]]
    if len(a)<3:raise ValueError(f'{path.name}: no usable wavelength/n/k numeric rows')
    x=a[:,0]*(1000 if unit=='um' else 1)
    if not zero_outside and (x[0]>w.min()+1e-5 or x[-1]<w.max()-1e-5):
        raise ValueError(f'{path.name}: data range {x[0]:g}-{x[-1]:g} nm does not cover {w.min():g}-{w.max():g} nm')
    return np.column_stack([np.interp(w,x,a[:,i],left=0 if zero_outside else None,right=0 if zero_outside else None) for i in range(1,a.shape[1])])

def load_data(step):
    mat=PROJECT/'optical_constants';atm=PROJECT/'Atmospheric data';vo=PROJECT/'VO2_temperature_data'
    w=np.arange(400.,2500.00001,step);w=np.unique(np.r_[w,400.,780.,2500.,430.,457.9,488.,514.5,650.,700.,720.,750.,900.,1550.,2100.])
    d={'w':w};names={'SiO2':'2SiO2.xlsx','ITO':ITO_FILE,'VI':'6VO2_I.xlsx','VM':'7VO2_M.xlsx',
                     'PMMA':'9PMMA.xlsx'}
    for key,file in names.items():
        a=datafile(mat/file,w,unit='um');n,k=a[:,0],a[:,1]
        if np.min(n)<=0 or np.min(k)<-1e-8:raise ValueError(f'{file}: invalid n/k (passive material required)')
        if (key=='SiO2' or (key=='PMMA' and LOSSLESS_PMMA)) and np.max(np.abs(k))>1e-3:print(f'WARNING: {key} imaginary index ignored for lossless external/embedding model')
        d[key]=n+1j*(np.maximum(k,0) if not (key=='SiO2' or (key=='PMMA' and LOSSLESS_PMMA)) else 0)
    d['solar']=datafile(vo/'ASTMG173.xlsx',w,(0,2),zero_outside=True)[:,0]
    vfile=atm/'CIE_Vlambda.xlsx'
    if not vfile.is_file():vfile=vo/'CIE_Vlambda.csv'
    if not vfile.is_file():raise FileNotFoundError('Need CIE_Vlambda.xlsx or CIE_Vlambda.csv')
    cols=(pd.read_excel(vfile,header=None,engine='openpyxl') if vfile.suffix.lower()=='.xlsx' else pd.read_csv(vfile,header=None)).shape[1]
    vcol=2 if cols>=4 else 1
    v=datafile(vfile,w,(0,vcol),zero_outside=True)[:,0]
    d65=atm/'CIE_std_illum_D65.csv'
    if not d65.is_file():d65=atm/'CIE_std_illum_D65.xlsx'
    d['lum']=v*datafile(d65,w,(0,1),zero_outside=True)[:,0]
    if DAST_SOURCE=='xlsx':
        a=datafile(mat/'10DAST.xlsx',w,unit='um');n=a[:,0];k=np.maximum(a[:,1],0)
    elif DAST_SOURCE=='hybrid':
        n=np.full(len(w),2.2);ix=w>650.;x=w[ix]
        ns=np.sqrt(2.078**2+1.645*533**2/(x*x-533**2));blend=.5-.5*np.cos(np.pi*np.clip((x-650)/100,0,1));n[ix]=(1-blend)*2.2+blend*ns
        alpha=np.full(len(w),1.9e5);ix=(w>430)&(w<=514.5)
        alpha[ix]=10**PchipInterpolator([430,457.9,488,514.5],np.log10([1.9e5,4.1e5,4.9e5,4.79e5]))(w[ix])
        ix=(w>514.5)&(w<=900);alpha[ix]=10**np.interp(w[ix],[514.5,700,750,900],np.log10([4.79e5,5.,1.8,.8]))
        ix=(w>900)&(w<=2100);x=w[ix]
        alpha[ix]=.7+np.exp(-.5*((x-1100)/55)**2)+1.5*np.exp(-.5*((x-1400)/70)**2)+3*np.exp(-.5*((x-1700)/95)**2)
        alpha[w>2100]=.7;k=alpha*w*1e-7/(4*np.pi)
    else:raise ValueError('DAST_SOURCE must be hybrid or xlsx')
    if np.min(n)<=0 or np.min(k)<0:raise ValueError('DAST n/k invalid')
    # r_eff wavelength interpolation: empirical scalar L1/E1 benchmark, NOT tensor crystal model
    r=PchipInterpolator([720.,750.,1550.],[530.,445.,200.])(np.clip(w,720.,1550.)) # pm/V
    d['DASTn']=n;d['DASTk']=k;d['beta']=EO_SIGN*.5*n**3*r*1e-3 # delta n = beta * V / H_DAST(nm)
    if TRAPZ(d['solar'],w)<=0 or TRAPZ(d['lum'],w)<=0:raise ValueError('Invalid solar or D65*V luminous tables')
    return d

def mask_square(width,nx,ny,device):
    def cover(n):
        x=np.arange(n)/n;a=(1-width/P)/2;b=(1+width/P)/2
        return np.maximum(0,np.minimum(x+1/n,b)-np.maximum(x,a))*n
    return torch.as_tensor(cover(nx)[:,None]*cover(ny)[None,:],device=device,dtype=DTYPE)


def channel_power(s,lam,nsub,order,direction,port,device):
    nout=(1. if direction=='forward' else nsub) if port=='reflection' else (nsub if direction=='forward' else 1.)
    if (order+1)*lam/P < nout-1e-9:raise ValueError('ORDER excludes a propagating diffraction order; increase ORDER.')
    orders=np.array([(i,j) for i in range(-order,order+1) for j in range(-order,order+1) if (i*i+j*j)*(lam/P)**2<nout*nout-1e-10],dtype=np.int64).reshape(-1,2)
    zero=np.flatnonzero((orders==0).all(axis=1))[0];powers=[];z=[]
    for incoming in ('p','s'):
        energy=torch.zeros(len(orders),device=device,dtype=torch.float64)
        for outgoing in ('p','s'):
            a=s.S_parameters(orders=orders,direction=direction,port=port,polarization=outgoing+incoming,ref_order=[0,0],power_norm=True)
            energy+=a.abs().square().reshape(-1).real
        e=energy.detach().cpu().numpy()
        if not np.isfinite(e).all():raise FloatingPointError('Nonfinite power-normalized p/s channel.')
        powers.append(e.sum());z.append(e[zero])
    return np.asarray(powers),np.asarray(z)  # TM/p and TE/s; never sum nonorthogonal x/y channels


def powers(s,lam,nsub,order,device):
    q=np.zeros((2,2,4)) # directions F/B, incident polarizations p/s, [R,T,A,T0]
    for di,direction in enumerate(('forward','backward')):
        q[di,:,0],_=channel_power(s,lam,nsub,order,direction,'reflection',device)
        q[di,:,1],q[di,:,3]=channel_power(s,lam,nsub,order,direction,'transmission',device)
        q[di,:,2]=1-q[di,:,0]-q[di,:,1]
    if not np.isfinite(q).all() or q.min()<-3e-5 or q.max()>1+3e-5:raise FloatingPointError(f'Passivity failure at {lam:g} nm: {q}; no clipping applied.')
    rec=np.max(abs(q[0,:,3]-q[1,:,3]))
    if rec>2e-6:raise FloatingPointError(f'Matched zero-order reciprocity failed: {rec:.3g}')
    return q


def safe_lambda(w,nsub,order):
    radii=np.array(sorted({i*i+j*j for i in range(-order,order+1) for j in range(-order,order+1) if i or j}))
    # Exact Rayleigh cutoff makes some library matrices singular; a tiny one-sided offset.
    if min(np.min(abs(radii*(w/P)**2-1)),np.min(abs(radii*(w/P)**2-nsub*nsub)))<1e-10:return w*(1+1e-7)
    return w


def solve_one(d,j,phase,bias,device,order,mask,geometry):
    mode,sp,width,height=geometry;w=float(d['w'][j]);ns=float(d['SiO2'][j].real);lam=safe_lambda(w,ns,order)
    pm=complex(d['PMMA'][j]**2);nv=d['VI'][j] if phase=='I' else d['VM'][j]
    nd=d['DASTn'][j]+d['beta'][j]*bias/H_DAST+1j*d['DASTk'][j]
    if nd.real<=0:raise ValueError('Nonpositive DAST real index.')
    sim=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    sim.add_output_layer(eps=ns*ns);sim.set_incident_angle(inc_ang=0.,azi_ang=0.)
    def layer(h,eps):
        if h>0:sim.add_layer(thickness=h,eps=eps)
    layer(H_PMMA_TOP,pm);layer(height,pm+mask*(complex(nv*nv)-pm))
    eps=complex(d['SiO2'][j]**2);layer(sp,pm+mask*(eps-pm) if mode=='patch' else eps)
    for h,n in [(H_ITO1,d['ITO'][j]),(H_DAST,nd),(H_ITO2,d['ITO'][j])]:layer(h,complex(n*n))
    sim.solve_global_smatrix();return powers(sim,lam,ns,order,device)

def metrics(d,q):
    u=q.mean(axis=2);w=d['w'];vis=(w>=400)&(w<=780)
    def avg(x,lum=False):
        mask=vis if lum else np.ones(len(w),dtype=bool);weights=d['lum'] if lum else d['solar']
        den=TRAPZ(weights[mask],w[mask])
        if den<=0:raise ValueError('Zero spectral normalization.')
        return float(TRAPZ(x[mask]*weights[mask],w[mask])/den)
    af,ab=avg(u[:,0,2]),avg(u[:,1,2]);tf,tb=avg(u[:,0,1],True),avg(u[:,1,1],True)
    t0f,t0b=avg(u[:,0,3],True),avg(u[:,1,3],True)
    return dict(DA=af-ab,AsF=af,AsB=ab,TlumF=tf,TlumB=tb,Tlum=min(tf,tb),Tlum0F=t0f,Tlum0B=t0b,Tlum0=min(t0f,t0b),TsolF=avg(u[:,0,1]),TsolB=avg(u[:,1,1]),DiffractedLumF=tf-t0f,DiffractedLumB=tb-t0b)


def checks(d,device,order):
    worst=0.;energy=0.;cross=0.;mask=mask_square(WIDTH,48,48,device)
    with torch.no_grad():
        for wl in (405.,415.,550.,900.,1550.):
            n=float(np.interp(wl,d['w'],d['SiO2'].real));o=min(order,2)
            for patterned in (False,True):
                s=torcwa.rcwa(freq=1/wl,order=[o,o],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
                s.add_output_layer(eps=n*n);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
                s.add_layer(thickness=20.,eps=1.48**2+mask*(1.72**2-1.48**2) if patterned else 1.)
                s.solve_global_smatrix();q=powers(s,wl,n,o,device);cross=max(cross,np.max(abs(q[0,:,3]-q[1,:,3])))
                if patterned:energy=max(energy,np.max(abs(q[:,:,2])))
                else:
                    re=((1-n)/(1+n))**2;te=4*n/(1+n)**2;worst=max(worst,np.max(abs(q[:,:,0]-re)),np.max(abs(q[:,:,1]-te)))
    print(f'CHECK Fresnel={worst:.3e}; lossless 2D energy={energy:.3e}; matched-zero-order reciprocity={cross:.3e}')
    if worst>1e-6 or energy>1e-6:raise RuntimeError('Solver checks failed.')
    return dict(Fresnel=worst,LosslessEnergy=energy,ZeroOrderReciprocity=cross)




def objectives(endpoints):
    out={};ts=[]
    for phase in ('I','M'):
        a,b=endpoints[(phase,-VEND)],endpoints[(phase,VEND)];pm=min(a['DA'],-b['DA']);mp=min(-a['DA'],b['DA'])
        out.update({f'DAminus_{phase}':a['DA'],f'DAplus_{phase}':b['DA'],f'Margin_{phase}':max(pm,mp),f'Orientation_{phase}':'+ to -' if pm>=mp else '- to +',f'VoltageSpanDA_{phase}':abs(b['DA']-a['DA'])})
        ts.extend([a['Tlum'],b['Tlum']])
    contrast=[endpoints[('I',v)]['TsolF']-endpoints[('M',v)]['TsolF'] for v in (-VEND,VEND)];ix=int(np.argmax(np.abs(contrast)))
    out.update(Tmin=min(ts),WorstMargin=min(out['Margin_I'],out['Margin_M']),PhaseTsolAbs=abs(contrast[ix]),PhaseTsolSigned=contrast[ix],PhaseTsolBias=(-VEND,VEND)[ix])
    out['PhaseAbsMargin']=max(max(min(endpoints[('I',v)]['DA'],-endpoints[('M',v)]['DA']),min(-endpoints[('I',v)]['DA'],endpoints[('M',v)]['DA'])) for v in (-VEND,VEND))
    out['Transparent']=out['Tmin']>=TLUM_FLOOR;out['BothReverse']=out['WorstMargin']>0
    out['Qualified']=out['Transparent'] and out['WorstMargin']>=MARGIN
    return out

def pareto_mask(table):
    x=table[['WorstMargin','PhaseTsolAbs']].to_numpy();ok=table.Transparent.to_numpy(bool);answer=np.zeros(len(x),bool)
    for i in np.flatnonzero(ok):answer[i]=not np.any(ok & np.all(x>=x[i],axis=1) & np.any(x>x[i],axis=1))
    return answer

def main():
    started=time.perf_counter();torch.set_num_threads(1)
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise ValueError('Use flag=0 CPU or flag=1 with available CUDA.')
    if H_DAST<=0 or VEND<=0 or VEND>10 or VEND*1000/H_DAST>E_LIMIT:raise ValueError('Voltage/field constraint violated.')
    if not 0<WIDTH<P or any(not 0<x<P for x in WIDTHS) or min(HEIGHTS+[H_VO2])<=0 or min(SPACERS)<0:raise ValueError('Invalid search dimensions.')
    if not SPACERS or STEP_NM<=0 or ORDER<1 or min(NX,NY)<=4*ORDER or DETAIL_POINTS<3:raise ValueError('Invalid numerical settings.')
    device=torch.device('cuda' if flag else 'cpu');torch.backends.cuda.matmul.allow_tf32=False;d=load_data(STEP_NM);checks(d,device,ORDER)
    print(f'NO Au; phase-dependent optical data, DAST L1/E1; P={P:g} nm; N={ORDER}, pixels={NX}x{NY}; step={STEP_NM:g} nm; V=+/-{VEND:g} V')
    print('Stage 1: two SiO2 shapes at equal height. Stage 2: VO2 width/height around two representatives. Fixed 50% floor is editable.')
    cache={};records={};masks={}
    def evaluate_state(g,phase,v):
        key=(g,phase,float(v))
        if key not in cache:
            width=g[2]
            if width not in masks:masks[width]=mask_square(width,NX,NY,device)
            with torch.no_grad():q=np.asarray([solve_one(d,j,phase,float(v),device,ORDER,masks[width],g) for j in range(len(d['w']))])
            cache[key]=(metrics(d,q),q)
        return cache[key]
    def evaluate(g,stage):
        if g in records:return records[g]
        endpoints={(ph,v):evaluate_state(g,ph,v)[0] for ph in ('I','M') for v in (-VEND,VEND)}
        row=dict(ID=len(records)+1,Stage=stage,Mode=g[0],Spacer_nm=g[1],Width_nm=g[2],VO2_nm=g[3],**objectives(endpoints));records[g]=row
        print(f"ID {row['ID']:02d} {stage}: {g}; minT={row['Tmin']:.3f}, margins I/M={row['Margin_I']:+.5f}/{row['Margin_M']:+.5f}, phase T contrast={row['PhaseTsolAbs']:.4f}; elapsed {time.perf_counter()-started:.1f}s",flush=True)
        return row
    stage1geom=[(mode,float(sp),WIDTH,H_VO2) for mode in ('patch','film') for sp in SPACERS]
    for g in stage1geom:evaluate(g,'shape')
    stage1=pd.DataFrame([records[g] for g in stage1geom]);pool=stage1[stage1.Transparent]
    if pool.empty:pool=stage1;print('No transparent stage-1 candidates: seeds are diagnostic, not qualified.')
    first=pool.sort_values(['WorstMargin','PhaseTsolAbs'],ascending=False).iloc[0];rest=pool[pool.ID!=first.ID]
    seeds=[first]
    if len(rest):seeds.append(rest.sort_values(['PhaseTsolAbs','WorstMargin'],ascending=False).iloc[0])
    seeds=pd.DataFrame(seeds);print('\nSTAGE-2 SEEDS\n'+seeds.to_string(index=False))
    if RUN_STAGE2:
        targets=list(dict.fromkeys((r.Mode,float(r.Spacer_nm),float(w),float(h)) for r in seeds.itertuples() for w in WIDTHS for h in HEIGHTS))
        print(f'Stage 2: {sum(g not in records for g in targets)} new geometries; repeated states reused.',flush=True)
        for g in targets:evaluate(g,'size')
    baseline=('patch',60.,180.,12.);evaluate(baseline,'baseline')
    table=pd.DataFrame(records.values());table['Pareto']=pareto_mask(table);bright=table[table.Transparent];eligible=bright[bright.Qualified]
    if len(eligible):chosen=eligible.sort_values(['PhaseTsolAbs','WorstMargin','Tmin'],ascending=False).iloc[0];why='meets both-phase margin and endpoint transparency; largest phase transmission contrast in this sampled pool'
    elif len(bright):chosen=bright.sort_values(['WorstMargin','PhaseTsolAbs'],ascending=False).iloc[0];why='no jointly qualified design; best worst-phase margin among transparent candidates (diagnostic)'
    else:chosen=table.sort_values(['Tmin','WorstMargin'],ascending=False).iloc[0];why='no transparent candidate; highest endpoint transparency (diagnostic)'
    g=(chosen.Mode,float(chosen.Spacer_nm),float(chosen.Width_nm),float(chosen.VO2_nm));detail=[]
    for ph in ('I','M'):
        for v in np.unique(np.r_[np.linspace(-VEND,VEND,DETAIL_POINTS),0.]):
            met,_=evaluate_state(g,ph,float(v));detail.append(dict(Phase=ph,Voltage=v,**met))
    detail=pd.DataFrame(detail);detailmin=float(detail.Tlum.min());print('\nALL CANDIDATES\n'+table.to_string(index=False));print('\nSELECTION: '+why+'\n'+chosen.to_string())
    print('\nSELECTED VOLTAGE SAMPLES\n'+detail.to_string(index=False));print(f'Min selected T over sampled voltages and both phases={detailmin:.6f}')
    if detailmin<TLUM_FLOOR:print('Selected endpoint transparency does NOT hold across all sampled interior voltages. Do not claim full-range transparency.')
    comparison=table[table.ID.isin([records[baseline]['ID'],chosen.ID])];print('\nBASELINE / SELECTED\n'+comparison.to_string(index=False))
    plt.rcParams.update({'font.size':FONT_SIZE,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in'})
    fig,axs=plt.subplots(2,4,figsize=(19,9),constrained_layout=True);ax=axs.ravel();cols={'patch':'tab:blue','film':'tab:orange'}
    for mode in ('patch','film'):
        t=stage1[stage1.Mode==mode].sort_values('Spacer_nm')
        for j,y in enumerate(['WorstMargin','Tmin','PhaseTsolAbs']):ax[j].plot(t.Spacer_nm,t[y],'-o',color=cols[mode],label=mode)
        t=table[table.Mode==mode];ax[3].scatter(t.WorstMargin,t.PhaseTsolAbs,color=cols[mode],marker='o',facecolors='none',label=mode)
        t=t[t.Pareto];ax[3].scatter(t.WorstMargin,t.PhaseTsolAbs,color=cols[mode],marker='s')
    ax[0].axhline(0,color='.6');ax[0].axhline(MARGIN,color='.6',ls=':');ax[1].axhline(TLUM_FLOOR,color='.6',ls=':')
    for j,(title,y) in enumerate([('Stage 1: weaker-phase margin','Worst-phase margin'),('Stage 1: endpoint transparency','Min luminous T'),('Stage 1: phase transmission contrast','Max endpoint |Tsol,I - Tsol,M|')]):ax[j].set(xlabel='SiO2 spacer thickness (nm)',ylabel=y);ax[j].set_title(f'({chr(97+j)}) '+title,loc='left')
    ax[3].plot(chosen.WorstMargin,chosen.PhaseTsolAbs,'r*',ms=14);ax[3].set(xlabel='Worst-phase margin',ylabel='Phase Tsol contrast');ax[3].set_title('(d) All trials; squares = feasible Pareto',loc='left')
    for ph,c in [('I','tab:blue'),('M','tab:red')]:
        t=detail[detail.Phase==ph];ax[4].plot(t.Voltage,t.DA,'-o',color=c,label=ph);ax[5].plot(t.Voltage,t.Tlum,'-o',color=c,label=ph);ax[6].plot(t.Voltage,t.TsolF,'-o',color=c,label=ph)
        for v,ls in [(-VEND,'--'),(VEND,'-')]:
            q=evaluate_state(g,ph,v)[1].mean(axis=2);ax[7].plot(d['w'],q[:,0,2]-q[:,1,2],ls,color=c,label=f'{ph} {v:+g}V')
    ax[4].axhline(0,color='.6');ax[4].axhline(MARGIN,color='.6',ls=':');ax[4].axhline(-MARGIN,color='.6',ls=':');ax[5].axhline(TLUM_FLOOR,color='.6',ls=':')
    for j,title,y in [(4,'Selected: directional response','Solar-weighted A_F - A_B'),(5,'Selected: transparency','Min luminous T'),(6,'Selected: solar transmission','Forward solar T')]:ax[j].set(xlabel='Voltage (V)',ylabel=y);ax[j].set_title(f'({chr(97+j)}) '+title,loc='left')
    ax[7].set(xlabel='Wavelength (nm)',ylabel='A_F - A_B');ax[7].set_title('(h) Selected: endpoint spectra',loc='left')
    for a in ax:a.legend(fontsize=8);a.grid(alpha=.15)
    print(f'\nSelected: {g}; other layers PMMA={H_PMMA_TOP:g}, ITO/DAST/ITO={H_ITO1:g}/{H_DAST:g}/{H_ITO2:g} nm. No Au.')
    print('PhaseTsolAbs is the largest absolute I-M difference at the SAME tested endpoint voltage; signed value and bias retained. Source side is air; reverse source is SiO2 half-space.')
    print('Transparency and reversal constraints are endpoint-only; plotted intermediate points are sampled, not continuous proofs. Both phases may reverse in different orientations.')
    print('Fixed phase datasets; no temperature dynamics or energy-saving prediction. Small discrete screen, not global optimum or converged performance. Same-height film/patch spacers have different SiO2 volumes.')
    print(f'Elapsed={time.perf_counter()-started:.1f}s; ONE figure; no files exported.')
    if SHOW_FIGURE:plt.show()
    return dict(Stage1Table=stage1,SeedTable=seeds,AllTable=table,ParetoTable=table[table.Pareto],Selected=chosen,DetailTable=detail,ComparisonTable=comparison,Figure=fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--project',type=Path);parser.add_argument('--no-show',action='store_true');parser.add_argument('--stage1-only',action='store_true');args=parser.parse_args()
    if args.project is not None:PROJECT=args.project.resolve()
    if args.no_show:SHOW_FIGURE=False
    if args.stage1_only:RUN_STAGE2=False
    Results=main();globals().update(Results)
