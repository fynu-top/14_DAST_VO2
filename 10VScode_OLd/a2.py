#!/usr/bin/env python3
"""Constrained 2D TORCWA search of PMMA/[VO2/Au squares]/ITO/DAST/ITO/SiO2 Janus.
One figure, no simulation exports, no MATLAB. This is exploratory, NOT a convergence study.
Place beside optical_constants, Atmospheric data, VO2_temperature_data.
"""
from pathlib import Path
import time
TSTART=time.perf_counter()
import numpy as np,pandas as pd,torch,torcwa,matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from matplotlib.patches import Rectangle
# ===== PARAMETERS: EDIT HERE =====
flag=0                            # 0 CPU, 1 CUDA
P=400.;H_CAP=20.;H_ITO_TOP=5.;H_ITO_BOTTOM=5.
BASE={'WIDTH':180.,'H_VO2':12.,'H_AU':12.,'H_DAST':95.}   # original transparent baseline
PRIOR={'WIDTH':260.,'H_VO2':18.,'H_AU':18.,'H_DAST':90.}  # previous EO-strong design
SEARCH={'WIDTH':(200.,220.,240.,260.),'H_VO2':(12.,15.,18.,21.),'H_AU':(9.,12.,15.,18.),'H_DAST':(90.,95.,100.,110.)}
SEED_K=2;TOP_K=4;LOCAL_POLISH=0  # 0 faster; 1 extra coordinate search
SCREEN_STEP=40.;REFINE_STEP=20.;ORDER=2;NX=NY=40
VOLTAGES=(-8.5,0.,8.5);PHASES=('I','M')
T_MIN=.50;J_TARGET=.06;EO_TARGET=.025;THERM_TARGET=.025
ITO_FILE='4TIO.xlsx'            # Check experimentally: file name does NOT establish ITO identity
DAST_SOURCE='hybrid';EO_SIGN=1.;E_LIMIT=100.;SHOW_FIGURE=1
PROJECT=Path(__file__).resolve().parent;DTYPE=torch.complex128
TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
def read_table(path,w,cols=(0,1,2),unit='nm',outside_zero=False):
    if not path.is_file():raise FileNotFoundError(f'Missing {path}; put this .py beside the original data folders')
    f=pd.read_csv(path,header=None,comment='#',on_bad_lines='skip') if path.suffix.lower()=='.csv' else pd.read_excel(path,header=None,engine='openpyxl')
    if f.shape[1]<=max(cols):raise ValueError(f'{path.name}: missing columns {cols}')
    a=f.iloc[:,list(cols)].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    a=a[np.isfinite(a).all(axis=1)];a=a[np.argsort(a[:,0])]
    if len(a)<3:raise ValueError(f'{path.name}: fewer than 3 numeric rows')
    a=a[np.unique(a[:,0],return_index=True)[1]];x=a[:,0]*(1000 if unit=='um' else 1)
    if not outside_zero and (x[0]>w[0]+1e-5 or x[-1]<w[-1]-1e-5):raise ValueError(f'{path.name}: wavelength range {x[0]:g}–{x[-1]:g} nm cannot cover {w[0]:g}–{w[-1]:g} nm')
    return np.column_stack([np.interp(w,x,a[:,k],left=0 if outside_zero else None,right=0 if outside_zero else None) for k in range(1,a.shape[1])])
def load_data(step):
    w=np.unique(np.r_[np.arange(400.,2500.0001,step),2500.,457.9,488.,514.5]);m=PROJECT/'optical_constants';at=PROJECT/'Atmospheric data';vdir=PROJECT/'VO2_temperature_data'
    d={'w':w};source={'PMMA':'9PMMA.xlsx','VO2_I':'6VO2_I.xlsx','VO2_M':'7VO2_M.xlsx','SiO2':'2SiO2.xlsx','Au':'11Au.xlsx','ITO':ITO_FILE}
    for key,file in source.items():
        a=read_table(m/file,w,unit='um');n,k=a[:,0],a[:,1]
        if n.min()<=0 or k.min()<-1e-8:raise ValueError(f'{file}: invalid n/k')
        if key in ('PMMA','SiO2') and np.max(k)>1e-3:print(f'NOTICE: {file} loss omitted to define lossless background/substrate')
        d[key]=n+1j*(np.maximum(k,0) if key not in ('PMMA','SiO2') else 0)
    if DAST_SOURCE=='xlsx':
        a=read_table(m/'10DAST.xlsx',w,unit='um');n=a[:,0];k=np.maximum(a[:,1],0)
    elif DAST_SOURCE=='hybrid':
        n=np.full(len(w),2.2);ix=w>650.;x=w[ix];ns=np.sqrt(2.078**2+1.645*533.**2/(x*x-533.**2));blend=.5-.5*np.cos(np.pi*np.clip((x-650)/100,0,1));n[ix]=(1-blend)*2.2+blend*ns
        alpha=np.full(len(w),1.9e5);ix=(w>430)&(w<=514.5)
        alpha[ix]=10**PchipInterpolator([430,457.9,488,514.5],np.log10([1.9e5,4.1e5,4.9e5,4.79e5]))(w[ix])
        ix=(w>514.5)&(w<=900);alpha[ix]=10**np.interp(w[ix],[514.5,700,750,900],np.log10([4.79e5,5.,1.8,.8]))
        ix=(w>900)&(w<=2100);x=w[ix];alpha[ix]=.7+np.exp(-.5*((x-1100)/55)**2)+1.5*np.exp(-.5*((x-1400)/70)**2)+3*np.exp(-.5*((x-1700)/95)**2)
        alpha[w>2100]=.7;k=alpha*w*1e-7/(4*np.pi)
    else:raise ValueError("DAST_SOURCE must be 'hybrid' or 'xlsx'")
    if n.min()<=0 or k.min()<0:raise ValueError('DAST n/k invalid')
    r=PchipInterpolator([720.,750.,1550.],[530.,445.,200.])(np.clip(w,720.,1550.))
    d['DASTn']=n;d['DASTk']=k;d['beta']=EO_SIGN*.5*n**3*r*1e-3 # delta n = beta*V/H_DAST(nm)
    d['solar']=read_table(vdir/'ASTMG173.xlsx',w,(0,2),outside_zero=True)[:,0]
    vp=at/'CIE_Vlambda.xlsx';vp=vp if vp.is_file() else vdir/'CIE_Vlambda.csv'
    if not vp.is_file():raise FileNotFoundError('CIE_Vlambda.xlsx or CIE_Vlambda.csv missing')
    vf=pd.read_excel(vp,header=None,engine='openpyxl') if vp.suffix.lower()=='.xlsx' else pd.read_csv(vp,header=None)
    cv=2 if vf.shape[1]>=4 else 1;vv=read_table(vp,w,(0,cv),outside_zero=True)[:,0]
    dp=at/'CIE_std_illum_D65.csv';dp=dp if dp.is_file() else at/'CIE_std_illum_D65.xlsx'
    d['lum']=vv*read_table(dp,w,(0,1),outside_zero=True)[:,0]
    if TRAPZ(d['solar'],w)<=0 or TRAPZ(d['lum'],w)<=0:raise ValueError('Invalid ASTM/D65 luminous weight tables')
    return d
# ===== VECTORIAL 2D TORCWA; PROPAGATING ORDERS ONLY =====
def square_mask(width,nx,ny,device):
    def fractional(n):
        x=np.arange(n)/n;lo=(1-width/P)/2;hi=(1+width/P)/2
        return np.clip(np.minimum(x+1/n,hi)-np.maximum(x,lo),0,None)*n
    return torch.as_tensor(fractional(nx)[:,None]*fractional(ny)[None,:],dtype=DTYPE,device=device)
def port_power(s,lam,nsub,order,direction,port,device):
    nout=(1. if direction=='forward' else nsub) if port=='reflection' else (nsub if direction=='forward' else 1.)
    ms=range(-order,order+1);orders=np.asarray([[i,j] for i in ms for j in ms if (i*i+j*j)*(lam/P)**2<nout*nout-1e-9],dtype=np.int64).reshape(-1,2)
    a=torch.zeros(len(orders),device=device,dtype=torch.float64)
    for pol in ('pp','sp'):
        z=s.S_parameters(orders=orders,direction=direction,port=port,polarization=pol,ref_order=[0,0],power_norm=True)
        a+=z.abs().square().reshape(-1).real
    b=a.detach().cpu().numpy()
    if not np.isfinite(b).all():raise FloatingPointError(f'Nonfinite power at {lam:g} nm')
    return float(b.sum()),float(b[np.flatnonzero((orders==[0,0]).all(axis=1))[0]])
def solve_one(d,j,phase,voltage,device,order,mask,g):
    w=float(d['w'][j]);lam=w+1e-5 if abs(w-P)<1e-7 else w
    nsub=float(d['SiO2'][j].real);pm=complex(d['PMMA'][j]**2)
    nd=d['DASTn'][j]+d['beta'][j]*voltage/g['H_DAST']+1j*d['DASTk'][j]
    if nd.real<=0:raise ValueError('Nonphysical DAST refractive index')
    s=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    s.add_output_layer(eps=nsub*nsub);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
    s.add_layer(thickness=H_CAP,eps=pm)
    for key,h in [('VO2_'+phase,g['H_VO2']),('Au',g['H_AU'])]:
        n=d[key][j];s.add_layer(thickness=h,eps=pm+mask*(complex(n*n)-pm))
    for n,h in [(d['ITO'][j],H_ITO_TOP),(nd,g['H_DAST']),(d['ITO'][j],H_ITO_BOTTOM)]:
        s.add_layer(thickness=h,eps=complex(n*n))
    s.solve_global_smatrix();q=np.zeros((2,4))
    for i,direction in enumerate(('forward','backward')):
        q[i,0],_=port_power(s,lam,nsub,order,direction,'reflection',device)
        q[i,1],q[i,3]=port_power(s,lam,nsub,order,direction,'transmission',device)
        q[i,2]=1-q[i,0]-q[i,1]
    if not np.isfinite(q).all() or q.min()<-3e-4 or q.max()>1+3e-4:raise FloatingPointError(f'Passivity failed at {w:g}nm, {phase}, {voltage:+g}V: {q}')
    return q
def simulate(d,g,device,order,mask):
    t=time.perf_counter();spec={};rows=[]
    with torch.no_grad():
        for ph in PHASES:
            for v in VOLTAGES:
                q=np.array([solve_one(d,j,ph,v,device,order,mask,g) for j in range(len(d['w']))]);spec[(ph,v)]=q
                rows.append(dict(Phase=ph,Voltage=v,**metrics(d,q)))
    if device.type=='cuda':torch.cuda.synchronize()
    tab=pd.DataFrame(rows);r=summary(tab);return dict(g=g.copy(),table=tab,spec=spec,**r,seconds=time.perf_counter()-t)
def metrics(d,q):
    w=d['w'];solar=d['solar'];lum=d['lum'];vis=(w>=400)&(w<=780)
    def weighted(y,weight,ix):return float(TRAPZ(y[ix]*weight[ix],w[ix])/TRAPZ(weight[ix],w[ix]))
    af=weighted(q[:,0,2],solar,slice(None));ab=weighted(q[:,1,2],solar,slice(None))
    tf=weighted(q[:,0,1],lum,vis);tb=weighted(q[:,1,1],lum,vis)
    return {'DA':af-ab,'AsF':af,'AsB':ab,'TlumF':tf,'TlumB':tb,'Tlum':min(tf,tb),
            'TsolF':weighted(q[:,0,1],solar,slice(None)),'TsolB':weighted(q[:,1,1],solar,slice(None))}
def validate(d,device,order,mask):
    worst=0.;loss=0.
    with torch.no_grad():
        for wl in (405.,415.,900.,1550.):
            n=float(np.interp(wl,d['w'],d['SiO2'].real));s=torcwa.rcwa(freq=1/wl,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
            s.add_output_layer(eps=n*n);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
            s.add_layer(thickness=12.,eps=1.);s.solve_global_smatrix()
            re=((1-n)/(1+n))**2;te=4*n/(1+n)**2
            for side in ('forward','backward'):
                r,_=port_power(s,wl,n,order,side,'reflection',device);t,_=port_power(s,wl,n,order,side,'transmission',device)
                worst=max(worst,abs(r-re),abs(t-te),abs(r+t-1))
            z=torcwa.rcwa(freq=1/wl,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
            z.add_output_layer(eps=n*n);z.set_incident_angle(inc_ang=0.,azi_ang=0.)
            z.add_layer(thickness=60.,eps=1.49**2+mask*(1.7**2-1.49**2));z.solve_global_smatrix()
            for side in ('forward','backward'):
                r,_=port_power(z,wl,n,order,side,'reflection',device);t,_=port_power(z,wl,n,order,side,'transmission',device);loss=max(loss,abs(r+t-1))
    print(f'Fresnel max error={worst:.3e}; lossless 2D max |R+T-1|={loss:.3e}')
    if worst>1e-6 or loss>1e-4:raise RuntimeError('2D solver checks failed; abort device calculation')
# ===== TRANSPARENCY-CONSTRAINED MULTI-OBJECTIVE OPTIMIZER =====
def summary(tab):
    i=tab[tab.Phase=='I'].set_index('Voltage');m=tab[tab.Phase=='M'].set_index('Voltage')
    J=float(tab.DA.abs().min());eo=max(float(i.DA.max()-i.DA.min()),float(m.DA.max()-m.DA.min()))
    therm=max(abs(float(m.loc[v,'DA']-i.loc[v,'DA'])) for v in VOLTAGES)
    T=float(tab.Tlum.min());rev=any(float(x.DA.iloc[0]*x.DA.iloc[-1])<0 for x in (i,m))
    quality=.10*min(T/.60,1.3)+.30*min(J/J_TARGET,1.5)+.35*min(eo/EO_TARGET,1.5)+.25*min(therm/THERM_TARGET,1.5)
    return dict(score=float(quality),Tlum=T,J=J,EO=eo,Thermal=therm,Feasible=bool(T>=T_MIN),Reversal=rev)
def gkey(g):return tuple(float(g[k]) for k in BASE)
def label(g):return ' '.join(f'{k}={g[k]:g}' for k in BASE)
def order_key(r):return (int(r['Feasible']),r['score'] if r['Feasible'] else r['score']-12*max(0.,T_MIN-r['Tlum']),r['Tlum'])
def print_rec(stage,n,r):
    ok='PASS' if r['Feasible'] else 'FAIL'
    print(f'{stage:>9} #{n:02d} {ok} score={r["score"]:.4f} T={r["Tlum"]:.4f} J={r["J"]:.4f} EO={r["EO"]:.4f} TH={r["Thermal"]:.4f} reversal={r["Reversal"]} {r["seconds"]:.1f}s  {label(r["g"])}',flush=True)
def plot_figure(d,base,old,best,screen):
    g=best['g'];w=d['w'];fig,ax=plt.subplots(2,4,figsize=(18.8,9.0),layout='constrained');a=ax.ravel();z=0;left=(P-g['WIDTH'])/2
    for nm,h,c,pattern in [('PMMA',H_CAP,'#e3f1eb',False),('VO2',g['H_VO2'],'#c99472',True),('Au',g['H_AU'],'#e9ca76',True),('ITO',H_ITO_TOP,'#8dc8d7',False),('DAST',g['H_DAST'],'#d9b5e6',False),('ITO',H_ITO_BOTTOM,'#8dc8d7',False)]:
        if pattern:a[0].add_patch(Rectangle((0,z),P,h,fc='#e3f1eb',ec='.6',lw=.4))
        a[0].add_patch(Rectangle((left if pattern else 0,z),g['WIDTH'] if pattern else P,h,fc=c,ec='.5',lw=.5))
        a[0].text(P/2,z+h/2,nm,ha='center',va='center',fontsize=8);z+=h
    a[0].set(xlim=(0,P),ylim=(z+14,-12),xlabel='x (nm)',ylabel='z (nm)',title='(a) Selected geometry')
    allr=list(screen.values());x=np.arange(1,len(allr)+1);good=[r['Feasible'] for r in allr]
    a[1].scatter(x,[r['score'] for r in allr],c=['tab:green' if v else 'tab:gray' for v in good],s=24)
    a[1].set(title='(b) Screening objective (green: T≥0.50)',xlabel='Candidate index',ylabel='Objective for ranking')
    sc=a[2].scatter([r['Tlum'] for r in allr],[r['J'] for r in allr],c=[r['EO'] for r in allr],cmap='viridis',s=32)
    a[2].axvline(T_MIN,ls='--',c='gray',lw=.9);a[2].scatter(best['Tlum'],best['J'],s=170,c='red',marker='*',label='Selected (refined)')
    fig.colorbar(sc,ax=a[2],label='EO tuning');a[2].set(title='(c) Transparency / Janus trade-off',xlabel='Minimum Tlum',ylabel='Minimum |solar ΔA|');a[2].legend(fontsize=8)
    models=(('Initial',base,'--'),('Previous',old,':'),('Selected',best,'-'))
    for ph,c in [('I','tab:blue'),('M','tab:red')]:
        for nm,r,ls in models:
            q=r['spec'][(ph,0.)];a[3].plot(w,q[:,0,2]-q[:,1,2],ls=ls,c=c,label=f'{nm} {ph}',lw=1.25)
            df=r['table'][r['table'].Phase==ph].sort_values('Voltage')
            a[4].plot(df.Voltage,df.DA,ls=ls,marker='o',ms=3,c=c,label=f'{nm} {ph}')
            a[5].plot(df.Voltage,df.Tlum,ls=ls,marker='o',ms=3,c=c,label=f'{nm} {ph}')
    a[3].axhline(0,color='.55',lw=.7);a[3].set(title='(d) Solar directional spectrum at 0 V',xlabel='Wavelength (nm)',ylabel='AF - AB');a[3].legend(fontsize=7,ncol=2)
    a[4].axhline(0,color='.55',lw=.7);a[4].set(title='(e) Electric tuning',xlabel='DAST voltage (V)',ylabel='Solar ΔA');a[4].legend(fontsize=7,ncol=2)
    a[5].axhline(T_MIN,ls='--',c='.4',lw=1,label='Tlum ≥ 0.50');a[5].set(title='(f) Luminous T (6 states)',xlabel='DAST voltage (V)',ylabel='Min(TF,TB)',ylim=(0,1));a[5].legend(fontsize=7,ncol=2)
    keys=['Tlum','J','EO','Thermal'];den=[T_MIN,J_TARGET,EO_TARGET,THERM_TARGET];z=np.arange(4)
    for k,(nm,r,ls) in enumerate(models):a[6].bar(z+(k-1)*.25,[r[t]/v for t,v in zip(keys,den)],width=.24,label=nm)
    a[6].set_xticks(z,['Tlum','Janus','EO','Thermal']);a[6].axhline(1,c='.5',ls=':');a[6].set(title='(g) Normalized performance',ylabel='Metric / exploratory target');a[6].legend(fontsize=8)
    for ph,c in [('I','tab:blue'),('M','tab:red')]:
        for nm,r,ls in models:a[7].plot(w,r['spec'][(ph,0.)][:,0,1],ls=ls,c=c,label=f'{nm} {ph}',lw=1.25)
    a[7].set(title='(h) Forward T spectrum at 0 V',xlabel='Wavelength (nm)',ylabel='Total T',ylim=(0,1));a[7].legend(fontsize=7,ncol=2)
    for v in a:v.grid(alpha=.15)
    fig.suptitle('Tlum-constrained 2D RCWA optimization: PMMA/[VO2/Au squares]/ITO/DAST/ITO/SiO2')
    return fig
def main():
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise RuntimeError('Check flag or CUDA installation')
    if ORDER<1 or min(NX,NY)<4*ORDER+4 or P<=max(SEARCH['WIDTH']):raise ValueError('Invalid Fourier order, real-space grid, or unit cell geometry')
    if set(PHASES)!={'I','M'} or len(VOLTAGES)!=3 or 0. not in VOLTAGES:raise ValueError('Require VO2 I/M and 3 voltages including zero')
    if SEED_K<1 or TOP_K<1 or min(SEARCH['H_DAST'])*E_LIMIT<max(abs(v) for v in VOLTAGES)*1000:raise ValueError('Seed count/DAST electric field exceeds limits')
    torch.set_num_threads(1);device=torch.device('cuda' if flag else 'cpu')
    if flag:torch.backends.cuda.matmul.allow_tf32=False
    print(f'\nDevice={device}; complex128; N={ORDER}, {NX}x{NY}; screen={SCREEN_STEP:g}nm; refine={REFINE_STEP:g}nm')
    print(f'Constrained search: min over all 6 states Tlum >= {T_MIN:.3f}, electric field <= {E_LIMIT:g} V/um')
    print('Geometry PMMA/[VO2/Au square]/ITO/DAST/ITO/SiO2; dimensions in nm')
    print('WARNING: ITO=4TIO.xlsx not verified; DAST=scalar hybrid approximation; exploratory non-converged RCWA')
    ds=load_data(SCREEN_STEP);dr=load_data(REFINE_STEP);validate(ds,device,ORDER,square_mask(BASE['WIDTH'],NX,NY,device))
    cache={};masks={};bad=set()
    def compute(g,data):
        wid=float(g['WIDTH']);masks.setdefault(wid,square_mask(wid,NX,NY,device))
        return simulate(data,g,device,ORDER,masks[wid])
    def screen(g):
        key=gkey(g)
        if key in bad:return None
        if key not in cache:
            try:r=compute(g,ds)
            except FloatingPointError as exc:print(f'REJECT due to passivity: {label(g)}: {exc}',flush=True);bad.add(key);return None
            cache[key]=r;print_rec('SCREEN',len(cache),r)
        return cache[key]
    screen(BASE);screen(PRIOR)
    # Stage 1: 4x4 scan of width/Au with reference VO2=12 nm, DAST=95 nm
    print('\nSTAGE 1: joint scan WIDTH x H_AU, fixed H_VO2=12, H_DAST=95')
    first=[]
    for width in SEARCH['WIDTH']:
        for au in SEARCH['H_AU']:
            g={'WIDTH':width,'H_VO2':12.,'H_AU':au,'H_DAST':95.};r=screen(g)
            if r is not None:first.append(r)
    if not first:raise RuntimeError('No acceptable initial candidates; inspect material data and Fourier truncation')
    seeds=[]
    for r in sorted(first,key=order_key,reverse=True):
        if r['g']['WIDTH'] not in [s['g']['WIDTH'] for s in seeds]:seeds.append(r)
        if len(seeds)>=SEED_K:break
    print('\nSTAGE 2 SEEDS:')
    for j,r in enumerate(seeds,1):print_rec('SEED',j,r)
    # Stage 2: VO2 x DAST joint scan near two promising WIDTH/Au combinations
    print('\nSTAGE 2: H_VO2 x H_DAST at selected WIDTH/Au seeds')
    for r in seeds:
        g0=r['g']
        for vo in SEARCH['H_VO2']:
            for da in SEARCH['H_DAST']:
                g={'WIDTH':g0['WIDTH'],'H_VO2':vo,'H_AU':g0['H_AU'],'H_DAST':da};screen(g)
    if LOCAL_POLISH:
        print('\nOPTIONAL POLISH: coordinate neighborhood about best screened design')
        for k in SEARCH:
            choices=[r for r in cache.values() if r['g']['WIDTH'] in SEARCH['WIDTH']]
            cur=max(choices,key=order_key)['g'].copy()
            for v in SEARCH[k]:
                g=cur.copy();g[k]=v;screen(g)
    ranked=sorted([r for r in cache.values() if r['g']['WIDTH'] in SEARCH['WIDTH']],key=order_key,reverse=True)
    print(f'\nSCREEN RESULTS: {len(cache)} valid, {sum(r["Feasible"] for r in ranked)} feasible within requested width range, {len(bad)} rejected')
    for j,r in enumerate(ranked[:min(10,len(ranked))],1):print_rec('TOP',j,r)
    # Refinement preserves diversity and a feasible reference; no unverified success claim
    candidates={gkey(BASE):BASE,gkey(PRIOR):PRIOR}
    for r in ranked[:TOP_K]:candidates[gkey(r['g'])]=r['g']
    print(f'\nREFINE: {len(candidates)} geometries at {REFINE_STEP:g} nm, full 6 states')
    refined=[]
    for j,g in enumerate(candidates.values(),1):
        try:r=compute(g,dr)
        except FloatingPointError as exc:print(f'REFINE REJECTED {label(g)}: {exc}');continue
        refined.append(r);print_rec('REFINE',j,r)
    base=next((r for r in refined if gkey(r['g'])==gkey(BASE)),None)
    prior=next((r for r in refined if gkey(r['g'])==gkey(PRIOR)),None)
    if base is None or prior is None:raise RuntimeError('Baseline or previous design failed refinement: cannot plot comparison')
    eligible=[r for r in refined if r['Feasible'] and r['g']['WIDTH'] in SEARCH['WIDTH']]
    if eligible:
        best=max(eligible,key=lambda r:r['score']);print('\nSUCCESS: Found an eligible Tlum >= 0.50 design in requested width range')
    else:
        best=max(refined,key=order_key)
        print('\nNO ELIGIBLE OPTIMUM in requested range at refined grid. Showing reference/nearest alternative, NOT a successful constrained optimum.')
    print('\nREFINED COMPARISON (same wavelength grid)')
    for nm,r in [('INITIAL',base),('PREVIOUS',prior),('SELECTED',best)]:
        print(f'{nm}: feasible={r["Feasible"]}; {label(r["g"])}; T={r["Tlum"]:.6f}, J={r["J"]:.6f}, EO={r["EO"]:.6f}, thermal={r["Thermal"]:.6f}, objective={r["score"]:.6f}')
        print(r['table'].to_string(index=False,float_format=lambda v:f'{v:.6f}'))
    print('Matched zero-order reciprocity max=',max(np.max(np.abs(q[:,0,3]-q[:,1,3])) for q in best['spec'].values()))
    print(f'Elapsed including imports/computation: {time.perf_counter()-TSTART:.1f}s; {len(cache)} screenings; {len(refined)} refinements; one figure, no exports')
    fig=plot_figure(dr,base,prior,best,cache)
    if SHOW_FIGURE:plt.show()
    return best,base,prior,pd.DataFrame([{**r['g'],**{k:r[k] for k in ('score','Tlum','J','EO','Thermal','Feasible')}} for r in cache.values()]),fig
if __name__=='__main__':Best,Baseline,Previous,ScreenTable,Figure=main()
