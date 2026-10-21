
#!/usr/bin/env python3
"""Local 2D TORCWA optimization around A/B for VO2/Au square Janus cavity.
Standalone, no MATLAB, no exports; place beside original optical data folders.
"""
from pathlib import Path
import time
TSTART=time.perf_counter()
import numpy as np,pandas as pd,torch,torcwa,matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from matplotlib.patches import Rectangle
# ===== PARAMETERS (0 CPU, 1 CUDA) =====
flag=0; P=400.; H_CAP=20.; H_ITO_TOP=5.; H_ITO_BOTTOM=5.
A={'WIDTH':200.,'H_VO2':18.,'H_AU':12.,'H_DAST':90.}
B={'WIDTH':200.,'H_VO2':21.,'H_AU':12.,'H_DAST':90.}
GRID={'WIDTH':(180.,190.,200.,210.,220.),'H_VO2':(15.,18.,21.,24.),'H_AU':(9.,10.5,12.,13.5),'H_DAST':(90.,92.5,95.,100.)}
PAIR_SCAN=1; TOP_K=6; SCREEN_STEP=40.; REFINE_STEP=20.; ORDER=2; NX=NY=40
VOLTAGES=(-8.5,0.,8.5); PHASES=('I','M')
T_MIN=.50; T_BUFFER=.003; J_TARGET=.04; EO_TARGET=.02; THERM_TARGET=.025
ITO_FILE='4TIO.xlsx'; DAST_SOURCE='hybrid'; EO_SIGN=1.; E_LIMIT=100.; SHOW_FIGURE=1
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
# ===== OBJECTIVES AND TWO-SEED LOCAL SEARCH =====
KEYS=('WIDTH','H_VO2','H_AU','H_DAST')
def gkey(g):return tuple(float(g[k]) for k in KEYS)
def label(g):return ' '.join(f'{k}={g[k]:g}' for k in KEYS)
def summary(tab):
    i=tab[tab.Phase=='I'].set_index('Voltage');m=tab[tab.Phase=='M'].set_index('Voltage')
    T=float(tab.Tlum.min());J=float(tab.DA.abs().min())
    EO=max(float(i.DA.max()-i.DA.min()),float(m.DA.max()-m.DA.min()))
    TH=max(abs(float(m.loc[v,'DA']-i.loc[v,'DA'])) for v in VOLTAGES)
    rev=any(float(x.DA.iloc[0]*x.DA.iloc[-1])<0 for x in (i,m))
    score=(.05*min(T/.55,1.5)+.20*min(J/J_TARGET,1.5)+.50*min(EO/EO_TARGET,1.5)
           +.25*min(TH/THERM_TARGET,1.5)+.12*min(max(T-T_MIN,0)/.02,1.))
    return dict(score=float(score),Tlum=T,J=J,EO=EO,Thermal=TH,Feasible=bool(T>=T_MIN-1e-12),Reversal=rev)
def rank(r):return (bool(r['Feasible']),r['score']+(min(r['Tlum']-T_MIN,T_BUFFER)*3 if r['Feasible'] else -15*(T_MIN-r['Tlum'])),r['Tlum'])
def report(tag,k,r):
    print(f'{tag:>8} #{k:02d} {"PASS" if r["Feasible"] else "FAIL"}  score={r["score"]:.4f}  T={r["Tlum"]:.4f}  J={r["J"]:.4f}  EO={r["EO"]:.4f}  TH={r["Thermal"]:.4f}  {r["seconds"]:.1f}s  {label(r["g"])}',flush=True)
def select_shortlist(cache):
    good=[r for r in cache.values() if r['Feasible']]; ordered=sorted(good or list(cache.values()),key=rank,reverse=True)
    chosen={gkey(A):A,gkey(B):B}
    for field in ('EO','Thermal','J'):
        r=max(good,key=lambda r:r[field]) if good else None
        if r is not None:chosen[gkey(r['g'])]=r['g']
    for r in ordered:
        if len(chosen)>=TOP_K:break
        chosen[gkey(r['g'])]=r['g']
    return list(chosen.values())[:TOP_K]
def make_figure(d,refA,refB,best,cache):
    w=d['w'];g=best['g'];fig,axes=plt.subplots(2,4,figsize=(19,9),layout='constrained');ax=axes.ravel();z=0;left=(P-g['WIDTH'])/2
    for nm,h,c,sq in [('PMMA',H_CAP,'#e3f1eb',False),('VO2',g['H_VO2'],'#c99472',True),('Au',g['H_AU'],'#e9ca76',True),('ITO',H_ITO_TOP,'#8dc8d7',False),('DAST',g['H_DAST'],'#d9b5e6',False),('ITO',H_ITO_BOTTOM,'#8dc8d7',False)]:
        if sq:ax[0].add_patch(Rectangle((0,z),P,h,fc='#e3f1eb',ec='.5',lw=.4))
        ax[0].add_patch(Rectangle((left if sq else 0,z),g['WIDTH'] if sq else P,h,fc=c,ec='.5',lw=.5))
        ax[0].text(P/2,z+h/2,nm,ha='center',va='center',fontsize=8);z+=h
    ax[0].set(xlim=(0,P),ylim=(z+15,-10),title='(a) Selected structure',xlabel='x (nm)',ylabel='z (nm)')
    rs=list(cache.values());ids=np.arange(1,len(rs)+1);ax[1].scatter(ids,[r['score'] for r in rs],c=['tab:green' if r['Feasible'] else 'gray' for r in rs],s=24)
    ax[1].set(title='(b) Local candidates (green: T≥0.50)',xlabel='Screened candidate',ylabel='Score')
    s=ax[2].scatter([r['Tlum'] for r in rs],[r['EO'] for r in rs],c=[r['Thermal'] for r in rs],cmap='viridis',s=36)
    ax[2].axvline(T_MIN,color='.5',ls='--');ax[2].scatter(best['Tlum'],best['EO'],marker='*',c='red',s=140,label='Refined choice')
    fig.colorbar(s,ax=ax[2],label='Thermal tuning');ax[2].set(title='(c) Transparency vs EO',xlabel='Minimum luminous T',ylabel='EO swing');ax[2].legend(fontsize=8)
    models=[('A',refA,':'),('B',refB,'--'),('Selected',best,'-')]
    for ph,c in [('I','tab:blue'),('M','tab:red')]:
        for tag,r,ls in models:
            q=r['spec'][(ph,0.)];ax[3].plot(w,q[:,0,2]-q[:,1,2],c=c,ls=ls,lw=1.15,label=f'{tag} {ph}')
            df=r['table'][r['table'].Phase==ph].sort_values('Voltage')
            ax[4].plot(df.Voltage,df.DA,c=c,ls=ls,marker='o',ms=3,label=f'{tag} {ph}')
            ax[5].plot(df.Voltage,df.Tlum,c=c,ls=ls,marker='o',ms=3,label=f'{tag} {ph}')
    ax[3].axhline(0,c='.5',lw=.7);ax[3].set(title='(d) Directional absorption at 0 V',xlabel='Wavelength (nm)',ylabel='AF - AB');ax[3].legend(fontsize=7,ncol=2)
    ax[4].axhline(0,c='.5',lw=.7);ax[4].set(title='(e) Electrical tuning',xlabel='DAST voltage (V)',ylabel='Solar DA');ax[4].legend(fontsize=7,ncol=2)
    ax[5].axhline(T_MIN,c='.4',ls='--');ax[5].set(title='(f) Luminous T (all states)',xlabel='DAST voltage (V)',ylabel='Min(TF,TB)',ylim=(.3,.65));ax[5].legend(fontsize=7,ncol=2)
    cols=('Tlum','J','EO','Thermal');targets=(T_MIN,J_TARGET,EO_TARGET,THERM_TARGET);positions=np.arange(4)
    for k,(tag,r,ls) in enumerate(models):ax[6].bar(positions+(k-1)*.25,[r[x]/v for x,v in zip(cols,targets)],width=.24,label=tag)
    ax[6].set_xticks(positions,['Tlum','Janus','EO','Thermal']);ax[6].axhline(1,color='.5',ls=':');ax[6].set(title='(g) Relative objectives',ylabel='Value / target');ax[6].legend(fontsize=8)
    for tag,r,ls in models:
        di=r['table'][r['table'].Phase=='I'].set_index('Voltage');dm=r['table'][r['table'].Phase=='M'].set_index('Voltage')
        ax[7].plot(VOLTAGES,[(dm.loc[v,'DA']-di.loc[v,'DA'])*100 for v in VOLTAGES],ls=ls,marker='o',ms=4,label=tag)
    ax[7].axhline(0,c='.5',lw=.7);ax[7].set(title='(h) VO2 thermal tuning',xlabel='DAST voltage (V)',ylabel='M - I solar DA (pp)');ax[7].legend(fontsize=8)
    for a in ax:a.grid(alpha=.15)
    fig.suptitle('A/B local search: 2D VO2/Au squares + ITO/DAST/ITO Janus (preliminary)')
    if not best['Feasible']:ax[0].text(.5,.93,'NO FEASIBLE DESIGN',ha='center',transform=ax[0].transAxes,color='red',weight='bold')
    return fig
def main():
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise RuntimeError('flag=1 requires working CUDA; use flag=0 for CPU')
    if P<=max(GRID['WIDTH']) or min(NX,NY)<4*ORDER+4 or ORDER<1:raise ValueError('Invalid period/grid/order settings')
    if set(PHASES)!={'I','M'} or tuple(VOLTAGES)!=(-8.5,0.,8.5):raise ValueError('This search expects two VO2 phases and -8.5,0,+8.5V')
    if min(GRID['H_DAST'])*E_LIMIT<1000*max(abs(v) for v in VOLTAGES):raise ValueError('DAST field limit exceeded')
    if TOP_K<2:raise ValueError('TOP_K must be at least 2 to include A and B')
    torch.set_num_threads(1);device=torch.device('cuda' if flag else 'cpu')
    if flag:torch.backends.cuda.matmul.allow_tf32=False
    print(f'Local A/B design; device={device}; N={ORDER} ({(2*ORDER+1)**2} modes), grid={NX}x{NY}; screening {SCREEN_STEP:g} nm, refinement {REFINE_STEP:g} nm')
    print(f'Hard constraint six-state min Tlum >= {T_MIN:.3f}; EO Emax={E_LIMIT:g} V/um; PMMA/[VO2/Au square]/ITO/DAST/ITO/SiO2')
    print('CAUTION: ITO file identity unverified; DAST scalar EO; N=2 is exploratory, not converged.')
    ds=load_data(SCREEN_STEP);dr=load_data(REFINE_STEP);maskcache={};cache={};invalid={}
    def compute(g,d):
        maskcache.setdefault(g['WIDTH'],square_mask(g['WIDTH'],NX,NY,device))
        return simulate(d,g,device,ORDER,maskcache[g['WIDTH']])
    validate(ds,device,ORDER,square_mask(A['WIDTH'],NX,NY,device))
    def screen(g):
        k=gkey(g)
        if k not in cache and k not in invalid:
            try:r=compute(g,ds)
            except FloatingPointError as e:invalid[k]=str(e);print(f'SKIP {label(g)}: {e}',flush=True);return None
            cache[k]=r;report('SCREEN',len(cache),r)
        return cache.get(k)
    print('\nSTAGE 1: A and B + one-parameter local neighbors')
    for seed in (A,B):
        screen(seed)
        for param,vals in GRID.items():
            for value in vals:
                g=seed.copy();g[param]=value;screen(g)
    if PAIR_SCAN:
        good=[r for r in cache.values() if r['Feasible']];anchor=max(good or list(cache.values()),key=rank)['g']
        print('\nSTAGE 2: one small WIDTH x AU interaction scan around best feasible screen result')
        ws=[v for v in GRID['WIDTH'] if abs(v-anchor['WIDTH'])<=10]
        aus=[v for v in GRID['H_AU'] if abs(v-anchor['H_AU'])<=1.5]
        for w in ws:
            for au in aus:
                g=anchor.copy();g.update(WIDTH=w,H_AU=au);screen(g)
    short=select_shortlist(cache)
    print(f'\nLOCAL SCREEN: {len(cache)} cases, {sum(r["Feasible"] for r in cache.values())} feasible, {len(invalid)} rejected. REFINE {len(short)} cases.')
    for i,g in enumerate(short,1):print(f'  Selected for refinement {i}: {label(g)}')
    fine=[]
    for i,g in enumerate(short,1):
        try:r=compute(g,dr)
        except FloatingPointError as e:print(f'REFINE SKIP {label(g)}: {e}',flush=True);continue
        fine.append(r);report('REFINE',i,r)
    if not fine:raise RuntimeError('No refinement succeeded')
    rA=next((r for r in fine if gkey(r['g'])==gkey(A)),None);rB=next((r for r in fine if gkey(r['g'])==gkey(B)),None)
    if rA is None or rB is None:raise RuntimeError('A/B reference missing after refinement')
    feasible=[r for r in fine if r['Feasible']]
    best=max(feasible,key=rank) if feasible else max(fine,key=rank)
    print('\n'+('SUCCESS: refined candidate meets Tlum>=0.50' if feasible else 'NO FEASIBLE REFINED CANDIDATE: showing a reference only'))
    print('\nREFINED COMPARISON (20 nm, consistent six states)')
    for tag,r in [('A',rA),('B',rB),('SELECTED',best)]:
        print(f'{tag}: {label(r["g"])}; T={r["Tlum"]:.6f}, J={r["J"]:.6f}, EO={r["EO"]:.6f}, Thermal={r["Thermal"]:.6f}, feasible={r["Feasible"]}')
        print(r['table'].to_string(index=False,float_format=lambda x:f'{x:.6f}'))
    print('Zero-order T reciprocity max =',max(float(np.max(np.abs(q[:,0,3]-q[:,1,3]))) for q in best['spec'].values()))
    print(f'COMPUTATION TIME = {time.perf_counter()-TSTART:.1f} s, coarse={len(cache)}, fine={len(fine)}; no files exported.')
    fig=make_figure(dr,rA,rB,best,cache)
    if SHOW_FIGURE:plt.show()
    return best,rA,rB,pd.DataFrame([{**r['g'],**{k:r[k] for k in ('score','Tlum','J','EO','Thermal','Feasible')}} for r in cache.values()]),fig
if __name__=='__main__':Best,ReferenceA,ReferenceB,ScreenTable,Figure=main()
