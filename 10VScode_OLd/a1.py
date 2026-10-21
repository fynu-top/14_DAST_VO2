
#!/usr/bin/env python3
"""Fast 2D TORCWA coordinate-search of PMMA/[VO2/Au squares]/ITO/DAST/ITO/SiO2 Janus.
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
flag=0                          # 0 CPU / 1 CUDA
P=400.;H_CAP=20.;H_ITO_TOP=5.;H_ITO_BOTTOM=5.
BASE={'WIDTH':180.,'H_VO2':12.,'H_AU':12.,'H_DAST':95.}
SEARCH={'WIDTH':(120.,160.,180.,220.,260.),'H_VO2':(6.,12.,18.),'H_AU':(6.,12.,18.),'H_DAST':(90.,95.,120.,150.)}
PASSES=1;TOP_K=3               # set PASSES=2 for second coordinate pass
SCREEN_STEP=40.;REFINE_STEP=20.;ORDER=2;NX=NY=40
VOLTAGES=(-8.5,0.,8.5);PHASES=('I','M')
T_TARGET=.60;J_TARGET=.12;EO_TARGET=.02;THERM_TARGET=.03
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
# ===== MULTI-OBJECTIVE SCREENING (NOT A GLOBAL OPTIMUM) =====
def summary(tab):
    i=tab[tab.Phase=='I'].set_index('Voltage');m=tab[tab.Phase=='M'].set_index('Voltage')
    J=float(tab.DA.abs().min());eo=max(float(i.DA.max()-i.DA.min()),float(m.DA.max()-m.DA.min()))
    therm=max(abs(float(m.loc[v,'DA']-i.loc[v,'DA'])) for v in VOLTAGES)
    T=float(tab.Tlum.min());sign=any(float(x.DA.iloc[0]*x.DA.iloc[-1])<0 for x in (i,m))
    score=.30*min(T/T_TARGET,1)+.28*min(J/J_TARGET,1)+.20*min(eo/EO_TARGET,1)+.22*min(therm/THERM_TARGET,1)
    score-=.25*max(0,.50-T)/.50+.20*max(0,.025-J)/.025
    if sign and J>=.005:score+=.05
    return dict(score=float(score),J=J,Tlum=T,EO=eo,Thermal=therm,Reversal=sign)
def gkey(g):return tuple(float(g[k]) for k in BASE)
def label(g):return ' '.join(f'{k}={g[k]:g}' for k in BASE)
def square_geom(g,device):return square_mask(g['WIDTH'],NX,NY,device)
def print_rec(stage,n,r):
    print(f'{stage:>8} #{n:02d}  score={r["score"]:+.4f} T={r["Tlum"]:.4f} J={r["J"]:.4f} EO={r["EO"]:.4f} TH={r["Thermal"]:.4f} rev={r["Reversal"]} {r["seconds"]:.1f}s  {label(r["g"])}',flush=True)
def draw(d,baseline,best,screen):
    g=best['g'];w=d['w'];fig,axs=plt.subplots(2,4,figsize=(19,9.2),layout='constrained');a=axs.ravel();z=0.;left=(P-g['WIDTH'])/2
    for nm,h,co,pattern in [('PMMA',H_CAP,'#dcece8',0),('VO2',g['H_VO2'],'#c99472',1),('Au',g['H_AU'],'#efce75',1),('ITO',H_ITO_TOP,'#8dc8d7',0),('DAST',g['H_DAST'],'#d6afe5',0),('ITO',H_ITO_BOTTOM,'#8dc8d7',0)]:
        if pattern:a[0].add_patch(Rectangle((0,z),P,h,fc='#dcece8',ec='.5',lw=.4))
        a[0].add_patch(Rectangle((left if pattern else 0,z),g['WIDTH'] if pattern else P,h,fc=co,ec='.5',lw=.5))
        a[0].text(P/2,z+h/2,nm,ha='center',va='center',fontsize=8);z+=h
    a[0].set(xlim=(0,P),ylim=(z+14,-12),xlabel='x (nm)',ylabel='z (nm)',title='(a) Selected square-array cross-section')
    s=list(screen.values());a[1].plot(range(1,len(s)+1),[v['score'] for v in s],'.-',lw=1);a[1].axhline(screen[gkey(BASE)]['score'],color='.5',ls=':',label='Baseline screening score')
    a[1].set(title='(b) Screening objective by evaluated candidate',xlabel='Candidate index',ylabel='Composite score');a[1].legend(fontsize=8)
    sc=a[2].scatter([v['Tlum'] for v in s],[v['J'] for v in s],c=[v['Thermal'] for v in s],cmap='viridis',s=40)
    a[2].scatter([screen[gkey(BASE)]['Tlum']],[screen[gkey(BASE)]['J']],facecolors='none',edgecolors='red',s=110,label='Baseline screening');a[2].scatter([best['Tlum']],[best['J']],marker='*',s=160,c='red',label='Selected refined')
    fig.colorbar(sc,ax=a[2],label='Thermal modulation');a[2].set(title='(c) Screening trade-off',xlabel='Min luminous T',ylabel='Min |signed solar ΔA|');a[2].legend(fontsize=8)
    for ph,cl in [('I','tab:blue'),('M','tab:red')]:
        for r,ls,nm in [(baseline,'--','Base'),(best,'-','Best')]:
            q=r['spec'][(ph,0.)];a[3].plot(w,q[:,0,2]-q[:,1,2],ls=ls,c=cl,label=f'{nm} {ph}')
            df=r['table'][r['table'].Phase==ph].sort_values('Voltage')
            a[4].plot(df.Voltage,df.DA,ls+'o',c=cl,label=f'{nm} {ph}')
            a[5].plot(df.Voltage,df.Tlum,ls+'o',c=cl,label=f'{nm} {ph}')
    a[3].axhline(0,color='.5',lw=.7);a[3].set(title='(d) Directional spectrum at 0 V',xlabel='Wavelength (nm)',ylabel='AF-AB');a[3].legend(fontsize=8)
    a[4].axhline(0,color='.5',lw=.7);a[4].set(title='(e) Electrical tuning',xlabel='Voltage (V)',ylabel='Signed solar ΔA');a[4].legend(fontsize=8)
    a[5].axhline(.6,color='.5',ls=':',lw=.7);a[5].set(title='(f) Min luminous T',xlabel='Voltage (V)',ylabel='Tlum',ylim=(0,1));a[5].legend(fontsize=8)
    names=['Tlum','Janus','EO','Thermal'];den=[T_TARGET,J_TARGET,EO_TARGET,THERM_TARGET];keys=['Tlum','J','EO','Thermal'];x=np.arange(4)
    a[6].bar(x-.18,[baseline[k]/v for k,v in zip(keys,den)],width=.36,label='Baseline');a[6].bar(x+.18,[best[k]/v for k,v in zip(keys,den)],width=.36,label='Selected')
    a[6].axhline(1,color='.5',ls=':',lw=.7);a[6].set_xticks(x,names);a[6].set(title='(g) Metrics / exploration targets',ylabel='Ratio to target');a[6].legend(fontsize=8)
    for ph,cl in [('I','tab:blue'),('M','tab:red')]:
        for r,ls,nm in [(baseline,'--','Base'),(best,'-','Best')]:
            a[7].plot(w,r['spec'][(ph,0.)][:,0,1],ls=ls,c=cl,label=f'{nm} {ph}')
    a[7].set(title='(h) Forward transmission at 0 V',xlabel='Wavelength (nm)',ylabel='Total T',ylim=(0,1));a[7].legend(fontsize=8)
    for ax in a:ax.grid(alpha=.14)
    fig.suptitle('Exploratory optimization: PMMA/[VO2/Au squares]/ITO/DAST/ITO/SiO2 (2D RCWA)',fontsize=13)
    return fig
def main():
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise RuntimeError('Check flag or CUDA availability')
    if ORDER<1 or min(NX,NY)<4*ORDER+4 or P<=max(SEARCH['WIDTH']):raise ValueError('Fourier/grid/period settings invalid')
    if 0. not in VOLTAGES or len(VOLTAGES)!=3 or set(PHASES)!={'I','M'}:raise ValueError('Require 3 voltages incl. 0 V, and phases I/M')
    for k,v in SEARCH.items():
        if min(v)<=0:raise ValueError(f'Nonpositive {k}')
    if min(SEARCH['H_DAST'])*E_LIMIT<max(abs(x) for x in VOLTAGES)*1000:raise ValueError('At least one DAST thickness exceeds field limit')
    torch.set_num_threads(1);device=torch.device('cuda' if flag else 'cpu')
    if flag:torch.backends.cuda.matmul.allow_tf32=False
    print(f'\nDevice={device}, complex128, TORCWA order={ORDER}, grid={NX}x{NY}; screen={SCREEN_STEP:g}nm, refine={REFINE_STEP:g}nm; passes={PASSES}')
    print('Geometry PMMA/[VO2/Au square]/ITO/DAST/ITO/SiO2; no internal SiO2 spacer, no MATLAB')
    print('WARNING: 4TIO.xlsx ITO identity is UNVERIFIED; PMMA/SiO2 k set to zero; DAST scalar hybrid EO is exploratory')
    ds=load_data(SCREEN_STEP);dr=load_data(REFINE_STEP);validate(ds,device,ORDER,square_geom(BASE,device))
    cache={};mask_cache={};history=[]
    def evaluate(g,data,stage):
        key=gkey(g);mkey=float(g['WIDTH'])
        if mkey not in mask_cache:mask_cache[mkey]=square_geom(g,device)
        return simulate(data,g,device,ORDER,mask_cache[mkey])
    def screen(g):
        key=gkey(g)
        if key not in cache:
            try:r=evaluate(g,ds,'screen')
            except FloatingPointError as e:
                print(f'REJECT: {label(g)}; {e}',flush=True);return None
            cache[key]=r;history.append(key);print_rec('SCREEN',len(history),r)
        return cache[key]
    cur=BASE.copy();initial=screen(cur)
    if initial is None:raise RuntimeError('Baseline passivity failed; increase ORDER and inspect materials before optimizing')
    for p in range(PASSES):
        for param,values in SEARCH.items():
            group=[screen(cur)]
            for value in values:
                trial=cur.copy();trial[param]=float(value);group.append(screen(trial))
            group=[x for x in group if x is not None];cur=max(group,key=lambda x:x['score'])['g'].copy()
            print(f'PASS {p+1} {param}: selected {label(cur)}',flush=True)
    print('\nTOP SCREEN CANDIDATES')
    ranked=sorted(cache.values(),key=lambda x:x['score'],reverse=True)
    for n,r in enumerate(ranked[:min(8,len(ranked))],1):print_rec('TOP',n,r)
    pick={gkey(BASE):BASE.copy()}
    for r in ranked[:TOP_K]:pick[gkey(r['g'])]=r['g'].copy()
    print(f'\nREFINEMENT: {len(pick)} candidate(s) at {REFINE_STEP:g} nm; both phases, all three voltages')
    refined=[]
    for n,g in enumerate(pick.values(),1):
        try:r=evaluate(g,dr,'refine')
        except FloatingPointError as e:print(f'REFINEMENT REJECTED: {label(g)}; {e}');continue
        refined.append(r);print_rec('REFINE',n,r)
    base=next((r for r in refined if gkey(r['g'])==gkey(BASE)),None)
    if base is None:raise RuntimeError('Refined baseline unavailable')
    best=max(refined,key=lambda x:x['score'])
    print('\nREFINED COMPARISON (all metrics computed on same wavelength grid)')
    for name,r in [('BASE',base),('SELECTED',best)]:
        print(f'{name}: {label(r["g"])}\n'+r['table'].to_string(index=False,float_format=lambda z:f'{z:.6f}'))
    print(f'\nSelected vs baseline: ΔTlum={best["Tlum"]-base["Tlum"]:+.5f}, ΔJ={best["J"]-base["J"]:+.5f}, ΔEO={best["EO"]-base["EO"]:+.5f}, ΔThermal={best["Thermal"]-base["Thermal"]:+.5f}')
    print('Zero-order transmission reciprocity max:',max(np.max(abs(q[:,0,3]-q[:,1,3])) for q in best['spec'].values()))
    print(f'SCREEN candidates={len(cache)}, refined={len(refined)}, wall time={time.perf_counter()-TSTART:.1f}s ({(time.perf_counter()-TSTART)/60:.2f}min); no exported simulation files')
    fig=draw(dr,base,best,cache)
    if SHOW_FIGURE:plt.show()
    return best,base,pd.DataFrame([{**r['g'],**{k:r[k] for k in ('score','J','Tlum','EO','Thermal')}} for r in cache.values()]),fig
if __name__=='__main__':Best,Baseline,ScreenTable,Figure=main()



