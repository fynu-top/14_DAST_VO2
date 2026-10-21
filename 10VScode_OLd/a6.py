
#!/usr/bin/env python3
"""Two-candidate limited numerical audit: true 2D aligned VO2/SiO2/Au square stack in PMMA, above continuous ITO/DAST/ITO.
Air / PMMA cap / patterned VO2 / patterned SiO2 / patterned Au / ITO / DAST /
ITO / semi-infinite lossless SiO2. Pattern gaps contain PMMA at EVERY patterned z.
Not the old 10.py geometry. Not yet optimized. No MATLAB, no automatic downloads,
no exports. Tested API: torcwa 0.1.4.2. Python 3.9-compatible; CPU default.
Dependencies: numpy scipy pandas openpyxl matplotlib torch torcwa.
Place beside original optical_constants, Atmospheric data, VO2_temperature_data.
"""
from pathlib import Path
import argparse,time,warnings
import numpy as np,pandas as pd,torch,torcwa
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
# ===== EDIT THESE SETTINGS; VS Code: Run Python File =====
flag=0                              # 0 CPU, 1 CUDA; small complex128 systems often favour CPU
P=400.; WIDTH=180.                   # common aligned square width; area fill=(WIDTH/P)^2
H_VO2=12.; H_SIO2_PATCH=60.; H_AU=12.; H_PMMA_TOP=15.
H_ITO1=5.; H_DAST=95.; H_ITO2=5.      # all three layers continuous, no bottom Au sheet
VOLTAGES=(-8.5,0.,8.5); PHASES=('I','M'); E_LIMIT=100.; MARGIN=.005
EO_SIGN=1.; DAST_SOURCE='hybrid'; ITO_FILE='4TIO.xlsx'
LOSSLESS_PMMA=True                   # same transparent-host assumption as 10.py; set False to retain k
RUN_OLD10=0
SHOW_FIGURE=1; PROJECT=Path(__file__).resolve().parent
DTYPE=torch.complex128; TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz

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
                     'PMMA':'9PMMA.xlsx','Au':'11Au.xlsx'}
    if RUN_OLD10: names['Ti']='14Ti.xlsx'
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


def solve_one(d,j,phase,bias,device,order,mask,geometry='pattern',eo=True):
    w=float(d['w'][j]);nsub=float(d['SiO2'][j].real);lam=safe_lambda(w,nsub,order)
    pm=complex(d['PMMA'][j]**2);nv=d['VI'][j] if phase=='I' else d['VM'][j]
    hd=95. if geometry=='old10' else H_DAST;nd=d['DASTn'][j]+d['beta'][j]*(bias if eo else 0)/hd+1j*d['DASTk'][j]
    if nd.real<=0:raise ValueError('DAST real index becomes nonpositive.')
    s=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    s.add_output_layer(eps=nsub*nsub);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
    def uniform(h,eps):
        if h>0:s.add_layer(thickness=h,eps=complex(eps))
    def patterned(h,eps):
        if h>0:s.add_layer(thickness=h,eps=pm+mask*(complex(eps)-pm))
    if geometry=='old10':
        # Original uploaded 10.py baseline, frozen geometry; same current weights/models.
        uniform(15.,pm);patterned(12.,nv*nv);uniform(23.,pm)
        for h,n in [(3.,d['Ti'][j]),(5.,d['ITO'][j]),(95.,nd),(5.,d['ITO'][j]),(25.,d['SiO2'][j]),(12.,d['Au'][j])]:uniform(h,n*n)
    else:
        uniform(H_PMMA_TOP,pm)
        if geometry in ('pattern','no_vo2'):
            if geometry=='no_vo2':uniform(H_VO2,pm)  # Replace VO2 by host, preserve all interface positions.
            else:patterned(H_VO2,nv*nv)
            patterned(H_SIO2_PATCH,d['SiO2'][j]**2);patterned(H_AU,d['Au'][j]**2)
        elif geometry in ('full_film','equal_volume'):
            f=1. if geometry=='full_film' else (WIDTH/P)**2
            # Center each equivalent film within its original z slab; host fills remainder.
            for h,n in [(H_VO2,nv),(H_SIO2_PATCH,d['SiO2'][j]),(H_AU,d['Au'][j])]:
                uniform((1-f)*h/2,pm);uniform(f*h,n*n);uniform((1-f)*h/2,pm)
        else:raise ValueError('Unknown geometry '+geometry)
        for h,n in [(H_ITO1,d['ITO'][j]),(H_DAST,nd),(H_ITO2,d['ITO'][j])]:uniform(h,n*n)
    s.solve_global_smatrix();return powers(s,lam,nsub,order,device)


def spectrum(d,phase,bias,device,order,mask,geometry='pattern',label='',eo=True):
    t=time.perf_counter();q=[];n=len(d['w'])
    with torch.no_grad():
        for j in range(n):
            q.append(solve_one(d,j,phase,bias,device,order,mask,geometry,eo))
            if n>50 and (j+1)%max(50,n//4)==0:print(f'  {label}: {j+1}/{n}, {time.perf_counter()-t:.1f}s',flush=True)
    if device.type=='cuda':torch.cuda.synchronize()
    print(f'{label or geometry}: {n} wavelengths, {time.perf_counter()-t:.1f}s',flush=True)
    return np.asarray(q)


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



# Four inexpensive configurations; independent one-factor comparisons, NOT a full convergence study.
CANDIDATES=[('A_Au4_Si20',4.,20.),('B_Au0_Si60',0.,60.)]
LEVELS=[('Base',40.,2,32),('Wave',20.,2,32),('Order',40.,3,32),('Pixel',40.,3,48)]
PAIRS=[('Base','Wave','wavelength'),('Base','Order','Fourier'),('Order','Pixel','pixels')]
VEND=8.5; DA_TOL=0.001; TLUM_TOL=0.005; FONT_SIZE=11

def summary_rows(result):
    rows=[]
    for (design,level,phase),g in result.groupby(['Design','Level','Phase'],sort=False):
        a=g.sort_values('Voltage').iloc[0];b=g.sort_values('Voltage').iloc[-1]
        pm=min(a.DA,-b.DA);mp=min(-a.DA,b.DA);margin=max(pm,mp)
        rows.append(dict(Design=design,Level=level,Phase=phase,DAminus=a.DA,DAplus=b.DA,
            Margin=margin,Orientation='+ to -' if pm>=mp else '- to +',Reversal=margin>0,
            MeetsMargin=margin>=MARGIN,Tlum_endpoint_min=min(a.Tlum,b.Tlum),Span=abs(b.DA-a.DA)))
    return pd.DataFrame(rows)

def main():
    global H_AU,H_SIO2_PATCH
    t0=time.perf_counter();torch.set_num_threads(1)
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise ValueError('Use flag=0 CPU or available flag=1 CUDA.')
    if VEND<=0 or VEND*1000/H_DAST>E_LIMIT:raise ValueError('Endpoint voltage violates field limit.')
    device=torch.device('cuda' if flag else 'cpu');torch.backends.cuda.matmul.allow_tf32=False
    print('Air/PMMA15/[VO2(12)/SiO2(s)/Au(a) aligned squares]/ITO5/DAST95/ITO5/SiO2 half-space')
    print(f'P={P:g}, width={WIDTH:g} nm; V=+/-{VEND:g} V; field={VEND*1000/H_DAST:.2f} V/um; DAST L1/E1; CPU threads=1')
    print('Configurations:',LEVELS,'; comparisons:',PAIRS,flush=True)
    ds={step:load_data(step) for _,step,_,_ in LEVELS};cache={};spec={};rows=[]
    checks(ds[LEVELS[0][1]],device,LEVELS[0][2])
    for name,au,sp in CANDIDATES:
        H_AU=au;H_SIO2_PATCH=sp
        for label,step,order,pixels in LEVELS:
            d=ds[step];mask=mask_square(WIDTH,pixels,pixels,device)
            for phase in ('I','M'):
                for v in (-VEND,VEND):
                    q=[];tic=time.perf_counter();nnew=0
                    with torch.no_grad():
                        for j,w in enumerate(d['w']):
                            key=(name,phase,v,order,pixels,round(float(w),8))
                            if key not in cache:cache[key]=solve_one(d,j,phase,v,device,order,mask);nnew+=1
                            q.append(cache[key])
                    q=np.asarray(q);spec[(name,label,phase,v)]=(d['w'],q)
                    row=dict(Design=name,Au_nm=au,SiO2_nm=sp,Level=label,Step_nm=step,N=order,Pixels=pixels,Phase=phase,Voltage=v,**metrics(d,q));rows.append(row)
                    print(f'{name} {label} {phase} {v:+g}V: DA={row["DA"]:+.6f}, T_min={row["Tlum"]:.4f}; new wavelengths={nnew}; {time.perf_counter()-tic:.1f}s',flush=True)
    result=pd.DataFrame(rows);summary=summary_rows(result);diff=[];status=[]
    for name,_,_ in CANDIDATES:
        for phase in ('I','M'):
            sub=summary[(summary.Design==name)&(summary.Phase==phase)].set_index('Level')
            for left,right,kind in PAIRS:
                a,b=sub.loc[left],sub.loc[right]
                x=result[(result.Design==name)&(result.Phase==phase)&(result.Level==left)].sort_values('Voltage')
                y=result[(result.Design==name)&(result.Phase==phase)&(result.Level==right)].sort_values('Voltage')
                da=max(abs(b.DAminus-a.DAminus),abs(b.DAplus-a.DAplus));dt=float(np.max(abs(y[['TlumF','TlumB']].to_numpy()-x[['TlumF','TlumB']].to_numpy())))
                diff.append(dict(Design=name,Phase=phase,Change=kind,From=left,To=right,MaxDeltaDA=da,DeltaMargin=abs(b.Margin-a.Margin),MaxDeltaTlum=dt,SignsSame=bool(np.sign(a.DAminus)==np.sign(b.DAminus) and np.sign(a.DAplus)==np.sign(b.DAplus)),ThresholdSame=bool(a.MeetsMargin==b.MeetsMargin),WithinScreenTolerance=da<=DA_TOL and dt<=TLUM_TOL))
            dd=[r for r in diff if r['Design']==name and r['Phase']==phase]
            status.append(dict(Design=name,Phase=phase,AllLevelsReverse=bool(sub.Reversal.all()),AnyLevelReverses=bool(sub.Reversal.any()),AllLevelsMeetMargin=bool(sub.MeetsMargin.all()),EndpointSignsStable=all(r['SignsSame'] for r in dd),AllPairsWithinTolerance=all(r['WithinScreenTolerance'] for r in dd),MinMarginAcrossLevels=sub.Margin.min(),MaxMarginAcrossLevels=sub.Margin.max()))
    differences=pd.DataFrame(diff);status=pd.DataFrame(status);thermal=[]
    for (name,level,v),g in result.groupby(['Design','Level','Voltage'],sort=False):
        g=g.set_index('Phase');a,b=g.loc['I','DA'],g.loc['M','DA'];m=max(min(a,-b),min(-a,b))
        thermal.append(dict(Design=name,Level=level,Voltage=v,DA_I=a,DA_M=b,PhaseMargin=m,OppositeSigns=m>0,MeetsMargin=m>=MARGIN,Tlum_min=g.Tlum.min()))
    thermal=pd.DataFrame(thermal)
    for title,t in [('ENDPOINT RESULTS',result),('REVERSAL SUMMARY',summary),('ONE-FACTOR DIFFERENCES',differences),('AUDIT STATUS',status),('FIXED-BIAS PHASE COMPARISON',thermal)]:print('\n'+title+'\n'+t.to_string(index=False,float_format=lambda x:f'{x:.6g}'))
    print(f'\nScreen tolerances: DA={DA_TOL:g}, Tlum={TLUM_TOL:g}. Stability refers ONLY to tested configurations, not proven convergence.')
    print('Wave checks N=2 at 20/40 nm; Order checks N=2/3 at 40 nm; Pixel checks 32/48 at N=3,40 nm. Combined N=3/20nm refinement is NOT tested.')
    print('Reverse source is the SiO2 half-space, not rear air. Au=0 removes that slab (same previous control convention). No thermal feedback or thickness-dependent nk.')
    plt.rcParams.update({'font.size':FONT_SIZE,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in'})
    fig,axs=plt.subplots(2,3,figsize=(16,9),constrained_layout=True);axes=axs.ravel();labels=[x[0] for x in LEVELS];colors={'I':'tab:blue','M':'tab:red'}
    for i,(name,au,sp) in enumerate(CANDIDATES):
        for phase in ('I','M'):
            g=summary[(summary.Design==name)&(summary.Phase==phase)].set_index('Level').loc[labels];x=np.arange(len(labels))
            axes[i].plot(x,g.DAminus,'o--',color=colors[phase],label=f'{phase}, -{VEND:g} V');axes[i].plot(x,g.DAplus,'s-',color=colors[phase],label=f'{phase}, +{VEND:g} V')
            axes[3+i].plot(x,g.Margin,'o-',color=colors[phase],label=phase)
            axes[2].plot(x,g.Tlum_endpoint_min,'o-' if i==0 else 's--',color=colors[phase],label=f'{name[0]} {phase}')
        axes[i].axhline(0,color='.5',lw=.8);axes[i].set_title(f'({chr(97+i)}) {name}: endpoint absorption difference',loc='left');axes[i].set_ylabel('Solar-weighted A_F - A_B')
        axes[3+i].axhline(0,color='.5',lw=.8);axes[3+i].axhline(MARGIN,color='.5',ls=':');axes[3+i].set_title(f'({chr(100+i)}) {name}: reversal margin',loc='left');axes[3+i].set_ylabel('Either-orientation margin')
    axes[2].set_title('(c) Minimum endpoint luminous T',loc='left');axes[2].set_ylabel('Min(F,B), D65 weighted')
    for name,_,_ in CANDIDATES:
        for phase in ('I','M'):
            g=differences[(differences.Design==name)&(differences.Phase==phase)];axes[5].plot(range(len(g)),g.MaxDeltaDA,'o-' if name.startswith('A') else 's--',color=colors[phase],label=f'{name[0]} {phase}')
    axes[5].axhline(DA_TOL,color='.5',ls=':');axes[5].set(xticks=range(3),xticklabels=['Wavelength','Fourier','Pixels'],ylabel='Max absolute endpoint DA change');axes[5].set_title('(f) One-factor changes',loc='left')
    for ax in axes[:5]:ax.set(xticks=range(len(labels)),xticklabels=labels)
    for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.18)
    print(f'Elapsed={time.perf_counter()-t0:.1f}s; ONE figure; no files exported.')
    if SHOW_FIGURE:plt.show()
    return dict(ResultTable=result,SummaryTable=summary,DifferenceTable=differences,StatusTable=status,PhaseTable=thermal,Spectra=spec,Figure=fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path);p.add_argument('--no-show',action='store_true');args=p.parse_args()
    if args.project is not None:PROJECT=args.project.resolve()
    if args.no_show:SHOW_FIGURE=0
    Results=main();globals().update(Results)
