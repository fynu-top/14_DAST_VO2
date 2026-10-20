
#!/usr/bin/env python3
"""True 2D aligned VO2/SiO2/Au square stack in PMMA, above continuous ITO/DAST/ITO.
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
QUICK=0                             # 1: functional check only (40 nm, N=1, 48x48)
STEP_NM=5.; ORDER=3; NX=NY=80         # preliminary working grid, NOT certified converged
P=400.; WIDTH=180.                   # common aligned square width; area fill=(WIDTH/P)^2
H_VO2=12.; H_SIO2_PATCH=60.; H_AU=12.; H_PMMA_TOP=15.
H_ITO1=5.; H_DAST=95.; H_ITO2=5.      # all three layers continuous, no bottom Au sheet
VOLTAGES=(-8.5,0.,8.5); PHASES=('I','M'); E_LIMIT=100.; MARGIN=.005
EO_SIGN=1.; DAST_SOURCE='hybrid'; ITO_FILE='4TIO.xlsx'
LOSSLESS_PMMA=True                   # same transparent-host assumption as 10.py; set False to retain k
RUN_CONTROLS=1; RUN_OLD10=1           # no-VO2 + full-coverage & equal-volume films + old10 at 0V
CHECK_FOURIER=1; CHECK_WAVELENGTH=0  # Fourier/pixel probes; optional extra full spectrum at STEP/2
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


def slice_data(d,indices):return {k:np.asarray(v)[indices] for k,v in d.items()}


def audit(d,spec,device,order,nx,ny):
    # Adaptive spectral probes are useful diagnostics, NOT full-spectrum certification.
    w=d['w'];idx=set(np.argmin(abs(w-x)) for x in [405,500,700,900,1200,1550,2000,2400])
    for q in spec.values():
        u=q.mean(axis=2)
        for y in (u[:,0,2],u[:,1,2],abs(u[:,0,2]-u[:,1,2])):idx.add(int(np.argmax(y)))
    ix=np.array(sorted(idx));small=slice_data(d,ix);rows=[]
    n2=order+1;mask1=mask_square(WIDTH,nx,ny,device);mask2=mask_square(WIDTH,2*nx,2*ny,device)
    with torch.no_grad():
        for (phase,v),q in spec.items():
            a=spectrum(small,phase,v,device,n2,mask1,label=f'Fourier probe {phase}/{v:g}');b=spectrum(small,phase,v,device,n2,mask2,label=f'Pixel probe {phase}/{v:g}')
            diff1=float(np.max(abs(a-q[ix])));diff2=float(np.max(abs(b-a)))
            rows.append(dict(Phase=phase,Voltage=v,ProbeCount=len(ix),Nbase=order,Ncheck=n2,MaxOrderRTA=diff1,MaxPixelRTA=diff2))
    table=pd.DataFrame(rows)
    if table[['MaxOrderRTA','MaxPixelRTA']].to_numpy().max()>1e-3:warnings.warn('Spectral probes differ by >1e-3: increase ORDER/NX/NY before final performance claims.')
    return table


def sketch(ax):
    z=0.;col='#dceeea'
    ax.add_patch(Rectangle((0,z),P,H_PMMA_TOP,facecolor=col,ec='k',lw=.4));z+=H_PMMA_TOP
    for name,h,c in [('VO2',H_VO2,'#d9a176'),('SiO2',H_SIO2_PATCH,'#b5cbe8'),('Au',H_AU,'#e9c76d')]:
        ax.add_patch(Rectangle((0,z),P,h,facecolor=col,ec='gray',lw=.3));ax.add_patch(Rectangle(((P-WIDTH)/2,z),WIDTH,h,facecolor=c,ec='k',lw=.5));ax.text(P/2,z+h/2,name,ha='center',va='center',fontsize=8);z+=h
    for name,h,c in [('ITO',H_ITO1,'#89c3d6'),('DAST',H_DAST,'#d6b1e7'),('ITO',H_ITO2,'#89c3d6')]:
        ax.add_patch(Rectangle((0,z),P,h,facecolor=c,ec='k',lw=.5));ax.text(P/2,z+h/2,name,ha='center',va='center',fontsize=8);z+=h
    ax.text(P/2,-10,'Air / F down',ha='center',fontsize=9);ax.text(P/2,z+13,'SiO2 half-space / B up',ha='center',fontsize=9)
    ax.text(P*.13,H_PMMA_TOP+H_VO2+H_SIO2_PATCH/2,'PMMA',ha='center',fontsize=8,rotation=90)
    ax.set(xlim=(0,P),ylim=(z+30,-25),xlabel='x (nm); square also periodic in y',ylabel='z (nm)',title='(a) Aligned patterned trilayer')


def draw(d,spec,table,control):
    fig,axs=plt.subplots(2,4,figsize=(20,10),constrained_layout=True);ax=axs.ravel();sketch(ax[0]);w=d['w'];colors={'I':'tab:blue','M':'tab:red'}
    for phase in PHASES:
        for v in (min(VOLTAGES),max(VOLTAGES)):
            q=spec[(phase,float(v))].mean(axis=2);ax[1].plot(w,q[:,0,2]-q[:,1,2],color=colors[phase],ls='--' if v<0 else '-',label=f'{phase}, {v:+g} V')
        q=spec[(phase,0.)].mean(axis=2)
        ax[2].plot(w,q[:,0,2],color=colors[phase],label=f'{phase}, F');ax[2].plot(w,q[:,1,2],'--',color=colors[phase],label=f'{phase}, B')
        ax[3].plot(w,q[:,0,1],color=colors[phase],label=f'{phase}, total F');ax[3].plot(w,q[:,0,3],'--',color=colors[phase],label=f'{phase}, zero order F')
        t=table[table.Phase==phase].sort_values('Voltage');ax[4].plot(t.Voltage,t.DA,'-o',color=colors[phase],label=phase)
        ax[5].plot(t.Voltage,t.Tlum,'-o',color=colors[phase],label=f'{phase}, total');ax[5].plot(t.Voltage,t.Tlum0,'--',color=colors[phase],label=f'{phase}, zero order')
    ax[1].set(title='(b) Directional absorptance',xlabel='Wavelength (nm)',ylabel='A_F - A_B');ax[1].axhline(0,color='.6',lw=.8)
    ax[2].set(title='(c) Forward/backward absorption, 0 V',xlabel='Wavelength (nm)',ylabel='A',ylim=(0,1))
    ax[3].set(title='(d) Transmission, 0 V',xlabel='Wavelength (nm)',ylabel='T',ylim=(0,1))
    ax[4].set(title='(e) Electrical response (sampled)',xlabel='DAST voltage (V)',ylabel='Solar-weighted A_F - A_B');ax[4].axhline(0,color='.6',lw=.8)
    ax[5].set(title='(f) Luminous transmittance',xlabel='DAST voltage (V)',ylabel='Min(F,B), D65 weighted',ylim=(0,1))
    i=table[table.Phase=='I'].set_index('Voltage');m=table[table.Phase=='M'].set_index('Voltage');vs=np.sort(np.array(VOLTAGES))
    ax[6].plot(vs,[m.loc[v,'DA']-i.loc[v,'DA'] for v in vs],'-o',label='M - I');ax[6].axhline(0,color='.6',lw=.8)
    ax[6].set(title='(g) Phase contrast at fixed voltage',xlabel='Voltage (V)',ylabel='DA_sol,M - DA_sol,I')
    for phase in PHASES:
        for _,r in control[control.Phase==phase].iterrows():
            ax[7].scatter(r.Tlum,abs(r.DA),color=colors[phase],marker='o' if r.Geometry=='pattern' else 'x');ax[7].annotate(phase+': '+r.Geometry,(r.Tlum,abs(r.DA)),fontsize=6,xytext=(3,3),textcoords='offset points')
    ax[7].set(title='(h) Zero-bias controls',xlabel='Min luminous total T',ylabel='|DA_sol| (not reversal)')
    for a in ax[1:7]:a.legend(fontsize=8);a.grid(alpha=.15)
    ax[7].grid(alpha=.15);fig.suptitle('VO2 / SiO2 / Au patterned trilayer + continuous ITO / DAST / ITO — exploratory',fontsize=13)
    return fig


def main():
    start=time.perf_counter()
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise ValueError('Set flag=0 for CPU; flag=1 requires CUDA.')
    if not 0<WIDTH<P or H_DAST<=0 or min(H_VO2,H_SIO2_PATCH,H_AU,H_ITO1,H_ITO2,H_PMMA_TOP)<0:raise ValueError('Invalid dimensions.')
    if STEP_NM<=0 or ORDER<1 or min(NX,NY)<=4*(ORDER+1):raise ValueError('Invalid wavelength/Fourier/pixel settings.')
    if tuple(PHASES)!=('I','M') or 0. not in VOLTAGES or len(set(VOLTAGES))!=len(VOLTAGES):raise ValueError('Use PHASES=(I,M), unique voltages including 0.')
    if max(map(abs,VOLTAGES))*1000/H_DAST>E_LIMIT+1e-6:raise ValueError('DAST field exceeds E_LIMIT.')
    device=torch.device('cuda' if flag else 'cpu');torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    step,order,nx,ny=(40.,1,48,48) if QUICK else (STEP_NM,ORDER,NX,NY)
    d=load_data(step);mask=mask_square(WIDTH,nx,ny,device)
    print(f'DEVICE={device}; torch={torch.__version__}; complex128; DAST field limit={E_LIMIT:g} V/um')
    print(f'NEW: air / PMMA({H_PMMA_TOP:g}) / aligned [VO2({H_VO2:g}), SiO2({H_SIO2_PATCH:g}), Au({H_AU:g})] squares / ITO({H_ITO1:g}) / DAST({H_DAST:g}) / ITO({H_ITO2:g}) / SiO2 half-space')
    print(f'P={P:g}, width={WIDTH:g} nm, area fill={(WIDTH/P)**2:.5f}; N={order}, Nx/Ny={nx}/{ny}, wavelengths={len(d["w"])}; PMMA lossless={LOSSLESS_PMMA}')
    print(f'ITO uses {ITO_FILE}, identity follows your data assignment; verify source. DAST={DAST_SOURCE}/E1 scalar. No thermal feedback.')
    validation=checks(d,device,order);spec={};rows=[];polrows=[];controls=[];cache={}
    for phase in PHASES:
        for v in sorted(VOLTAGES):
            q=spectrum(d,phase,float(v),device,order,mask,label=f'Pattern {phase} {v:+g}V');spec[(phase,float(v))]=q
            row=dict(Phase=phase,Voltage=float(v),**metrics(d,q));rows.append(row)
            for j,pol in enumerate(('TM','TE')):polrows.append(dict(Phase=phase,Voltage=float(v),Polarization=pol,**metrics(d,q[:,:,j:j+1,:])))
            print(pd.DataFrame([row]).to_string(index=False,float_format=lambda x:f'{x:.6f}'),flush=True)
    table=pd.DataFrame(rows);summ=[]
    for phase in PHASES:
        t=table[table.Phase==phase].sort_values('Voltage');da=t.DA.to_numpy();a,b=da[0],da[-1]
        summ.append(dict(Phase=phase,DAmin_sample=da.min(),DAmax_sample=da.max(),VoltageRangeDA=np.ptp(da),EndpointChange=b-a,EndpointReversal=bool(a*b<0),SampleSignChange=bool(da.min()<0<da.max()),Mrev_plus_to_minus=min(a,-b),Mrev_minus_to_plus=min(-a,b),MinTlum_sample=t.Tlum.min(),MinTlum0_sample=t.Tlum0.min()))
        controls.append(dict(Phase=phase,Geometry='pattern',**metrics(d,spec[(phase,0.)])))
        for v in VOLTAGES:cache[(phase,v)]=dict(Phase=phase,Voltage=v,Geometry='EO_off',**metrics(d,spec[(phase,0.)]))
    if RUN_CONTROLS:
        q=spectrum(d,'I',0.,device,order,mask,'no_vo2',label='No VO2, 0V')
        for phase in PHASES:controls.append(dict(Phase=phase,Geometry='no_vo2',**metrics(d,q)))
        for geometry in ('full_film','equal_volume'):
            for phase in PHASES:
                q=spectrum(d,phase,0.,device,order,mask,geometry,label=f'{geometry} {phase} 0V');controls.append(dict(Phase=phase,Geometry=geometry,**metrics(d,q)))
    if RUN_OLD10:
        oldmask=mask_square(180.,nx,ny,device)
        if P!=400.:warnings.warn('Old10 comparison disabled because P changed from 400 nm.')
        else:
            for phase in PHASES:
                q=spectrum(d,phase,0.,device,order,oldmask,'old10',label=f'Original10 {phase} 0V');controls.append(dict(Phase=phase,Geometry='old10',**metrics(d,q)))
    control=pd.DataFrame(controls);summary=pd.DataFrame(summ);zero=table[table.Voltage==0].set_index('Phase')
    audit_table=audit(d,spec,device,order,nx,ny) if CHECK_FOURIER else pd.DataFrame()
    wave_rows=[]
    if CHECK_WAVELENGTH:
        # One potentially most interesting state: largest |DA|; no claim for unchecked states.
        chosen=table.iloc[np.argmax(abs(table.DA.to_numpy()))];p=chosen.Phase;v=float(chosen.Voltage);fine=load_data(step/2)
        q=spectrum(fine,p,v,device,order,mask,label=f'Wavelength check {p}/{v:g}');a=metrics(d,spec[(p,v)]);b=metrics(fine,q)
        wave_rows=[dict(Phase=p,Voltage=v,Metric=k,Coarse=a[k],Fine=b[k],AbsDiff=abs(a[k]-b[k])) for k in a]
    for name,tab in [('RESULTS',table),('STATE SUMMARY',summary),('ZERO-BIAS CONTROLS',control),('FOURIER/PIXEL PROBES',audit_table),('WAVELENGTH CHECK',pd.DataFrame(wave_rows))]:
        if not tab.empty:print('\n'+name+'\n'+tab.to_string(index=False,float_format=lambda x:f'{x:.6g}'))
    phase_table=table.pivot(index='Voltage',columns='Phase',values='DA');phase_table['M_minus_I']=phase_table['M']-phase_table['I'];print('\nPHASE CONTRAST\n',phase_table)
    rec=max(np.max(abs(q[:,0,:,3]-q[:,1,:,3])) for q in spec.values());pol=max(np.max(abs(q[:,:,0,:]-q[:,:,1,:])) for q in spec.values())
    print(f'Max matched zero-order T reciprocity={rec:.3e}; C4 TE/TM power difference={pol:.3e}')
    if pol>1e-5:warnings.warn('Unexpected TE/TM difference for scalar C4 geometry; check discretization.')
    print('EO_off controls reuse 0V spectra exactly (same n when beta=0). No-VO2 preserves patch-zone thickness by replacing VO2 with PMMA.')
    print('Equal-volume films occupy fraction f of each slab height, centered in PMMA; this is a specified reference, not a unique homogenization.')
    print('Tlum0 is ZERO DIFFRACTION ORDER, not zero voltage. Backward illumination originates in the SiO2 half-space; no rear air/glass interface is included.')
    print('Probe convergence is not a whole-spectrum certificate. 3 voltage points cannot establish continuous range or minimum switching voltage.')
    print(f'Elapsed={time.perf_counter()-start:.1f}s; one figure; no file exports.')
    fig=draw(d,spec,table,control)
    if SHOW_FIGURE:plt.show()
    return dict(ResultTable=table,SummaryTable=summary,ControlTable=control,EOOffTable=pd.DataFrame(cache.values()),PolarizationTable=pd.DataFrame(polrows),AuditTable=audit_table,WavelengthTable=pd.DataFrame(wave_rows),PhaseTable=phase_table,Spectra=spec,Figure=fig,Checks=validation)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--project',type=Path);parser.add_argument('--quick',action='store_true');parser.add_argument('--no-show',action='store_true')
    args=parser.parse_args()
    if args.project is not None:PROJECT=args.project.resolve()
    if args.quick:QUICK=1
    if args.no_show:SHOW_FIGURE=0
    Results=main();globals().update(Results)
