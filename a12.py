
#!/usr/bin/env python3
"""Wu, Meng & Liu, Materials 2020, 13, 3751, doi:10.3390/ma13173751, Fig.2(a,b).
Normal incidence; x-polarized TM; F=air/Ti side, B=SiO2/Au side.
2D RCWA alternative to the paper's FDTD. No fitting, no exports, one figure.
Requires numpy pandas openpyxl matplotlib torch torcwa and optical-constant tables.
Paper cites Palik 1998 for Ti/Au/Al2O3 but does not supply numeric tables.
Use matching Palik tables for quantitative reproduction. No fabricated metal fallback.
"""
from pathlib import Path
import argparse,time,re
import numpy as np,pandas as pd,torch,torcwa
import matplotlib.pyplot as plt
flag=0; ORDER=5; NX=NY=128; STEP_NM=5.; FONT_SIZE=12; SHOW_FIGURE=True
P=400.; WIDTH=180.; H_TI=30.; H_AL=60.; H_AU=40.; H_PMMA=150.
PMMA_MODE='total' # total: PMMA extends 150 nm above substrate (20 nm top cap); cap: extra 150 nm
PROJECT=Path(__file__).resolve().parent; DTYPE=torch.complex128
# Set an explicit filename/path to override auto-discovery; files in optical_constants or project root.
FILES={'Ti':None,'Au':None,'Al2O3':None,'SiO2':None,'PMMA':None}
UNITS={k:'auto' for k in FILES} # auto, nm, um; columns: wavelength, n, k (k optional for dielectrics)
N_SIO2=1.45; N_PMMA=1.49 # explicit approximations ONLY if their tables are absent; paper does not give values
POLARIZATION='p' # normal-incidence p corresponds to x polarization; 'unpolarized' averages p and s
AUDIT=False; AUDIT_ORDERS=[3,5,7]; AUDIT_WAVELENGTHS=[500.,750.,890.,913.5,1100.]
TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz

def locate(name):
    if FILES[name]:
        p=Path(FILES[name]);p=p if p.is_absolute() else PROJECT/p
        if not p.is_file():raise FileNotFoundError(str(p))
        return p
    found=[]
    for folder in (PROJECT/'optical_constants',PROJECT):
        if not folder.is_dir():continue
        for p in folder.iterdir():
            token=re.sub(r'^[0-9_ -]+','',p.stem).lower()
            if p.suffix.lower() in ('.csv','.xlsx','.txt') and token in (name.lower(),name.lower()+'_palik',name.lower()+'_nk'):found.append(p)
    if len(found)>1:raise ValueError(f'Multiple {name} tables: {found}; select FILES[{name!r}] explicitly.')
    return found[0] if found else None

def read_nk(path,w,unit,name):
    f=pd.read_excel(path,header=None) if path.suffix.lower()=='.xlsx' else pd.read_csv(path,header=None,comment='#',sep=None,engine='python')
    a=f.iloc[:,:3].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    if a.shape[1]<2:raise ValueError(f'{path}: need wavelength,n,k columns')
    if a.shape[1]==2:
        if name in ('Ti','Au'):raise ValueError(f'{path}: metal k column missing')
        a=np.column_stack([a,np.zeros(len(a))])
    a=a[np.isfinite(a).all(axis=1)];a=a[np.argsort(a[:,0])];a=a[np.unique(a[:,0],return_index=True)[1]]
    if len(a)<3:raise ValueError(f'{path}: fewer than 3 usable rows')
    if unit=='auto':unit='um' if np.median(a[:,0])<50 else 'nm'
    if unit not in ('um','nm'):raise ValueError('UNITS must be auto/um/nm')
    x=a[:,0]*(1000 if unit=='um' else 1)
    if x[0]>w.min() or x[-1]<w.max():raise ValueError(f'{path}: {x[0]:g}-{x[-1]:g} nm does not cover {w.min():g}-{w.max():g} nm; no extrapolation')
    if np.any(a[:,1]<=0) or np.any(a[:,2]<0):raise ValueError(f'{path}: require passive n>0,k>=0')
    n=np.interp(w,x,a[:,1]);k=np.interp(w,x,a[:,2]);print(f'{name}: {path}; wavelength unit={unit}')
    if name in ('PMMA','SiO2'):
        if k.max()>1e-6:print(f'NOTICE {name}: ignoring k for lossless coating/exterior model')
        k=np.zeros_like(k)
    return n+1j*k

def load_data():
    w=np.unique(np.r_[np.arange(400.,1300.00001,STEP_NM),400.,750.,760.,890.,913.5,1300.,AUDIT_WAVELENGTHS if AUDIT else []])
    paths={k:locate(k) for k in FILES};missing=[k for k in ('Ti','Au','Al2O3') if paths[k] is None]
    if missing:raise FileNotFoundError('Missing material tables: '+', '.join(missing)+'. Put Ti.csv, Au.csv, Al2O3.csv (wavelength,n,k) in optical_constants, or set FILES. Paper source: Palik 1998. No substitute metal model is used.')
    d={'w':w}
    for k,path in paths.items():
        if path:d[k]=read_nk(path,w,UNITS[k],k)
        else:
            n=N_SIO2 if k=='SiO2' else N_PMMA;d[k]=np.full(w.size,n+0j);print(f'APPROXIMATION: {k} constant n={n}; supply a table to replace it.')
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



def solve_one(d,j,order,mask,device):
    w=float(d['w'][j]);ns=float(d['SiO2'][j].real);lam=safe_lambda(w,ns,order);pm=complex(d['PMMA'][j]**2)
    cap=H_PMMA-H_TI-H_AL-H_AU if PMMA_MODE=='total' else H_PMMA
    if cap<0:raise ValueError('Total PMMA thickness is smaller than the MIM stack.')
    sim=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    sim.add_output_layer(eps=ns*ns);sim.set_incident_angle(inc_ang=0.,azi_ang=0.)
    if cap>0:sim.add_layer(thickness=cap,eps=pm)
    for h,key in [(H_TI,'Ti'),(H_AL,'Al2O3'),(H_AU,'Au')]:sim.add_layer(thickness=h,eps=pm+mask*(complex(d[key][j]**2)-pm))
    sim.solve_global_smatrix();return powers(sim,lam,ns,order,device)

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





def main():
    started=time.perf_counter();torch.set_num_threads(1)
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise ValueError('flag=1 requires CUDA; use flag=0 for CPU')
    if PMMA_MODE not in ('total','cap') or POLARIZATION not in ('p','unpolarized'):raise ValueError('Invalid PMMA_MODE or POLARIZATION')
    if ORDER<1 or STEP_NM<=0 or not 0<WIDTH<P or min(NX,NY)<=4*max([ORDER]+(AUDIT_ORDERS if AUDIT else [])):raise ValueError('Invalid geometry or numerical settings')
    device=torch.device('cuda' if flag else 'cpu');torch.backends.cuda.matmul.allow_tf32=False;d=load_data();checks(d,device,ORDER)
    print(f'Wu Fig.2: P={P:g}, w={WIDTH:g}; Ti/Al2O3/Au={H_TI:g}/{H_AL:g}/{H_AU:g} nm; PMMA mode={PMMA_MODE}; N={ORDER}; pixels={NX}x{NY}.')
    mask=mask_square(WIDTH,NX,NY,device);spectra=[]
    with torch.no_grad():
        for j,w in enumerate(d['w']):
            spectra.append(solve_one(d,j,ORDER,mask,device))
            if j%10==0 or j==len(d['w'])-1:print(f'{j+1}/{len(d["w"])} wavelengths; {w:g} nm; elapsed {time.perf_counter()-started:.1f}s',flush=True)
    q=np.asarray(spectra);u=q[:,:,0,:] if POLARIZATION=='p' else q.mean(axis=2);w=d['w'];rows=[]
    for di,label in enumerate(('Front (air/Ti)','Back (SiO2/Au)')):
        nir=w>=750.;ids=np.flatnonzero(nir);ia=ids[np.argmax(u[nir,di,2])];ir=ids[np.argmax(u[nir,di,0])]
        vis=(w>=400)&(w<=750);tv=TRAPZ(u[vis,di,1],w[vis])/(750-400)
        rows.append(dict(Direction=label,PeakA_NIR=u[ia,di,2],LambdaA_nm=w[ia],PeakR_NIR=u[ir,di,0],LambdaR_nm=w[ir],MeanT_400_750=tv))
    summary=pd.DataFrame(rows);print('\n'+summary.to_string(index=False));print(f'Max matched-zero-order |TF0-TB0|={np.max(abs(u[:,0,3]-u[:,1,3])):.3e}')
    print('Visible mean is unweighted wavelength mean, NOT D65 luminous transmission. Total T may differ in the diffractive visible range; zero-order T obeys reciprocity.')
    audit=[]
    if AUDIT:
        with torch.no_grad():
            for n in AUDIT_ORDERS:
                for lam in AUDIT_WAVELENGTHS:
                    j=int(np.argmin(abs(w-lam)));z=q[j] if n==ORDER else solve_one(d,j,n,mask,device);v=z[:,0,:] if POLARIZATION=='p' else z.mean(axis=1)
                    audit.append(dict(N=n,Wavelength=w[j],RF=v[0,0],TF=v[0,1],AF=v[0,2],RB=v[1,0],TB=v[1,1],AB=v[1,2],MaxDiffFromMain=np.max(abs(v[:,:3]-u[j,:,:3]))))
        print('\nORDER SPOT CHECK (not full convergence)\n'+pd.DataFrame(audit).to_string(index=False))
    plt.rcParams.update({'font.family':'Times New Roman','font.size':FONT_SIZE,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in'})
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    for di,a in enumerate(ax):
        for col,color,label in [(0,'tab:blue','R'),(1,'tab:green','T'),(2,'red','A')]:a.plot(w,u[:,di,col],color=color,lw=1.6,label=label)
        a.set(xlim=(400,1300),ylim=(0,1),xlabel='Wavelength (nm)',ylabel='R / T / A');a.set_title(('(a) Front illumination','(b) Back illumination')[di],loc='left',fontsize=FONT_SIZE+2);a.legend();a.grid(alpha=.15)
    print(f'Elapsed={time.perf_counter()-started:.1f}s. One figure, no file exports. RCWA agreement with published FDTD is not guaranteed before material matching and convergence.')
    if SHOW_FIGURE:plt.show()
    return dict(Wavelength=w,RTA=u[:,:,:3],Tzero=u[:,:,3],PolarizedSpectra=q,Summary=summary,Audit=pd.DataFrame(audit),Figure=fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path);p.add_argument('--no-show',action='store_true');p.add_argument('--audit',action='store_true');a=p.parse_args()
    if a.project:PROJECT=a.project.resolve()
    if a.no_show:SHOW_FIGURE=False
    if a.audit:AUDIT=True
    Results=main();globals().update(Results)
