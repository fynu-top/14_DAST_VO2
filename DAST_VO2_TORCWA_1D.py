#!/usr/bin/env python3
"""Pure-Python TORCWA simulation: Air/ITO/DAST/ITO/SiO2/VO2-strips/Air.

Read ORIGINAL measured/tabulated Excel/CSV data directly from the local project;
MATLAB, Janus.m, and Janus.py are NOT called. No synthetic optical constants.
    python -m pip install numpy scipy pandas openpyxl matplotlib torch torcwa
    python DAST_VO2_TORCWA_1D.py --project "D:/09_PCloud/14_DAST_VO2" --quick
    python DAST_VO2_TORCWA_1D.py --project "D:/09_PCloud/14_DAST_VO2"

Optional --torcwa-path points at a checked-out torcwa repository root. Source
complex indices from spreadsheets are n - i*k (MATLAB notation), converted to
n + i*k as required by TORCWA exp(-i*omega*t).

Trial geometry, not optimized; normal incidence, 1D x-periodic VO2 strips,
TE=Ey / TM=Ex, scalar DAST E1/L1 model, fixed insulating/metallic VO2 states.
A single 2x3 figure is displayed; no figures or data are exported. The --quick
mode is for functional testing only; NOT adequate for quantitative reporting.
"""
from pathlib import Path
import argparse, importlib, sys, time, warnings
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
import torch
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

# Editable starting geometry; all lengths in nm. VO2 occupies f of each period.
D_ITO1, D_DAST, D_ITO2, D_SIO2 = 30., 85., 30., 175.
PERIOD, FILL, H_VO2 = 300., 0.50, 40.
EMAX, VHARD, MARGIN, EO_SIGN = 100., 10., 0.005, 1.
DTYPE = torch.complex128
METRICS = ['DA_sol','Tlum_F','Tlum_B','Tlum0order_F','Tlum0order_B',
           'Tsol_F','Tsol_B','Asol_F','Asol_B']
TRAPZ = np.trapezoid if hasattr(np,'trapezoid') else np.trapz


def arguments():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--project',type=Path,default=Path(__file__).resolve().parent,help='Folder with optical_constants/ and VO2_temperature_data/')
    p.add_argument('--torcwa-path',type=Path,help='Downloaded repository root containing torcwa/')
    p.add_argument('--device',choices=['auto','cpu','cuda'],default='auto')
    p.add_argument('--step',type=float,default=10.,help='Initial wavelength spacing, nm')
    p.add_argument('--order',type=int,default=7,help='x Fourier truncation N: 2N+1 orders; y order=0')
    p.add_argument('--check-order',type=int,default=11)
    p.add_argument('--nx',type=int,default=512,help='Real-space x grid; y invariant')
    p.add_argument('--nv',type=int,default=9,help='Odd number of voltage samples including zero')
    p.add_argument('--no-check',action='store_true',help='Skip endpoint convergence audit only')
    p.add_argument('--quick',action='store_true',help='Fast functionality test: coarse grid; disables convergence audit')
    p.add_argument('--no-show',action='store_true',help='Create one figure without blocking display')
    a=p.parse_args()
    if a.quick: a.step,a.order,a.check_order,a.nx,a.nv,a.no_check=50.,3,5,96,3,True
    if a.step<=0 or a.order<1 or a.check_order<=a.order or a.nx<=4*a.check_order or a.nv<3 or a.nv%2!=1:
        p.error('Require step>0, order>=1, check-order>order, nx>4*check-order and odd nv>=3.')
    return a


def get_torcwa(path):
    if path is not None:
        path=path.resolve()
        if not (path/'torcwa'/'__init__.py').is_file(): raise FileNotFoundError(f'Not a torcwa repository root: {path}')
        sys.path.insert(0,str(path))
    try: return importlib.import_module('torcwa')
    except ImportError as e:
        raise RuntimeError('Install torcwa or pass --torcwa-path "D:/.../torcwa-main".') from e


def find_project(folder):
    """Locate project using file presence, not the current working directory."""
    root=folder.expanduser().resolve()
    choices=[root,root/'02_MATLAB_All_Fig',root/'02_Python',*list(root.parents)[:3]]
    for candidate in choices:
        if (candidate/'optical_constants'/'4TIO.xlsx').is_file() and (candidate/'VO2_temperature_data'/'ASTMG173.xlsx').is_file():
            return candidate
    raise FileNotFoundError(f'Cannot find project material folders under {root}. Use --project with your 14_DAST_VO2 directory.')


def table_xy(file,cols,min_rows=3):
    """Unique, sorted wavelength + data columns; original wavelength unit retained."""
    if not file.is_file(): raise FileNotFoundError(f'Missing real source data: {file}')
    if file.suffix.lower()=='.csv': df=pd.read_csv(file,header=None,comment='#',on_bad_lines='skip')
    else: df=pd.read_excel(file,header=None,engine='openpyxl')
    if df.shape[1]<=max(cols): raise ValueError(f'{file.name}: expected at least {max(cols)+1} columns')
    a=df.iloc[:,list(cols)].apply(pd.to_numeric,errors='coerce').to_numpy(dtype=float)
    a=a[np.all(np.isfinite(a),axis=1)]; a=a[np.argsort(a[:,0])]
    if len(a)<min_rows: raise ValueError(f'{file.name}: insufficient numeric rows ({len(a)})')
    a=a[np.unique(a[:,0],return_index=True)[1]]
    if len(a)<min_rows: raise ValueError(f'{file.name}: insufficient distinct wavelengths')
    return a


def interp_real(file,w_nm,cols,unit='nm',outside='error',min_rows=3):
    a=table_xy(file,cols,min_rows); x=a[:,0]*(1000 if unit=='um' else 1)
    if outside=='error' and (x[0]>w_nm.min()+1e-6 or x[-1]<w_nm.max()-1e-6):
        raise ValueError(f'{file.name}: measured range {x[0]:g}–{x[-1]:g} nm does not cover {w_nm.min():g}–{w_nm.max():g} nm; no invented extrapolation.')
    return np.column_stack([np.interp(w_nm,x,a[:,j],left=0 if outside=='zero' else None,right=0 if outside=='zero' else None) for j in range(1,a.shape[1])])


def load_user_data(a):
    """Match Janus.m init(0.15625,true,false) without any MATLAB execution."""
    root=find_project(a.project); w=wavelengths(a.step)
    mat=root/'optical_constants'; env=root/'VO2_temperature_data'; atmo=root/'Atmospheric data'
    # nk_data.m originally loads 4TIO.xlsx in material index 4 (called ITO).
    # Janus.init(...,physicalSiO2=true) explicitly overwrites index 2 with 2SiO2.xlsx.
    optical={k:mat/f for k,f in [('nITO','4TIO.xlsx'),('nSiO2','2SiO2.xlsx'),('nVI','6VO2_I.xlsx'),('nVM','7VO2_M.xlsx')]}
    d={'w':w}
    for k,file in optical.items():
        nk=interp_real(file,w,(0,1,2),unit='um')
        if np.any(nk[:,0]<=0) or np.any(nk[:,1]<-1e-9): raise ValueError(f'{file.name}: nonpassive or nonpositive n,k values')
        d[k]=nk[:,0]-1j*np.maximum(nk[:,1],0) # original MATLAB n - i*k
    solar=interp_real(env/'ASTMG173.xlsx',w,(0,2),outside='zero')[:,0]
    photofile=atmo/'CIE_Vlambda.xlsx'
    if not photofile.is_file(): photofile=env/'CIE_Vlambda.csv'
    if not photofile.is_file(): raise FileNotFoundError('Missing CIE_Vlambda.xlsx or CIE_Vlambda.csv')
    photo_df=pd.read_excel(photofile,header=None,engine='openpyxl') if photofile.suffix.lower()=='.xlsx' else pd.read_csv(photofile,header=None)
    col=2 if photo_df.shape[1]>=4 else 1 if photo_df.shape[1]==2 else None
    if col is None: raise ValueError(f'{photofile.name}: expected [w,V] or [w,xbar,ybar,zbar]')
    photopic=interp_real(photofile,w,(0,col),outside='zero',min_rows=100)[:,0]
    d65file=atmo/'CIE_std_illum_D65.csv'
    if not d65file.is_file(): d65file=atmo/'CIE_std_illum_D65.xlsx'
    d65=interp_real(d65file,w,(0,1),outside='zero',min_rows=100)[:,0]
    d.update(solar=solar,photopic=photopic,d65=d65)
    print(f'Read original material/solar/CIE/D65 tables directly with Python: {root}',flush=True)
    print('Material index 4 in nk_data.m uses 4TIO.xlsx; verify it represents your intended ITO experimentally.',flush=True)
    return validate_data(d)


def validate_data(raw):
    keys=['w','nITO','nSiO2','nVI','nVM','solar','photopic','d65']
    missing=[k for k in keys if k not in raw]
    if missing: raise ValueError('Missing material vectors: '+', '.join(missing))
    d={k:np.asarray(raw[k]).ravel().copy() for k in keys}; n=len(d['w'])
    if n<3 or any(len(v)!=n or not np.all(np.isfinite(v)) for v in d.values()): raise ValueError('Material vectors have invalid lengths or NaN/Inf.')
    if np.any(np.diff(d['w'])<=0) or d['w'][0]>400 or d['w'][-1]<2500: raise ValueError('w must increase and cover 400–2500 nm.')
    for key in ['nITO','nSiO2','nVI','nVM']:
        if np.min(d[key].real)<=0 or np.max(d[key].imag)>1e-10: raise ValueError(f'{key}: expected passive MATLAB n-ik, not n+ik.')
        d[key]=np.conj(d[key])  # TORCWA exp(-i omega t): passive index is n+ik.
    for key in ['solar','photopic','d65']:
        if np.iscomplexobj(d[key]) and np.any(np.abs(d[key].imag)>1e-12): raise ValueError(f'{key} must be real.')
        d[key]=d[key].real
        if np.any(d[key]<0): raise ValueError(f'{key} contains negative spectral weights.')
    vis=(d['w']>=400)&(d['w']<=780)
    if TRAPZ(d['solar'],d['w'])<=0 or TRAPZ(d['photopic'][vis]*d['d65'][vis],d['w'][vis])<=0:
        raise ValueError('Zero solar / D65*photopic normalization.')
    return d


def wavelengths(step):
    # Include material transitions and photopic integration endpoints exactly.
    return np.unique(np.r_[np.arange(400.,2500.,step),400.,430.,457.9,488.,514.5,650.,700.,720.,750.,780.,900.,1550.,2100.,2500.])


def optical_data(raw,w):
    d={k:np.interp(w,raw['w'],v) for k,v in raw.items() if k!='w'}; d['w']=w
    n=np.full_like(w,2.2); mask=w>650; x=w[mask]
    ns=np.sqrt(2.078**2+1.645*533**2/(x*x-533**2)); s=.5*(1-np.cos(np.pi*np.clip((x-650)/100,0,1)))
    n[mask]=(1-s)*2.2+s*ns
    alpha=np.full_like(w,1.9e5); mask=(w>430)&(w<=514.5)
    alpha[mask]=10**PchipInterpolator([430,457.9,488,514.5],np.log10([1.9e5,4.1e5,4.9e5,4.79e5]))(w[mask])
    mask=(w>514.5)&(w<=900); alpha[mask]=10**np.interp(w[mask],[514.5,700,750,900],np.log10([4.79e5,5,1.8,.8]))
    mask=(w>900)&(w<=2100); x=w[mask]
    alpha[mask]=.7+np.exp(-.5*((x-1100)/55)**2)+1.5*np.exp(-.5*((x-1400)/70)**2)+3*np.exp(-.5*((x-1700)/95)**2)
    alpha[w>2100]=.7; r=PchipInterpolator([720,750,1550],[530,445,200])(np.clip(w,720,1550))
    d.update(n0=n,kDAST=alpha*w*1e-7/(4*np.pi),beta=EO_SIGN*.5*n**3*r*1e-3)
    d['lum']=d['photopic']*d['d65']; return d


class Solver:
    def __init__(self,tc,device,period=PERIOD,fill=FILL,height=H_VO2):
        self.tc,self.device,self.period,self.fill,self.height=tc,device,period,fill,height
        self.masks={}; self.max_t0_error=0.; self.max_cross_power=0.; self.cutoff_shifts=0
    def mask(self,nx):
        if nx not in self.masks:
            # Cell-average filling: preserves f exactly, converges to sharp strips.
            l=np.arange(nx)/nx; r=l+1/nx; a=(1-self.fill)/2; b=(1+self.fill)/2
            rho=np.maximum(0,np.minimum(r,b)-np.maximum(l,a))*nx
            self.masks[nx]=torch.as_tensor(rho[:,None],dtype=DTYPE,device=self.device)
        return self.masks[nx]
    def one(self,lam,nlist,order,nx,uniform=None):
        solve_lam=float(lam)
        # Avoid singular exactly grazing external orders; report all such offsets.
        m=np.arange(1,order+1)
        if np.any(np.abs(m*lam/self.period-1)<1e-10):
            solve_lam=lam*(1+1e-7); self.cutoff_shifts+=1
        s=self.tc.rcwa(freq=1/solve_lam,order=[order,0],L=[self.period,self.period],dtype=DTYPE,device=self.device,stable_eig_grad=False)
        s.set_incident_angle(inc_ang=0.,azi_ang=0.)  # Air on both sides by default.
        for h,n in zip([D_ITO1,D_DAST,D_ITO2,D_SIO2],nlist[:4]): s.add_layer(thickness=h,eps=complex(n*n))
        h=self.height if uniform is None else uniform
        if h>0:
            eps=complex(nlist[4]**2) if uniform is not None else 1+self.mask(nx)*(complex(nlist[4]**2)-1)
            s.add_layer(thickness=h,eps=eps)
        s.solve_global_smatrix()
        orders=np.array([[j,0] for j in range(-order,order+1) if abs(j*solve_lam/self.period)<1-1e-12],dtype=int)
        if order*solve_lam/self.period<1: raise ValueError('Fourier order too small to include all propagating orders.')
        zero=np.flatnonzero(orders[:,0]==0)[0]; result=np.zeros((2,2,4))
        # Output axes: direction(F,B), polarization(TE,TM), quantities(R,T,A,T0).
        for di,direction in enumerate(['forward','backward']):
            for pi,inc in enumerate(['y','x']):
                for qi,port in enumerate(['reflection','transmission']):
                    power=np.zeros(len(orders))
                    for out in ['x','y']:
                        amp=s.S_parameters(orders=orders,direction=direction,port=port,polarization=out+inc,ref_order=[0,0],power_norm=True)
                        v=amp.abs().square().detach().cpu().numpy().ravel()
                        if not np.all(np.isfinite(v)): raise FloatingPointError('Nonfinite S parameter.')
                        if out!=inc: self.max_cross_power=max(self.max_cross_power,float(v.sum()))
                        power+=v
                    result[di,pi,qi]=power.sum()
                    if qi==1: result[di,pi,3]=power[zero]
                result[di,pi,2]=1-result[di,pi,0]-result[di,pi,1]
        if np.min(result)<-2e-6 or np.max(result)>1+2e-6: raise FloatingPointError(f'Passivity failed at {lam:g} nm: check signs/order/materials.')
        # Reciprocity of matched zero-order channels, NOT total multimode T.
        self.max_t0_error=max(self.max_t0_error,float(np.max(np.abs(result[0,:,3]-result[1,:,3]))))
        if self.max_t0_error>2e-6 or self.max_cross_power>2e-6: raise FloatingPointError('Zero-order reciprocity / 1D polarization-decoupling check failed.')
        return result
    def spectrum(self,d,phase,voltage,order,nx):
        nD=d['n0']+d['beta']*voltage/D_DAST+1j*d['kDAST']
        if np.any(nD.real<=0): raise ValueError('DAST real index nonpositive under this voltage.')
        nv=d['nVI'] if phase=='I' else d['nVM']; output=[]
        for j,lam in enumerate(d['w']):
            ns=[d['nITO'][j],nD[j],d['nITO'][j],d['nSiO2'][j],nv[j]]
            output.append(self.one(lam,ns,order,nx))
        return np.asarray(output)


def weighted(d,spec):
    u=spec.mean(axis=2); w=d['w']; vis=(w>=400)&(w<=780)
    def avg(x,lum=False):
        ii=vis if lum else np.ones(len(w),dtype=bool); weight=d['lum'] if lum else d['solar']
        norm=TRAPZ(weight[ii],w[ii])
        if norm<=0: raise ValueError('Zero photopic or solar normalization.')
        return TRAPZ(x[ii]*weight[ii],w[ii])/norm
    af,ab=avg(u[:,0,2]),avg(u[:,1,2])
    return np.array([af-ab,avg(u[:,0,1],True),avg(u[:,1,1],True),avg(u[:,0,3],True),avg(u[:,1,3],True),avg(u[:,0,1]),avg(u[:,1,1]),af,ab])


def flat_spectrum(d,phase,voltage,height):
    ns=[d['nITO'],d['n0']+d['beta']*voltage/D_DAST+1j*d['kDAST'],d['nITO'],d['nSiO2'],d['nVI'] if phase=='I' else d['nVM']]
    ds=[D_ITO1,D_DAST,D_ITO2,D_SIO2,height]; npts=len(d['w']); spec=np.zeros((npts,2,2,4))
    for di in range(2):
        a=np.ones(npts,dtype=complex); b=np.zeros(npts,dtype=complex); c=b.copy(); e=a.copy()
        for n,h in zip(ns if di==0 else ns[::-1],ds if di==0 else ds[::-1]):
            z=2*np.pi*n*h/d['w']; u=np.cos(z); v=-1j*np.sin(z)/n; t=-1j*n*np.sin(z)
            a,b,c,e=a*u+b*t,a*v+b*u,c*u+e*t,c*v+e*u
        den=a+b+c+e; R=np.abs((a+b-c-e)/den)**2; T=np.abs(2/den)**2
        spec[:,di,:,:]=np.stack([R,T,1-R-T,T],axis=-1)[:,None,:]
    return spec


def validate_flat_limit(solver,raw,order,nx):
    d=optical_data(raw,np.array([411.,557.,733.,991.,1537.,2203.]))
    worst=0.
    for phase in ['I','M']:
        target=flat_spectrum(d,phase,0.,solver.height)
        for j,lam in enumerate(d['w']):
            ns=[d['nITO'][j],d['n0'][j]+1j*d['kDAST'][j],d['nITO'][j],d['nSiO2'][j],d['nVI' if phase=='I' else 'nVM'][j]]
            actual=solver.one(lam,ns,order,nx,uniform=solver.height)
            worst=max(worst,float(np.max(np.abs(actual-target[j]))))
    print(f'Uniform-VO2 RCWA versus independent planar TMM: max |difference|={worst:.3e}',flush=True)
    if worst>2e-6: raise RuntimeError('Planar-limit validation failed.')
    return worst


def crossing_estimates(v,y,target=0.):
    roots=[]; z=y-target
    for j in range(len(v)-1):
        if z[j]*z[j+1]<0: roots.append(v[j]-z[j]*(v[j+1]-v[j])/(z[j+1]-z[j]))
    return roots  # Linear interpolation only; not exact root solves or all-root proof.


def make_figure(d,V,spectra,curves,controls,solver):
    fig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained'); axes=axs.ravel(); colors=['tab:blue','tab:red']
    ax=axes[0]; z=0.
    for name,h,color in zip(['ITO1','DAST','ITO2','SiO2'],[D_ITO1,D_DAST,D_ITO2,D_SIO2],['#7cb4c6','#daa3dc','#7cb4c6','#d6dfec']):
        ax.add_patch(Rectangle((0,z),solver.period,h,facecolor=color,edgecolor='k',lw=.6));ax.text(solver.period/2,z+h/2,name,ha='center',va='center'); z+=h
    left=(1-solver.fill)*solver.period/2
    ax.add_patch(Rectangle((left,z),solver.fill*solver.period,solver.height,facecolor='#dc9365',edgecolor='k'))
    ax.text(solver.period/2,z+solver.height/2,'VO2',ha='center',va='center'); ax.text(solver.period/2,-18,'F: air',ha='center'); ax.text(solver.period/2,z+solver.height+25,'B: air',ha='center')
    ax.set(xlim=(0,solver.period),ylim=(z+solver.height+45,-35),xlabel='x (nm), periodic',ylabel='z (nm)',title='(a) One unit cell; invariant along y')
    for pi,p in enumerate(['I','M']):
        for j,ls in [(0,'--'),(len(V)-1,'-')]:
            u=spectra[p][j].mean(axis=2); label=f'{p}, {V[j]:+.2g} V'
            axes[1].plot(d['w'],u[:,0,1],ls,color=colors[pi],label=label)
            axes[2].plot(d['w'],u[:,0,2]-u[:,1,2],ls,color=colors[pi],label=label)
        y=curves[p]; axes[3].plot(V,y[:,0],'-o',ms=3,color=colors[pi],label=p)
        axes[4].plot(V,np.minimum(y[:,1],y[:,2]),'-',color=colors[pi],label=f'{p}, total min(F,B)')
        axes[4].plot(V,np.minimum(y[:,3],y[:,4]),'--',color=colors[pi],label=f'{p}, zero-order min(F,B)')
    axes[1].set(title='(b) Forward total transmission',xlabel='Wavelength (nm)',ylabel='T, TE/TM mean',ylim=(0,1))
    axes[2].set(title='(c) Directional absorption difference',xlabel='Wavelength (nm)',ylabel='A_F - A_B');axes[2].axhline(0,color='gray',lw=.7)
    axes[3].set(title='(d) Voltage response',xlabel='DAST voltage (V)',ylabel='Solar-weighted A_F - A_B')
    for lev in [-MARGIN,0,MARGIN]:axes[3].axhline(lev,color='gray',ls=':',lw=.8)
    axes[4].set(title='(e) Luminous transmission',xlabel='DAST voltage (V)',ylabel='D65 x V(lambda) weighted T')
    labels=['Grating','Film h','Film f h']; xx=np.arange(3)
    for pi,p in enumerate(['I','M']):
        y=[controls[(p,k)][1:3].min() for k in labels]
        axes[5].bar(xx+(pi-.5)*.32,y,width=.32,color=colors[pi],label=p)
    axes[5].set(xticks=xx,xticklabels=labels,ylim=(0,1),title='(f) Zero-bias geometry controls',ylabel='Minimum luminous T of both directions')
    for ax in axes[1:]:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Air / ITO / DAST / ITO / SiO2 / VO2 grating / Air — exploratory model',fontsize=13)
    return fig


def main(a=None):
    a=arguments() if a is None else a; tc=get_torcwa(a.torcwa_path)
    device=torch.device('cuda' if (a.device=='cuda' or (a.device=='auto' and torch.cuda.is_available())) else 'cpu')
    if device.type=='cuda' and not torch.cuda.is_available(): raise RuntimeError('CUDA requested but unavailable.')
    torch.set_num_threads(1); torch.backends.cuda.matmul.allow_tf32=False
    if not 0<FILL<1 or min(D_ITO1,D_DAST,D_ITO2,D_SIO2,H_VO2,PERIOD)<=0: raise ValueError('Invalid geometry.')
    raw=load_user_data(a); w=raw['w']; d=optical_data(raw,w); solver=Solver(tc,device)
    vmax=min(VHARD,EMAX*D_DAST/1000); V=np.linspace(-vmax,vmax,a.nv)
    print(f'\nAir/ITO({D_ITO1:g})/DAST({D_DAST:g})/ITO({D_ITO2:g})/SiO2({D_SIO2:g})/VO2 grating/Air (nm)')
    print(f'Period={PERIOD:g} nm, width={PERIOD*FILL:g} nm, height={H_VO2:g} nm; same-volume film={FILL*H_VO2:g} nm')
    print(f'DAST=L1/E1 scalar model; |V|<={vmax:g} V; Emax={EMAX:g} V/um; device={device}; complex128')
    print(f'Exploratory grid: {len(w)} wavelengths, N={a.order}, Nx={a.nx}, {len(V)} biases per phase; no thermal evolution.')
    if a.quick: print('QUICK MODE: functional check only; coarse grids are NOT converged.',flush=True)
    spectra={}; curves={}; rows=[]; controls={}; audits=[]; t=time.time()
    with torch.no_grad():
        flat_error=validate_flat_limit(solver,raw,a.order,a.nx)
        for p in ['I','M']:
            spectra[p]=[]
            for j,v in enumerate(V):
                q=solver.spectrum(d,p,v,a.order,a.nx);spectra[p].append(q); y=weighted(d,q)
                rows.append(dict(Phase=p,Voltage_V=v,**dict(zip(METRICS,y))))
                print(f'{p}: {j+1}/{len(V)} V={v:+.3f}; DA={y[0]:+.6f}; Tlum_min={min(y[1:3]):.4f}; elapsed={time.time()-t:.0f}s',flush=True)
            curves[p]=np.array([weighted(d,q) for q in spectra[p]])
            y=curves[p][:,0]; mr=min(y[0],-y[-1]); zr=crossing_estimates(V,y)
            print(f'{p}: Mrev(+ to -)={mr:+.7f}; zero-crossing linear estimates={zr}; voltage grid not a proof of all roots.')
            controls[(p,'Grating')]=curves[p][len(V)//2]
            for name,h in [('Film h',H_VO2),('Film f h',FILL*H_VO2)]: controls[(p,name)]=weighted(d,flat_spectrum(d,p,0.,h))
        if not a.no_check:
            fine=optical_data(raw,wavelengths(a.step/2))
            for p in ['I','M']:
                for j in [0,len(V)//2,len(V)-1]:
                    v=V[j]; base=spectra[p][j]; high=solver.spectrum(d,p,v,a.check_order,a.nx)
                    spatial=solver.spectrum(d,p,v,a.check_order,2*a.nx)
                    spectral=solver.spectrum(fine,p,v,a.check_order,2*a.nx)
                    err_order=float(np.max(abs(high-base))); err_space=float(np.max(abs(spatial-high)))
                    mo=np.max(abs(weighted(d,high)-weighted(d,base))); ms=np.max(abs(weighted(d,spatial)-weighted(d,high)))
                    mw=np.max(abs(weighted(fine,spectral)-weighted(d,spatial)))
                    audits.append(dict(Phase=p,Voltage_V=v,OrderSpectralDiff=err_order,PixelSpectralDiff=err_space,OrderMetricDiff=mo,PixelMetricDiff=ms,WavelengthMetricDiff=mw))
                    print(f'AUDIT {p}, V={v:+.3f}: order={mo:.2e}, pixels={ms:.2e}, wavelength={mw:.2e} (weighted metrics)',flush=True)
    CurveTable=pd.DataFrame(rows); ControlTable=pd.DataFrame([dict(Phase=p,Geometry=k,**dict(zip(METRICS,v))) for (p,k),v in controls.items()])
    PrecisionTable=pd.DataFrame(audits); summaries=[]; polrows=[]
    for p in ['I','M']:
        y=curves[p]; summaries.append(dict(Phase=p,Mrev=min(y[0,0],-y[-1,0]),DAminus=y[0,0],DAplus=y[-1,0],Tlum_zero_min=min(y[len(V)//2,1:3]),Tlum_sample_min=y[:,1:3].min()))
        for j in [0,len(V)//2,len(V)-1]:
            for pi,pol in enumerate(['TE','TM']):
                z=weighted(d,spectra[p][j][:,:,pi:pi+1,:]);polrows.append(dict(Phase=p,Voltage_V=V[j],Polarization=pol,**dict(zip(METRICS,z))))
    ResultTable=pd.DataFrame(summaries);PolarizationTable=pd.DataFrame(polrows)
    pd.set_option('display.width',180);pd.set_option('display.max_columns',20)
    print('\nResultTable\n',ResultTable.to_string(index=False));print('\nZero-bias controls\n',ControlTable.to_string(index=False))
    if audits:
        print('\nPrecisionTable\n',PrecisionTable.to_string(index=False))
        if PrecisionTable[['OrderMetricDiff','PixelMetricDiff','WavelengthMetricDiff']].to_numpy().max()>1e-4 or PrecisionTable[['OrderSpectralDiff','PixelSpectralDiff']].to_numpy().max()>1e-3:
            warnings.warn('Convergence audit exceeds exploratory tolerances. Increase order/nx and reduce step; do not report these as converged results.')
    print(f'\nMaximum matched-zero-order T reciprocity error={solver.max_t0_error:.3e}; cross-polar power={solver.max_cross_power:.3e}; cutoff offsets={solver.cutoff_shifts}')
    print('Total T_F need not equal total T_B when higher diffraction orders propagate. No clipping of negative absorption is used.')
    print('Mrev tests the designated + to - endpoint orientation only. All voltage curves are sampled; no exact minimum switching swing is claimed.')
    fig=make_figure(d,V,spectra,curves,controls,solver)
    print('ONE figure; no exported files. Tables/spectra are returned by main() and retained in interactive execution.')
    if not a.no_show:plt.show()
    return dict(ResultTable=ResultTable,CurveTable=CurveTable,ControlTable=ControlTable,PrecisionTable=PrecisionTable,PolarizationTable=PolarizationTable,spectra=spectra,wavelength_nm=w,voltage_V=V,figure=fig,flat_limit_error=flat_error)


if __name__=='__main__':
    results=main(); globals().update(results)
