#!/usr/bin/env python3
"""Wu, Meng & Liu, Materials 13 (2020) 3751, Fig. 1/2-like 2D RCWA reproduction.
Paper: DOI 10.3390/ma13173751. Exact geometry: square Ti/Al2O3/Au MIM 180 x 180 nm2,
400 x 400 nm2 periodicity, thickness 30/60/40 nm, PMMA 150 nm, SiO2 half-space.
Uses 14Ti.xlsx, 11Au.xlsx, 9PMMA.xlsx, 2SiO2.xlsx from optical_constants/.
Al2O3 is Boidin et al. 2016 amorphous-film experimental n (embedded CC0 table), k=0;
NOT the unpublished exact tabulation employed by Wu et al. Paper used FDTD; this uses 2D RCWA.
PMMA geometry: 130 nm MIM embedded at top of a 150 nm PMMA coating, 20 nm PMMA below.
This is one physically plausible interpretation of paper Fig. 1, not uniquely specified.
Run via VS Code Python interpreter, no MATLAB, no exports, one 2x3 figure.
"""
from pathlib import Path
import time
T0=time.perf_counter()
import numpy as np, pandas as pd, torch, torcwa
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
# --- TOP SETTINGS (edit here; press Run Python File in VS Code) ---
flag=0                # 0=CPU, 1=CUDA; complex128 GPU is often slower on RTX 3050 Laptop
QUICK=0               # 1: exploratory 20 nm / lower Fourier order; 0: main 5 nm
STEP_NM=5.; ORDER=3; NX=NY=80   # 2D Fourier modes (2*ORDER+1)^2; convergence NOT assumed
RUN_MAPS=0            # Diagnose main spectrum first; set 1 after R+T+A checks pass
TEST_ONLY=0           # 1=three numerical self-checks only; 0=full 2D paper spectrum
MAP_STEP_NM=20.; MAP_H=(20.,30.,40.,50.,60.,70.,80.)
PERIOD=400.; WIDTH=180.; TI_H=30.; AL_H=60.; AU_H=40.; PMMA_H=150.
PROJECT=Path(__file__).resolve().parent; SHOW_FIGURE=1
DTYPE=torch.complex128; TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
# Al2O3: Boidin et al., Ceramics International 42 (2016) 1177-1182, tabulated n;
# https://refractiveindex.info/?shelf=main&book=Al2O3&page=Boidin (CC0 database).
# Values every 20 nm from 400 to 1300 nm; assume negligible k in this range.
BOIDIN_N=np.array([1.70185,1.69791,1.69454,1.69164,1.68912,1.68691,1.68496,1.68324,1.68169,1.68031,
                   1.67906,1.67792,1.67689,1.67595,1.67507,1.67427,1.67353,1.67283,1.67219,1.67158,
                   1.67101,1.67047,1.66996,1.66948,1.66902,1.66858,1.66816,1.66776,1.66737,1.66699,
                   1.66663,1.66628,1.66594,1.66561,1.66528,1.66496,1.66465,1.66435,1.66405,1.66376,
                   1.66347,1.66318,1.66290,1.66262,1.66234,1.66207])

def material_excel(file,w,lossless=False):
    path=PROJECT/'optical_constants'/file
    if not path.is_file(): raise FileNotFoundError(f'Missing {path}; place this script beside optical_constants/')
    df=pd.read_excel(path,header=None,engine='openpyxl')
    if df.shape[1]<3: raise ValueError(f'{file}: require three columns [wavelength,n,k]')
    a=df.iloc[:,:3].apply(pd.to_numeric,errors='coerce').to_numpy(dtype=float)
    a=a[np.isfinite(a).all(axis=1)];a=a[np.argsort(a[:,0])];
    if len(a)<3: raise ValueError(f'{file}: insufficient valid data')
    a[:,0]*=(1000 if np.nanmedian(a[:,0])<20 else 1)
    x,idx=np.unique(a[:,0],return_index=True);n=a[idx,1];k=a[idx,2]
    if x[0]>w[0]+1e-8 or x[-1]<w[-1]-1e-8:raise ValueError(f'{file}: {x[0]:g}–{x[-1]:g} nm does not cover {w[0]:g}–{w[-1]:g} nm')
    if np.min(n)<=0 or np.min(k)<-1e-8:raise ValueError(f'{file}: invalid n/k values')
    return np.interp(w,x,n)+1j*(0 if lossless else np.maximum(np.interp(w,x,k),0))

def read_data(w):
    d={'w':w,'Ti':material_excel('14Ti.xlsx',w),'Au':material_excel('11Au.xlsx',w),
       'PMMA':material_excel('9PMMA.xlsx',w,lossless=True),'SiO2':material_excel('2SiO2.xlsx',w,lossless=True)}
    for fname in ('Al2O3.xlsx','15Al2O3.xlsx'):
        if (PROJECT/'optical_constants'/fname).is_file():
            d['Al2O3']=material_excel(fname,w);al_source=fname+' (user dataset)';break
    else:
        d['Al2O3']=np.interp(w,np.arange(400.,1300.1,20.),BOIDIN_N).astype(complex)
        al_source='Boidin 2016 measured amorphous Al2O3 n; k set to zero (not Wu original data)'
    print('Material sources: Ti=14Ti.xlsx, Au=11Au.xlsx, PMMA=9PMMA.xlsx (lossless), SiO2=2SiO2.xlsx (lossless), Al2O3='+al_source)
    return d

def square_mask(device,nx,ny):
    # Cell-average square geometry (not 1D strip). 1D edge overlaps multiplied in x,y.
    def a(n):
        edges=np.arange(n)/n;u=(1-WIDTH/PERIOD)/2;v=(1+WIDTH/PERIOD)/2
        return np.clip(np.minimum(edges+1/n,v)-np.maximum(edges,u),0,None)*n
    return torch.as_tensor(a(nx)[:,None]*a(ny)[None,:],dtype=DTYPE,device=device)

def channel_power(s,lam,order,nsub,direction,port,device):
    """Sum *orthogonal* p/s power-normalized diffraction channels.
    At normal incidence with azimuth=0, incident p corresponds to x polarization.
    For arbitrary 2D orders, pp/sp are orthogonal outgoing plane-wave polarizations;
    summing xx+yx is not generally a valid 2D power flux decomposition.
    """
    nout=(1. if direction=='forward' else nsub) if port=='reflection' else (nsub if direction=='forward' else 1.)
    m=range(-order,order+1); orders=np.array([(i,j) for i in m for j in m if (i*i+j*j)*(lam/PERIOD)**2 < nout*nout-1e-10],dtype=int).reshape(-1,2)
    if not len(orders):raise ValueError('No propagating orders')
    total=torch.zeros(len(orders),device=device,dtype=torch.float64)
    for pol in ('pp','sp'):
        amp=s.S_parameters(orders=orders,direction=direction,port=port,polarization=pol,ref_order=[0,0],power_norm=True)
        total+=amp.abs().square().reshape(-1).real
    values=total.detach().cpu().numpy();zero=np.flatnonzero((orders==[0,0]).all(axis=1))
    if not np.isfinite(values).all():raise FloatingPointError(f'Nonfinite {direction} {port} power at {lam:g} nm')
    return float(values.sum()),float(values[zero[0]])

def solve(tc,d,j,h,dev,mask,order):
    lam=float(d['w'][j]);nsub=float(d['SiO2'][j].real)
    if abs(lam-PERIOD)<1e-9:lam+=1e-5  # avoid air first-order Rayleigh singularity at 400 nm
    s=tc.rcwa(freq=1/lam,order=[order,order],L=[PERIOD,PERIOD],dtype=DTYPE,device=dev,stable_eig_grad=False)
    s.add_output_layer(eps=nsub**2);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
    pm=complex(d['PMMA'][j]**2)
    for name,thick in [('Ti',TI_H),('Al2O3',h),('Au',AU_H)]:
        eps=pm+mask*(complex(d[name][j]**2)-pm)
        s.add_layer(thickness=thick,eps=eps)
    leftover=PMMA_H-(TI_H+h+AU_H)
    if leftover<0:raise ValueError('PMMA_H must >= TI_H + Al2O3 thickness + AU_H')
    if leftover>0:s.add_layer(thickness=leftover,eps=pm)
    s.solve_global_smatrix()
    result=np.zeros((2,4))
    for di,direction in enumerate(('forward','backward')):
        result[di,0],_=channel_power(s,lam,order,nsub,direction,'reflection',dev)
        result[di,1],result[di,3]=channel_power(s,lam,order,nsub,direction,'transmission',dev)
        result[di,2]=1-result[di,0]-result[di,1]
    if np.min(result)<-2e-4 or np.max(result)>1+2e-4 or not np.isfinite(result).all():
        q2=np.array([a*a+b*b for a in range(-order,order+1) for b in range(-order,order+1)])
        anomaly=np.min(np.abs(nsub*nsub-q2*(lam/PERIOD)**2))
        raise FloatingPointError(f'Passivity failed at {lam:g} nm, N={order}, Nx={NX}: {result}; cutoff-distance={anomaly:.3g}. '+
                                   'Using pp/sp power-normalized channels; investigate material signs and Fourier convergence, do not clip absorption.')
    return result

def validate_interface(tc,d,dev,order):
    """Air / (12 nm air spacer) / SiO2 Fresnel check; explicit film avoids empty-stack ambiguity."""
    worst=0.; print('CHECK 1: analytic Fresnel vs TORCWA (explicit 12 nm lossless air spacer)')
    with torch.no_grad():
        for lam in (405.,415.,500.,890.,1300.):
            n=float(np.interp(lam,d['w'],d['SiO2'].real));s=tc.rcwa(freq=1/lam,order=[order,order],L=[PERIOD,PERIOD],dtype=DTYPE,device=dev,stable_eig_grad=False)
            s.add_output_layer(eps=n*n);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
            s.add_layer(thickness=12.,eps=1.0);s.solve_global_smatrix()
            exactR=((n-1)/(n+1))**2;exactT=4*n/(n+1)**2
            for direction in ('forward','backward'):
                r,_=channel_power(s,lam,order,n,direction,'reflection',dev)
                t,_=channel_power(s,lam,order,n,direction,'transmission',dev)
                delta=max(abs(r-exactR),abs(t-exactT),abs(r+t-1));worst=max(worst,delta)
                if delta>1e-6:
                    def xy_power(port):
                        return sum(float(s.S_parameters(orders=[0,0],direction=direction,port=port,polarization=pol,ref_order=[0,0],power_norm=True).abs().square().sum().item()) for pol in ('xx','yx'))
                    print(f'    diagnostic xy zero-order: R={xy_power("reflection"):.7f}, T={xy_power("transmission"):.7f}',flush=True)
                print(f'  {lam:6.1f} nm {direction:8s}: R={r:.7f}, T={t:.7f}; exact R={exactR:.7f}, T={exactT:.7f}; err={delta:.2e}',flush=True)
    print(f'Fresnel max error={worst:.3e} (target <=1e-6)',flush=True)
    if worst>1e-6:raise RuntimeError('Fresnel failed: power-basis, port or installed TORCWA convention requires investigation. No paper spectra computed.')

def validate_lossless_pattern(tc,d,dev,order,mask):
    """Check R+T=1 for a lossless 2D square grating and a SiO2 substrate, all open orders included."""
    worst=0.;print('CHECK 2: lossless patterned 2D grating / SiO2 energy conservation')
    with torch.no_grad():
        for lam in (405.,415.,500.,890.,1300.):
            n=float(np.interp(lam,d['w'],d['SiO2'].real));s=tc.rcwa(freq=1/lam,order=[order,order],L=[PERIOD,PERIOD],dtype=DTYPE,device=dev,stable_eig_grad=False)
            s.add_output_layer(eps=n*n);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
            s.add_layer(thickness=60.,eps=(1.49**2+mask*(1.72**2-1.49**2)));s.solve_global_smatrix()
            for direction in ('forward','backward'):
                r,_=channel_power(s,lam,order,n,direction,'reflection',dev)
                t,_=channel_power(s,lam,order,n,direction,'transmission',dev)
                err=abs(r+t-1);worst=max(worst,err)
                print(f'  {lam:6.1f} nm {direction:8s}: R={r:.7f}, T={t:.7f}, R+T-1={r+t-1:+.3e}',flush=True)
    print(f'Lossless 2D max energy error={worst:.3e} (target <=1e-3)',flush=True)
    if worst>1e-3:raise RuntimeError('Lossless 2D grating failed energy conservation: check outgoing power and order truncation before MIM calculation.')

def run_spectrum(d,h,device,mask,order,label):
    start=time.perf_counter();res=[]
    with torch.no_grad():
        for j in range(len(d['w'])):
            res.append(solve(torcwa,d,j,h,device,mask,order))
            if (j+1)%max(10,len(d['w'])//4)==0:
                if flag:torch.cuda.synchronize()
                print(f' {label}: {j+1}/{len(d["w"])} wavelengths; elapsed {time.perf_counter()-start:.1f}s',flush=True)
    if flag:torch.cuda.synchronize()
    q=np.array(res);print(f' {label}: completed {len(q)} wavelengths in {time.perf_counter()-start:.1f}s',flush=True)
    return q

def draw(d,q,dh,hs,afmaps,rbmaps):
    w=d['w'];A=q[:,:,2].T;fig,axes=plt.subplots(2,3,figsize=(16.8,9.2),layout='constrained');ax=axes.ravel()
    y=0.;xl=(PERIOD-WIDTH)/2
    for name,h,c in [('Ti',TI_H,'#aaaab0'),('Al2O3',AL_H,'#dba4ad'),('Au',AU_H,'#e9cc79')]:
        ax[0].add_patch(Rectangle((xl,y),WIDTH,h,facecolor=c,edgecolor='k',lw=.6));ax[0].text(PERIOD/2,y+h/2,name,ha='center',va='center',fontsize=9);y+=h
    remaining=PMMA_H-y
    ax[0].add_patch(Rectangle((0,0),xl,PMMA_H,facecolor='#dcefed',alpha=.9))
    ax[0].add_patch(Rectangle((xl+WIDTH,0),xl,PMMA_H,facecolor='#dcefed',alpha=.9))
    if remaining>0:ax[0].add_patch(Rectangle((xl,y),WIDTH,remaining,facecolor='#dcefed',edgecolor='k',lw=.5))
    ax[0].axhline(PMMA_H,color='k',lw=.6);ax[0].text(PERIOD/2,PMMA_H+13,'SiO2 substrate',ha='center',fontsize=9)
    ax[0].set(xlim=(0,PERIOD),ylim=(PMMA_H+26,-16),title='(a) Cross-section; square in x-y',xlabel='x (nm)',ylabel='Depth (nm)')
    for a,di,name in [(ax[1],0,'Front: air → Ti'),(ax[2],1,'Back: SiO2 → Au')]:
        for j,k in enumerate(['R','T','A']):a.plot(w,q[:,di,j],label=k,lw=1.5)
        a.set(title=f'({"b" if di==0 else "c"}) {name}',xlabel='Wavelength (nm)',ylabel='Power fraction',ylim=(-.02,1.02));a.legend()
    if afmaps is not None:
        for k,(data,title) in enumerate([(afmaps,'Front absorptance A'),(rbmaps,'Back reflectance R')]):
            im=ax[3+k].imshow(data,aspect='auto',origin='lower',extent=(dh[0],dh[-1],hs[0],hs[-1]),vmin=0,vmax=1,cmap='inferno');fig.colorbar(im,ax=ax[3+k],fraction=.05)
            ax[3+k].set(title=f'({"d" if k==0 else "e"}) {title} vs Al2O3 thickness',xlabel='Wavelength (nm)',ylabel='Al2O3 thickness (nm)')
    else:
        for j in (3,4):ax[j].text(.5,.5,'Set RUN_MAPS=1 for thickness map',ha='center',va='center',transform=ax[j].transAxes);ax[j].set_axis_off()
    ax[5].plot(w,A[0]-A[1],label=r'$A_F-A_B$',color='darkred')
    ax[5].axhline(0,color='gray',lw=.8);ax[5].set(title='(f) Directional absorption difference',xlabel='Wavelength (nm)',ylabel=r'$A_F-A_B$');ax[5].legend()
    for a in ax:
        if a.axison:a.grid(alpha=.13)
    fig.suptitle('Wu et al. (2020) Ti / Al2O3 / Au, 2D square-array RCWA reproduction attempt',fontsize=13)
    return fig

def main():
    if flag not in (0,1):raise ValueError('flag must be 0 or 1')
    if flag and not torch.cuda.is_available():raise RuntimeError('CUDA unavailable: use flag=0')
    if QUICK: step,order=20.,2;mapstep=40.
    else:step,order=STEP_NM,ORDER;mapstep=MAP_STEP_NM
    if (2*order+1)**2>NX*NY or min(NX,NY)<4*order+4:raise ValueError('Increase Nx/Ny or reduce Fourier order')
    if PMMA_H<TI_H+max((AL_H,*MAP_H) if RUN_MAPS else (AL_H,))+AU_H:raise ValueError('Increase PMMA_H or reduce maximum map thickness for embedded model')
    device=torch.device('cuda' if flag else 'cpu');torch.set_num_threads(1)
    if flag:torch.backends.cuda.matmul.allow_tf32=False
    w=np.arange(400.,1300.01,step);w=np.unique(np.r_[w,1300.]);d=read_data(w);mask=square_mask(device,NX,NY)
    print(f'\nWu 2020 2D square MIM reproduction; device={device}; PyTorch={torch.__version__}; 2D N={order}, modes={(2*order+1)**2}; Nx×Ny={NX}×{NY}; λ={step:g} nm')
    print(f'P={PERIOD:g}, square={WIDTH:g}, Ti={TI_H:g}, Al2O3={AL_H:g}, Au={AU_H:g}, PMMA={PMMA_H:g} nm; maps={bool(RUN_MAPS)}',flush=True)
    validate_interface(torcwa,d,device,order)
    validate_lossless_pattern(torcwa,d,device,order,mask)
    print('CHECK 3: absorptive Ti/Al2O3/Au structure, passivity must hold')
    if TEST_ONLY:
        print('Two preliminary checks PASSED. TEST_ONLY=1: set TEST_ONLY=0 for full spectrum; no result figure yet.',flush=True)
        return None,None
    q=run_spectrum(d,AL_H,device,mask,order,'BASE')
    v=(w>=400)&(w<=750);visible_avg=TRAPZ(q[v,0,1],w[v])/(w[v][-1]-w[v][0])
    jaf=np.argmax(q[:,0,2]);jrb=np.argmax(q[:,1,0]);k=np.argmin(abs(w-900))
    print(f'Front peak A={q[jaf,0,2]:.4f} at {w[jaf]:g} nm; back peak R={q[jrb,1,0]:.4f} at {w[jrb]:g} nm')
    print(f'At ~{w[k]:g} nm: front [R,T,A]={q[k,0]}, back [R,T,A]={q[k,1]}')
    print(f'Visible (400-750 nm) unweighted front average T={visible_avg:.4f}; zero-order reciprocity max |T0F-T0B|={np.max(np.abs(q[:,0,3]-q[:,1,3])):.2e}')
    print('At normal incidence, the matched zero-order transmittance must agree for reciprocal materials. Side-dependent TOTAL power accounts for all propagating substrate orders.')
    print('Paper benchmark (NOT computed): Ti/Al2O3/Au front NIR absorption >0.90 near ~890nm, back reflectance >0.70 near ~910nm, average visible transmission >0.60.')
    afmaps=rbmaps=hs=None;dh=None
    if RUN_MAPS:
        dh=np.arange(400.,1300.01,mapstep);dh=np.unique(np.r_[dh,1300.]);ds=read_data(dh)
        hs=np.array(MAP_H);afmaps=np.zeros((len(hs),len(dh)));rbmaps=afmaps.copy()
        for i,h in enumerate(hs):
            result=run_spectrum(ds,float(h),device,mask,order,f'MAP Al2O3={h:g}nm')
            afmaps[i]=result[:,0,2];rbmaps[i]=result[:,1,0]
    fig=draw(d,q,dh,hs,afmaps,rbmaps)
    print(f'ALL startup + solve + figures: {time.perf_counter()-T0:.1f}s ({(time.perf_counter()-T0)/60:.2f}min); 1 figure, NO exports.',flush=True)
    if SHOW_FIGURE:plt.show()
    return q,fig

if __name__=='__main__':
    Spectra,Figure=main()
