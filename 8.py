#!/usr/bin/env python3
"""DAST/VO2 voltage-and-phase tunable asymmetric 1D MIM-inspired metasurface.
Place alongside optical_constants/, VO2_temperature_data/, Atmospheric data/.
Run directly in VS Code using your xhh31132 Python interpreter. No MATLAB, no exports.
Ref: Wu, Meng, Liu, Materials 13 (2020) 3751; this is a NEW 1D candidate,
not a reproduction of their 2D Ti/Al2O3/Au structure or a verified device.
"""
from pathlib import Path
import time, warnings
T_START=time.perf_counter()
import numpy as np, pandas as pd, torch, torcwa
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ===== EDIT HERE; save and click VS Code 'Run Python File' =====
flag=0                     # 0 CPU; 1 CUDA (RTX 3050 may be slower at complex128)
QUICK=0                    # 1 small functionality test; 0 = 5 nm candidate sweep
STEP_NM=5.0; ORDER=9; NX=256
PERIOD=400.; FILL=.45       # 1D strip width = 180 nm, not the 2D 180x180 square of Wu et al.
H_PMMA=150.; H_VO2=20.; H_ITO1=8.; H_DAST=85.; H_ITO2=8.; H_AU=40.
USE_TI=1; H_TI=30.         # 1 uses optical_constants/14Ti.xlsx; 0 removes Ti (comparison control)
VOLTAGES=(-8.5,0.,8.5)     # V across DAST, scalar L1/E1 model
EO_SIGN=1; PHASES=('I','M'); FIELD_LIMIT=100.  # V/um, design screening only
PROJECT=Path(__file__).resolve().parent
SHOW_FIGURE=1; DTYPE=torch.complex128
TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz


def read_table(path,cols,unit='nm',outside='error',min_rows=3):
    if not path.is_file(): raise FileNotFoundError(f'Missing original material/weight file: {path}')
    if path.suffix.lower()=='.csv': df=pd.read_csv(path,header=None,comment='#',on_bad_lines='skip')
    else: df=pd.read_excel(path,header=None,engine='openpyxl')
    if df.shape[1]<=max(cols): raise ValueError(f'{path.name}: expected at least {max(cols)+1} columns')
    a=df.iloc[:,list(cols)].apply(pd.to_numeric,errors='coerce').to_numpy(dtype=float)
    a=a[np.all(np.isfinite(a),axis=1)];a=a[np.argsort(a[:,0])]
    a=a[np.unique(a[:,0],return_index=True)[1]]
    if len(a)<min_rows: raise ValueError(f'{path.name}: too few numeric data rows')
    a[:,0]*=(1000. if unit=='um' else 1.)
    return a


def interpolate(path,w,cols,unit='nm',outside='error',min_rows=3):
    a=read_table(path,cols,unit,outside,min_rows)
    if outside=='error' and (a[0,0]>w[0]+1e-7 or a[-1,0]<w[-1]-1e-7):
        raise ValueError(f'{path.name}: data cover {a[0,0]:g}--{a[-1,0]:g} nm, but simulation needs {w[0]:g}--{w[-1]:g} nm')
    return np.stack([np.interp(w,a[:,0],a[:,j],left=0 if outside=='zero' else None,
                               right=0 if outside=='zero' else None) for j in range(1,a.shape[1])],axis=1)


def load_data():
    t=time.perf_counter();root=PROJECT
    if not (root/'optical_constants'/'11Au.xlsx').is_file():
        raise FileNotFoundError(f'Place this script next to optical_constants/, VO2_temperature_data/ and Atmospheric data/: {root}')
    w=np.unique(np.r_[np.arange(400,2500+1e-9,STEP_NM if not QUICK else 100.),400.,430.,457.9,488.,514.5,
                        650.,700.,720.,750.,780.,900.,1550.,2100.,2500.]);w.sort()
    m=root/'optical_constants'; a=root/'Atmospheric data'; s=root/'VO2_temperature_data'
    names={'ITO':'4TIO.xlsx','SiO2':'2SiO2.xlsx','VI':'6VO2_I.xlsx','VM':'7VO2_M.xlsx',
           'PMMA':'9PMMA.xlsx','Au':'11Au.xlsx'}
    if USE_TI: names['Ti']='14Ti.xlsx' # Dataset shown in project: wavelength 0.28--10 um, n, k
    d={'w':w}
    for name,file in names.items():
        nk=interpolate(m/file,w,(0,1,2),unit='um')
        if np.any(nk[:,0]<=0) or np.any(nk[:,1]<-1e-10): raise ValueError(f'{file}: invalid passive n,k')
        d[name]=nk[:,0]+1j*np.maximum(nk[:,1],0) # TORCWA: exp(-i omega t), passive n+i*k
    d['solar']=interpolate(s/'ASTMG173.xlsx',w,(0,2),outside='zero')[:,0]
    vf=a/'CIE_Vlambda.xlsx'
    if not vf.exists(): vf=s/'CIE_Vlambda.csv'
    if not vf.exists(): raise FileNotFoundError('Missing CIE_Vlambda.xlsx or CIE_Vlambda.csv')
    src=pd.read_excel(vf,header=None) if vf.suffix.lower()=='.xlsx' else pd.read_csv(vf,header=None)
    col=2 if src.shape[1]>=4 else 1 if src.shape[1]==2 else None
    if col is None:raise ValueError('CIE_Vlambda expects [nm,V] or [nm,xbar,ybar,zbar]')
    v=interpolate(vf,w,(0,col),outside='zero',min_rows=100)[:,0]
    file=a/'CIE_std_illum_D65.csv'
    if not file.exists(): file=a/'CIE_std_illum_D65.xlsx'
    d['lum']=v*interpolate(file,w,(0,1),outside='zero',min_rows=100)[:,0]
    if TRAPZ(d['solar'],w)<=0 or TRAPZ(d['lum'][(w>=400)&(w<=780)],w[(w>=400)&(w<=780)])<=0:
        raise ValueError('Bad spectral weighting data')
    # Retain existing Janus.py / previous TORCWA scalar DAST L1/E1 hybrid model.
    n=np.full(len(w),2.2);j=w>650;x=w[j]; ns=np.sqrt(2.078**2+1.645*533**2/(x*x-533**2))
    blend=.5-.5*np.cos(np.pi*np.clip((x-650)/100,0,1));n[j]=(1-blend)*2.2+blend*ns
    alpha=np.full(len(w),1.9e5);j=(w>430)&(w<=514.5)
    alpha[j]=10**PchipInterpolator([430,457.9,488,514.5],np.log10([1.9e5,4.1e5,4.9e5,4.79e5]))(w[j])
    j=(w>514.5)&(w<=900);alpha[j]=10**np.interp(w[j],[514.5,700,750,900],np.log10([4.79e5,5,1.8,.8]))
    j=(w>900)&(w<=2100);xx=w[j]
    alpha[j]=.7+np.exp(-.5*((xx-1100)/55)**2)+1.5*np.exp(-.5*((xx-1400)/70)**2)+3*np.exp(-.5*((xx-1700)/95)**2)
    alpha[w>2100]=.7
    r=PchipInterpolator([720,750,1550],[530,445,200])(np.clip(w,720,1550))
    d['DAST_n']=n;d['DAST_k']=alpha*w*1e-7/(4*np.pi)
    d['DAST_beta']=EO_SIGN*.5*n**3*r*1e-3  # refr-index correction beta*V/H_DAST, H in nm
    print(f'Data: {root}; {len(w)} wavelengths; load {time.perf_counter()-t:.2f} s')
    print('Material map: ITO=4TIO.xlsx (verify material identity); Ti=14Ti.xlsx ('+('ON' if USE_TI else 'OFF')+').')
    if USE_TI:
        for nm in (550.,900.,1550.):
            nkti=np.interp(nm,w,d['Ti'].real)+1j*np.interp(nm,w,d['Ti'].imag)
            print(f' Ti at {nm:g} nm: n={nkti.real:.4g}, k={nkti.imag:.4g} (interpolated)')
    return d


def get_mask(nx,device):
    x=np.arange(nx)/nx;left=x;right=x+1/nx;lo=(1-FILL)/2;hi=(1+FILL)/2
    rho=np.maximum(0,np.minimum(right,hi)-np.maximum(left,lo))*nx
    return torch.as_tensor(rho[:,None],dtype=DTYPE,device=device)


def solve_point(tc,d,j,phase,bias,device,mask):
    lam=float(d['w'][j]);nsub=float(d['SiO2'][j].real)
    # Semi-infinite SiO2 substrate is approximated as lossless in the external port.
    if abs(d['SiO2'][j].imag)>2e-3: raise ValueError('Significant SiO2 substrate absorption: lossless external-port approximation invalid.')
    nDA=d['DAST_n'][j]+d['DAST_beta'][j]*bias/H_DAST+1j*d['DAST_k'][j]
    if nDA.real<=0:raise ValueError('DAST index became nonphysical under applied bias')
    s=tc.rcwa(freq=1/lam,order=[ORDER,0],L=[PERIOD,PERIOD],dtype=DTYPE,device=device,stable_eig_grad=False)
    s.add_output_layer(eps=nsub*nsub);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
    pm=complex(d['PMMA'][j]**2);s.add_layer(thickness=H_PMMA,eps=pm)
    layers=[]
    if USE_TI:layers.append(('Ti',H_TI,d['Ti'][j]))
    layers += [('VO2',H_VO2,d['VI'][j] if phase=='I' else d['VM'][j]),
               ('ITO1',H_ITO1,d['ITO'][j]),('DAST',H_DAST,nDA),('ITO2',H_ITO2,d['ITO'][j]),('Au',H_AU,d['Au'][j])]
    for _,h,nmat in layers:s.add_layer(thickness=h,eps=pm+mask*(complex(nmat*nmat)-pm))
    s.solve_global_smatrix();result=np.zeros((2,2,4)) # direction F/B, TE/TM, R,T,A,T0
    for di,direction in enumerate(('forward','backward')):
        for pi,inc in enumerate(('y','x')):
            for qi,port in enumerate(('reflection','transmission')):
                nout=(1. if (direction=='forward')==(port=='reflection') else nsub)
                mm=[m for m in range(-ORDER,ORDER+1) if abs(m*lam/PERIOD)<nout-1e-7]
                ords=np.array([[m,0] for m in mm],dtype=np.int64)
                sums=torch.zeros(len(mm),dtype=torch.float64,device=device)
                for out in ('x','y'):
                    amp=s.S_parameters(orders=ords,direction=direction,port=port,polarization=out+inc,ref_order=[0,0],power_norm=True)
                    sums+=amp.abs().square().reshape(-1)
                vals=sums.detach().cpu().numpy();result[di,pi,qi]=vals.sum()
                if qi==1: result[di,pi,3]=vals[mm.index(0)]
            result[di,pi,2]=1-result[di,pi,0]-result[di,pi,1]
    if not np.all(np.isfinite(result)) or np.min(result)<-5e-5 or np.max(result)>1+5e-5:
        raise FloatingPointError(f'Nonpassive or nonfinite RCWA result at {lam:g} nm, {phase}, {bias:+g} V: {result}')
    return result


def weighted(d,q):
    w=d['w'];vis=(w>=400)&(w<=780);u=q.mean(axis=2) # unpolarized TE/TM average
    def avg(x,weight,sel):return float(TRAPZ(x[sel]*weight[sel],w[sel])/TRAPZ(weight[sel],w[sel]))
    allw=np.ones(len(w),dtype=bool);fs=avg(u[:,0,2],d['solar'],allw);bs=avg(u[:,1,2],d['solar'],allw)
    lf=avg(u[:,0,1],d['lum'],vis);lb=avg(u[:,1,1],d['lum'],vis)
    return dict(DA=fs-bs,AsF=fs,AsB=bs,TlumF=lf,TlumB=lb,Tlum=min(lf,lb),
                TsolF=avg(u[:,0,1],d['solar'],allw),TsolB=avg(u[:,1,1],d['solar'],allw))


def one_figure(d,spectra,stats,biases):
    fig,axs=plt.subplots(2,4,figsize=(20,9),layout='constrained');ax=axs.ravel();z=0
    def layer(name,h,col,pattern=False):
        nonlocal z
        left=(PERIOD*(1-FILL))/2 if pattern else 0;wid=PERIOD*FILL if pattern else PERIOD
        ax[0].add_patch(Rectangle((left,z),wid,h,facecolor=col,edgecolor='black',lw=.6))
        ax[0].text(PERIOD/2,z+h/2,name,ha='center',va='center',fontsize=9);z+=h
    layer('PMMA',H_PMMA,'#dceeea')
    if USE_TI:layer('Ti',H_TI,'#b5b5ba',True)
    for name,h,col in [('VO2',H_VO2,'#d69b76'),('ITO1',H_ITO1,'#82bace'),
                       ('DAST',H_DAST,'#cfafe1'),('ITO2',H_ITO2,'#82bace'),('Au',H_AU,'#efd179')]:layer(name,h,col,True)
    ax[0].text(PERIOD/2,-16,'F: air ↓',ha='center');ax[0].text(PERIOD/2,z+16,'B: SiO2 ↑',ha='center')
    ax[0].set(xlim=(0,PERIOD),ylim=(z+35,-32),xlabel='x (nm)',ylabel='z (nm)',title='(a) Asymmetric 1D patterned unit')
    w=d['w'];base=len(biases)//2
    for phase,col in zip(PHASES,('tab:blue','tab:red')):
        for k,ls in ((0,'--'),(len(biases)-1,'-')):
            u=spectra[phase][k].mean(axis=2)
            ax[1].plot(w,u[:,0,2]-u[:,1,2],ls,color=col,label=f'{phase}, {biases[k]:+g} V')
        u=spectra[phase][base].mean(axis=2)
        ax[2].plot(w,u[:,0,2],'-',color=col,label=f'{phase}, F')
        ax[2].plot(w,u[:,1,2],'--',color=col,label=f'{phase}, B')
        ax[3].plot(w,u[:,0,1],'-',color=col,label=f'{phase}, F')
        ax[3].plot(w,u[:,1,1],'--',color=col,label=f'{phase}, B')
        rows=[r for r in stats if r['Phase']==phase]
        ax[4].plot(biases,[r['DA'] for r in rows],'-o',color=col,label=phase)
        ax[5].plot(biases,[r['Tlum'] for r in rows],'-o',color=col,label=phase)
        ax[6].plot(biases,[r['DA']-next(t['DA'] for t in stats if t['Phase']=='I' and t['Voltage']==r['Voltage']) for r in rows],'-o',color=col,label=phase)
    ax[1].axhline(0,color='gray',lw=.8)
    ax[1].set(title='(b) Spectral directional absorption',xlabel='Wavelength (nm)',ylabel='A_F - A_B')
    ax[2].set(title='(c) Absorption at 0 V',xlabel='Wavelength (nm)',ylabel='A',ylim=(0,1))
    ax[3].set(title='(d) Transmission at 0 V',xlabel='Wavelength (nm)',ylabel='T',ylim=(0,1))
    ax[4].axhline(0,color='gray',lw=.8);ax[4].set(title='(e) Voltage tuning of solar DA',xlabel='DAST voltage (V)',ylabel='Solar-weighted A_F - A_B')
    ax[5].set(title='(f) Luminous transmission',xlabel='DAST voltage (V)',ylabel='D65 × V(λ) weighted min(T_F,T_B)',ylim=(0,1))
    ax[6].axhline(0,color='gray',lw=.8);ax[6].set(title='(g) Thermal shift from insulating VO2',xlabel='DAST voltage (V)',ylabel='DA_phase - DA_insulating')
    yy=np.arange(2);bw=.34
    for k,phase in enumerate(PHASES):
        rr=[r for r in stats if r['Phase']==phase];ax[7].bar(yy+(k-.5)*bw,[rr[0]['DA'],rr[-1]['DA']],bw,label=phase)
    ax[7].axhline(0,color='gray',lw=.8);ax[7].set(xticks=yy,xticklabels=[f'{biases[0]:+g} V',f'{biases[-1]:+g} V'],title='(h) Directional solar absorption endpoints',ylabel='A_sol,F - A_sol,B')
    for a in ax[1:]:
        a.grid(alpha=.2)
        if a.get_legend_handles_labels()[1]:a.legend(fontsize=7)
    fig.suptitle(('Ti / ' if USE_TI else '')+'VO2 / ITO / DAST / ITO / Au Janus grating — preliminary 1D RCWA',fontsize=13)
    return fig


def main():
    if flag not in (0,1):raise ValueError('flag must be 0 (CPU) or 1 (CUDA)')
    if not (0<FILL<1 and NX>4*ORDER and STEP_NM>0):raise ValueError('Invalid periodic geometry or grid')
    if any(h<=0 for h in (PERIOD,H_PMMA,H_VO2,H_ITO1,H_DAST,H_ITO2,H_AU,*([H_TI] if USE_TI else []))):raise ValueError('Thicknesses must be positive')
    if max(abs(v) for v in VOLTAGES)*1000/H_DAST>FIELD_LIMIT:raise ValueError('Bias exceeds selected DAST field screening limit')
    device=torch.device('cuda' if flag else 'cpu')
    if flag and not torch.cuda.is_available():raise RuntimeError('CUDA unavailable; set flag=0')
    torch.set_num_threads(1);torch.backends.cuda.matmul.allow_tf32=False
    d=load_data();mask=get_mask(NX,device);specs={};stats=[];t=time.perf_counter()
    print(f'DEVICE={device}; CUDA available={torch.cuda.is_available()}, complex128; TORCWA N={ORDER}, Nx={NX}, lambda step={STEP_NM:g} nm')
    print(f'1D unit: PMMA / '+('Ti / ' if USE_TI else '')+'VO2 / ITO / DAST / ITO / Au / semi-infinite SiO2')
    print(f'P={PERIOD:g} nm, width={PERIOD*FILL:g} nm; voltages={VOLTAGES}; phases={PHASES}')
    print('All quoted optical results are NEW numerical predictions, not Wu et al. published performance.')
    with torch.no_grad():
        for p in PHASES:
            specs[p]=[]
            for v in VOLTAGES:
                st=time.perf_counter();q=np.stack([solve_point(torcwa,d,j,p,v,device,mask) for j in range(len(d['w']))]);specs[p].append(q)
                row=dict(Phase=p,Voltage=float(v),**weighted(d,q));stats.append(row)
                print(f'{p}, {v:+g} V: ΔA_sol={row["DA"]:+.6f}, A_F={row["AsF"]:.4f}, A_B={row["AsB"]:.4f}, Tlum={row["Tlum"]:.4f}; {time.perf_counter()-st:.1f}s',flush=True)
    if flag:torch.cuda.synchronize()
    recerr=max(float(np.max(np.abs(q[:,0,:,3]-q[:,1,:,3]))) for group in specs.values() for q in group)
    calc=time.perf_counter()-t;df=pd.DataFrame(stats)
    print('\nRESULTS\n'+df.to_string(index=False,float_format=lambda x:f'{x:.6f}'))
    for p in PHASES:
        r=df[df.Phase==p];da=r.DA.to_numpy();print(f'{p}: endpoint sign switch={bool(da[0]*da[-1]<0)}, Mrev(+ to -)={min(da[0],-da[-1]):+.6f}')
    print(f'Matched zero-order reciprocity max |T0_F-T0_B| = {recerr:.3e}; total T may differ when diffraction channels propagate.')
    print(f'Runtime (RCWA sweep): {calc:.2f}s ({calc/60:.2f} min); including import/data so far: {time.perf_counter()-T_START:.2f}s')
    fig=one_figure(d,specs,stats,VOLTAGES)
    print(f'Full startup + simulation + figure: {time.perf_counter()-T_START:.2f}s. One figure, no exports.')
    if SHOW_FIGURE:plt.show()
    return df,specs,fig

if __name__=='__main__':
    ResultTable,Spectra,Figure=main()
