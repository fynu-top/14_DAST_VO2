
#!/usr/bin/env python3
"""DAST + VO2 2D metasurface + Janus: exploratory, electrically/thermally tunable solar optics.
Inspired by Wang et al., Results in Materials 17 (2023) 100363 (thin asymmetric MIM,
FP cavity, metal duty-cycle vs transparency), and Wu et al., Materials 13 (2020) 3751.
NEW candidate, NOT a reproduced or already optimized device. Fixed VO2 insulating/metallic
states, scalar DAST EO/loss hybrid. Uses local original Excel/CIE/ASTM tables, no MATLAB.
Save beside optical_constants/, VO2_temperature_data/, Atmospheric data/ and click Run.
One 2x4 plot, no file/image/data exports. Dependencies: torch,torcwa,numpy,scipy,pandas,
openpyxl,matplotlib. 2023 paper's Ag/SiO2/Ti baseline is NOT claimed for this Au design.
"""
from pathlib import Path
import time
START=time.perf_counter()
import numpy as np, pandas as pd, torch, torcwa
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
# ========== USER SETTINGS (edit here, then VS Code: Run Python File) ==========
flag=0                 # 0: CPU; 1: CUDA; complex128 may be faster on CPU for small RCWA
QUICK=0                # 1: coarse functional test; 0: 5-nm working grid
SCAN=0                 # 1: ALSO compare three candidate VO2 patch sizes (slow)
RUN_CONTROL=1          # 1: calculate zero-bias reference without VO2 (one extra spectrum)
STEP_NM=5.; ORDER=3; NX=NY=80          # 2D order=(N,N), (2N+1)^2 modes
P=400.; W_VO2=180.; H_VO2=12.        # Square VO2 nanoislands: areal fill=(W/P)^2
H_PMMA_TOP=15.; H_PMMA_BOTTOM=23.     # PMMA total = 15+12+23=50 nm
H_TI=3.; H_ITO1=5.; H_DAST=95.; H_ITO2=5.; H_SIO2=25.; H_AU=12.
VOLTAGES=(-8.5,0.,8.5); PHASES=('I','M')
EO_SIGN=1.; E_LIMIT=100.               # V/um, scalar DAST L1/E1 exploratory
ITO_FILE='4TIO.xlsx'                 # IMPORTANT: verify if this is really ITO, not TiO!
DAST_SOURCE='hybrid'                 # 'hybrid'=previous scalar model; 'xlsx'=10DAST.xlsx
PROJECT=Path(__file__).resolve().parent
SHOW_FIGURE=1; DTYPE=torch.complex128
TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
# Optional search: three designs at 20 nm, N=2. Full-resolution plots always show base settings.
SEARCH_CANDIDATES=((120.,8.),(180.,12.),(240.,16.))

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
    w=np.arange(400.,2500.00001,step);w=np.unique(np.r_[w,2500.,457.9,488.,514.5])
    d={'w':w};names={'SiO2':'2SiO2.xlsx','ITO':ITO_FILE,'VI':'6VO2_I.xlsx','VM':'7VO2_M.xlsx',
                     'PMMA':'9PMMA.xlsx','Au':'11Au.xlsx','Ti':'14Ti.xlsx'}
    for key,file in names.items():
        a=datafile(mat/file,w,unit='um');n,k=a[:,0],a[:,1]
        if np.min(n)<=0 or np.min(k)<-1e-8:raise ValueError(f'{file}: invalid n/k (passive material required)')
        if key in ('SiO2','PMMA') and np.max(np.abs(k))>1e-3:print(f'WARNING: {key} imaginary index ignored for lossless external/embedding model')
        d[key]=n+1j*(np.maximum(k,0) if key not in ('SiO2','PMMA') else 0)
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

def square_mask(width,nx,ny,device):
    def cover(n):
        x=np.arange(n)/n;lo=(1-width/P)/2;hi=(1+width/P)/2
        return np.clip(np.minimum(x+1/n,hi)-np.maximum(x,lo),0,None)*n
    return torch.as_tensor(cover(nx)[:,None]*cover(ny)[None,:],device=device,dtype=DTYPE)

def channel_power(s,lam,nsub,order,direction,port,device):
    # Verified p/s orthonormal channels for 2D: works at arbitrary diffraction azimuth
    nout=(1. if direction=='forward' else nsub) if port=='reflection' else (nsub if direction=='forward' else 1.)
    m=range(-order,order+1)
    orders=np.array([[i,j] for i in m for j in m if (i*i+j*j)*(lam/P)**2<nout*nout-1e-9],dtype=np.int64).reshape(-1,2)
    v=torch.zeros(len(orders),device=device,dtype=torch.float64)
    for pol in ('pp','sp'):
        a=s.S_parameters(orders=orders,direction=direction,port=port,polarization=pol,ref_order=[0,0],power_norm=True)
        v+=a.abs().square().reshape(-1).real
    v=v.detach().cpu().numpy()
    if not np.isfinite(v).all():raise FloatingPointError('Nonfinite scattering amplitude')
    return float(v.sum()),float(v[np.where((orders==0).all(axis=1))[0][0]])

def solve_one(d,j,phase,bias,device,order,mask,patch_h,with_vo2=True):
    w=float(d['w'][j]);lam=w+1e-5 if abs(w-P)<1e-7 else w
    nsub=float(d['SiO2'][j].real);pm=complex(d['PMMA'][j]**2)
    n_da=d['DASTn'][j]+d['beta'][j]*bias/H_DAST+1j*d['DASTk'][j]
    if n_da.real<=0:raise ValueError('DAST effective n becomes nonphysical; lower bias')
    s=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    s.add_output_layer(eps=nsub*nsub);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
    s.add_layer(thickness=H_PMMA_TOP,eps=pm)
    nv=d['VI'][j] if phase=='I' else d['VM'][j]
    s.add_layer(thickness=patch_h,eps=pm+mask*(complex(nv*nv)-pm) if with_vo2 else pm)
    s.add_layer(thickness=H_PMMA_BOTTOM+H_VO2-patch_h,eps=pm) # constant 50-nm PMMA-zone total in scans
    for material,h,n in [('Ti',H_TI,d['Ti'][j]),('ITO1',H_ITO1,d['ITO'][j]),('DAST',H_DAST,n_da),
                          ('ITO2',H_ITO2,d['ITO'][j]),('SiO2',H_SIO2,d['SiO2'][j]),('Au',H_AU,d['Au'][j])]:
        if h>0:s.add_layer(thickness=h,eps=complex(n*n))
    s.solve_global_smatrix();q=np.zeros((2,4)) # directions F/B; columns R,T,A,T0
    for di,direction in enumerate(('forward','backward')):
        q[di,0],_=channel_power(s,lam,nsub,order,direction,'reflection',device)
        q[di,1],q[di,3]=channel_power(s,lam,nsub,order,direction,'transmission',device)
        q[di,2]=1-q[di,0]-q[di,1]
    if not np.isfinite(q).all() or q.min()<-3e-4 or q.max()>1+3e-4:
        raise FloatingPointError(f'PASSIVITY failure at {w:g}nm {phase} {bias:+g}V: {q}; DO NOT clip absorption')
    return q

def spectrum(d,phase,bias,device,order,mask,patch_h,vo2=True,label=''):
    t=time.perf_counter();q=[]
    with torch.no_grad():
        for j in range(len(d['w'])):
            q.append(solve_one(d,j,phase,bias,device,order,mask,patch_h,vo2))
            if (j+1)%max(25,len(d['w'])//4)==0:
                if device.type=='cuda':torch.cuda.synchronize()
                print(f'  {label}: {j+1}/{len(d["w"])}; {time.perf_counter()-t:.1f}s',flush=True)
    if device.type=='cuda':torch.cuda.synchronize()
    q=np.array(q);print(f'{label}: {len(q)} wavelengths, {time.perf_counter()-t:.1f}s',flush=True)
    return q

def metrics(d,q):
    w=d['w'];v=(w>=400)&(w<=780);s=d['solar'];lum=d['lum']
    def avg(x,weight,which):return float(TRAPZ(x[which]*weight[which],w[which])/TRAPZ(weight[which],w[which]))
    af=avg(q[:,0,2],s,slice(None));ab=avg(q[:,1,2],s,slice(None))
    tf=avg(q[:,0,1],lum,v);tb=avg(q[:,1,1],lum,v)
    return {'DA':af-ab,'AsF':af,'AsB':ab,'TlumF':tf,'TlumB':tb,'Tlum':min(tf,tb),
            'TsolF':avg(q[:,0,1],s,slice(None)),'TsolB':avg(q[:,1,1],s,slice(None))}

def checks(d,device,order):
    # Bare-interface Fresnel reference, one lossless air film prevents empty-stack ambiguity
    worst=0.;w=d['w'];
    for wl in (405.,550.,900.,1550.,2400.):
        n=float(np.interp(wl,w,d['SiO2'].real));s=torcwa.rcwa(freq=1/wl,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
        s.add_output_layer(eps=n*n);s.set_incident_angle(inc_ang=0.,azi_ang=0.);s.add_layer(thickness=10.,eps=1.);s.solve_global_smatrix()
        re=((1-n)/(1+n))**2;te=4*n/(1+n)**2
        for side in ('forward','backward'):
            rr,_=channel_power(s,wl,n,order,side,'reflection',device);tt,_=channel_power(s,wl,n,order,side,'transmission',device)
            worst=max(worst,abs(rr-re),abs(tt-te),abs(rr+tt-1))
    print(f'Fresnel interface max difference={worst:.2e}')
    if worst>1e-6:raise RuntimeError('Fresnel validation failed; stop before scientific claims')
    # 2D lossless patterned dielectric energy test, including all propagating orders
    mask=square_mask(W_VO2,40,40,device);worst2=0.
    for wl in (405.,415.,500.,890.,1300.):
        n=float(np.interp(wl,w,d['SiO2'].real));s=torcwa.rcwa(freq=1/wl,order=[min(order,2),min(order,2)],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
        s.add_output_layer(eps=n*n);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
        s.add_layer(thickness=20.,eps=1.48**2+mask*(1.72**2-1.48**2));s.solve_global_smatrix()
        for side in ('forward','backward'):
            r,_=channel_power(s,wl,n,min(order,2),side,'reflection',device)
            t,_=channel_power(s,wl,n,min(order,2),side,'transmission',device)
            worst2=max(worst2,abs(r+t-1))
    print(f'2D lossless grating max |R+T-1|={worst2:.2e}')
    if worst2>1e-4:raise RuntimeError('Lossless 2D passivity failed; abort device simulation')

def sketch(ax,width,patch_h):
    z=0
    def layer(name,h,c,patch=False):
        nonlocal z
        if h<=0:return
        if patch:
            ax.add_patch(Rectangle(((P-width)/2,z),width,h,color=c,ec='k',lw=.4))
        else:ax.add_patch(Rectangle((0,z),P,h,color=c,ec='k',lw=.4))
        ax.text(P/2,z+h/2,name,ha='center',va='center',fontsize=8);z+=h
    layer('PMMA',H_PMMA_TOP,'#dceeea')
    ax.add_patch(Rectangle((0,z),P,patch_h,color='#dceeea',ec='k',lw=.4))
    layer('VO2 square in PMMA',patch_h,'#d9a176',True)
    layer('PMMA',H_PMMA_BOTTOM,'#dceeea')
    for name,h,c in [('Ti',H_TI,'#99999e'),('ITO1',H_ITO1,'#89c3d6'),('DAST',H_DAST,'#d6b1e7'),
                     ('ITO2',H_ITO2,'#89c3d6'),('SiO2',H_SIO2,'#cfdcea'),('Au',H_AU,'#e9c76d')]:layer(name,h,c)
    ax.text(P/2,z+10,'SiO2 substrate  ↑ B',fontsize=9,ha='center');ax.text(P/2,-12,'Air  F ↓',fontsize=9,ha='center')
    ax.set(xlim=(0,P),ylim=(z+25,-25),xlabel='x (nm)',ylabel='z (nm)',title='(a) Periodic VO2 square + planar Janus cavity')

def figure(d,spec,df,control):
    w=d['w'];fig,axs=plt.subplots(2,4,figsize=(20,9.5),layout='constrained');ax=axs.ravel();sketch(ax[0],W_VO2,H_VO2)
    styles={'I':'tab:blue','M':'tab:red'}
    for p in PHASES:
        for v in (VOLTAGES[0],VOLTAGES[-1]):
            q=spec[(p,v)];ax[1].plot(w,q[:,0,2]-q[:,1,2],color=styles[p],ls='--' if v<0 else '-',label=f'{p}, {v:+g} V')
    ax[1].axhline(0,color='.5',lw=.8);ax[1].set(title='(b) Spectral directional absorption',xlabel='Wavelength (nm)',ylabel=r'$A_F-A_B$');ax[1].legend(fontsize=8)
    for p in PHASES:
        q=spec[(p,0.)];ax[2].plot(w,q[:,0,2],color=styles[p],label=f'{p}, F');ax[2].plot(w,q[:,1,2],ls='--',color=styles[p],label=f'{p}, B')
        ax[3].plot(w,q[:,0,1],color=styles[p],label=f'{p}, F');ax[3].plot(w,q[:,1,1],ls='--',color=styles[p],label=f'{p}, B')
    ax[2].set(title='(c) Two-side absorptance at 0 V',xlabel='Wavelength (nm)',ylabel='A',ylim=(0,1));ax[2].legend(fontsize=8)
    ax[3].set(title='(d) Two-side transmittance at 0 V',xlabel='Wavelength (nm)',ylabel='T',ylim=(0,1));ax[3].legend(fontsize=8)
    for p in PHASES:
        z=df[df.Phase==p].sort_values('Voltage');ax[4].plot(z.Voltage,z.DA,'-o',color=styles[p],label=p)
        ax[5].plot(z.Voltage,z.Tlum,'-o',color=styles[p],label=p)
    ax[4].axhline(0,lw=.8,color='.6');ax[4].set(title='(e) Electrical tunability',xlabel='DAST voltage (V)',ylabel=r'$\Delta A_{sol}=A_F-A_B$');ax[4].legend()
    ax[5].axhline(.6,lw=.8,ls=':',color='.6',label='Tlum=0.6 target');ax[5].set(title='(f) D65 luminous transmittance',xlabel='DAST voltage (V)',ylabel='Min(Tlum F,B)',ylim=(0,1));ax[5].legend(fontsize=8)
    di=df[df.Phase=='I'].set_index('Voltage');dm=df[df.Phase=='M'].set_index('Voltage');vs=np.array(VOLTAGES)
    ax[6].plot(vs,[dm.loc[v,'DA']-di.loc[v,'DA'] for v in vs],'-o',color='darkred');ax[6].axhline(0,color='.6',lw=.8)
    ax[6].set(title='(g) VO2 phase response',xlabel='DAST voltage (V)',ylabel=r'$\Delta A_{sol,M}-\Delta A_{sol,I}$')
    texts=['I / 0 V','M / 0 V'];values=[di.loc[0.,'DA'],dm.loc[0.,'DA']]
    if control is not None:texts+=['No VO2 / 0 V'];values+=[control['DA']]
    ax[7].bar(texts,values,color=['tab:blue','tab:red','gray'][:len(values)]);ax[7].set(title='(h) Solar directional absorption, controls',ylabel='Signed solar-weighted difference');ax[7].tick_params(axis='x',labelrotation=12)
    for a in ax:a.grid(alpha=.15)
    fig.suptitle('Exploratory Ti / ITO / DAST / ITO / SiO2 / Au Janus cavity + 2D VO2 metasurface (NOT yet optimized)',fontsize=12)
    return fig

def main():
    if flag not in (0,1):raise ValueError('flag must be 0 or 1')
    if flag and not torch.cuda.is_available():raise RuntimeError('CUDA unavailable; set flag=0')
    if H_DAST<=0 or max(map(abs,VOLTAGES))*1000/H_DAST>E_LIMIT+1e-6:raise ValueError('DAST electric field exceeds E_LIMIT V/um')
    if not 0<W_VO2<P or STEP_NM<=0 or ORDER<1 or NX<4*ORDER+4 or NY<4*ORDER+4:raise ValueError('Check P, W_VO2, ORDER and NX/NY')
    torch.set_num_threads(1);device=torch.device('cuda' if flag else 'cpu');
    if flag:torch.backends.cuda.matmul.allow_tf32=False
    if QUICK:step,order,nx,ny=40.,2,56,56
    else:step,order,nx,ny=STEP_NM,ORDER,NX,NY
    d=load_data(step);mask=square_mask(W_VO2,nx,ny,device)
    print(f'\nDEVICE={device}, PyTorch={torch.__version__}, CUDA available={torch.cuda.is_available()}, complex128')
    print(f'2D square VO2/P= {W_VO2:g}/{P:g}nm, h={H_VO2:g}nm; 2D order={order}, modes={(2*order+1)**2}, Nx×Ny={nx}×{ny}, {len(d["w"])} wavelengths')
    print('Geometry air / PMMA(VO2 square) / Ti / ITO / DAST / ITO / SiO2 / Au / semi-infinite SiO2')
    print(f'DAST_SOURCE={DAST_SOURCE}, scalar EO; Emax={max(map(abs,VOLTAGES))*1000/H_DAST:.2f} V/um; ITO={ITO_FILE} UNVERIFIED')
    checks(d,device,order);t=time.perf_counter();spec={};rows=[]
    for phase in PHASES:
        for voltage in VOLTAGES:
            key=(phase,float(voltage));q=spectrum(d,phase,float(voltage),device,order,mask,H_VO2,label=f'{phase} {voltage:+g} V');spec[key]=q
            rows.append(dict(Phase=phase,Voltage=voltage,**metrics(d,q)))
    table=pd.DataFrame(rows);ref=None
    if RUN_CONTROL:
        qc=spectrum(d,'I',0.,device,order,mask,H_VO2,False,label='No VO2 @0V');ref=metrics(d,qc)
        print(f'No-VO2 0V: DA={ref["DA"]:+.5f}, Tlum={ref["Tlum"]:.4f}')
    if SCAN:
        print('\nOPTIONAL COARSE DESIGN SCREEN: not a converged optimization')
        dc=load_data(20.);nn=2;maskdim=56
        for width,h in SEARCH_CANDIDATES:
            m=square_mask(width,maskdim,maskdim,device);tmp={}
            for ph in PHASES:
                for v in (VOLTAGES[0],VOLTAGES[-1]):
                    z=spectrum(dc,ph,v,device,nn,m,h,label=f'SCREEN {width:g}/{h:g} {ph} {v:+g}V');tmp[(ph,v)]=metrics(dc,z)
            ea=abs(tmp[('I',VOLTAGES[0])]['DA']-tmp[('I',VOLTAGES[-1])]['DA'])
            em=abs(tmp[('M',VOLTAGES[0])]['DA']-tmp[('M',VOLTAGES[-1])]['DA'])
            th=max(abs(tmp[('M',v)]['DA']-tmp[('I',v)]['DA']) for v in (VOLTAGES[0],VOLTAGES[-1]))
            tl=min(x['Tlum'] for x in tmp.values())
            print(f'SCREEN width={width:g} h={h:g}: |DA(I,-V)|={abs(tmp[("I",VOLTAGES[0])]["DA"]):.4f}, EO_I={ea:.4f}, EO_M={em:.4f}, thermal_max={th:.4f}, min Tlum={tl:.4f}',flush=True)
    worst=max(np.max(np.abs(q[:,0,3]-q[:,1,3])) for q in spec.values())
    print('\nRESULT TABLE\n'+table.to_string(index=False,float_format=lambda x:f'{x:.6f}'))
    for ph in PHASES:
        z=table[table.Phase==ph].set_index('Voltage');da0=z.loc[VOLTAGES[0],'DA'];da1=z.loc[VOLTAGES[-1],'DA']
        print(f'{ph}: signed DA endpoints {da0:+.5f}, {da1:+.5f}; EO swing={abs(da1-da0):.5f}; directional sign reversal={bool(da0*da1<0)}')
    print(f'Matched zero-order T reciprocity max difference={worst:.3e}; cross-order total T may differ when multiple ports open.')
    print(f'RCWA+diagnostics runtime={time.perf_counter()-t:.1f}s; total incl startup={time.perf_counter()-START:.1f}s; NO exports')
    fig=figure(d,spec,table,ref)
    if SHOW_FIGURE:plt.show()
    return table,spec,fig

if __name__=='__main__':ResultTable,Spectra,Figure=main()
