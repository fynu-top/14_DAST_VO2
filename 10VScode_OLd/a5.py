#!/usr/bin/env python3
"""Fast Au/SiO2 thickness scan: true 2D aligned VO2/SiO2/Au square stack in PMMA, above continuous ITO/DAST/ITO.
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
QUICK=0                             # 1: fewer wavelengths (40 nm); retains N=2
STEP_NM=40.; ORDER=2; NX=NY=32        # fast screening only; not convergence-certified
P=400.; WIDTH=180.                   # common aligned square width; area fill=(WIDTH/P)^2
H_VO2=12.; H_SIO2_PATCH=60.; H_AU=12.; H_PMMA_TOP=15.
H_ITO1=5.; H_DAST=95.; H_ITO2=5.      # all three layers continuous, no bottom Au sheet
VOLTAGES=(-8.5,0.,8.5); PHASES=('I','M'); E_LIMIT=100.; MARGIN=.005
EO_SIGN=1.; DAST_SOURCE='hybrid'; ITO_FILE='4TIO.xlsx'
LOSSLESS_PMMA=True                   # same transparent-host assumption as 10.py; set False to retain k
RUN_CONTROLS=0; RUN_OLD10=0           # no-VO2 + full-coverage & equal-volume films + old10 at 0V
CHECK_FOURIER=0; CHECK_WAVELENGTH=0  # Fourier/pixel probes; optional extra full spectrum at STEP/2
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


# Small discrete screen; Au=0 is a no-metal control, not an Au-containing design.
AU_LIST=[0.,4.,8.,12.]; SIO2_LIST=[20.,60.,100.,140.]
TLUM_FLOOR=.45; TARGET_PHASE='M'; DETAIL_POINTS=9
BASE_AU=12.; BASE_SIO2=60.

def main():
    global H_AU,H_SIO2_PATCH
    t0=time.perf_counter();torch.set_num_threads(1)
    if flag and not torch.cuda.is_available():raise RuntimeError('CUDA unavailable; use flag=0.')
    device=torch.device('cuda' if flag else 'cpu');torch.backends.cuda.matmul.allow_tf32=False
    if TARGET_PHASE not in ('I','M'):raise ValueError('TARGET_PHASE must be I or M.')
    if min(AU_LIST)<0 or min(SIO2_LIST)<0 or len(set(AU_LIST))!=len(AU_LIST) or len(set(SIO2_LIST))!=len(SIO2_LIST):raise ValueError('Use unique nonnegative thicknesses.')
    vmax=min(8.5,10.,E_LIMIT*H_DAST/1000);d=load_data(STEP_NM);mask=mask_square(WIDTH,NX,NY,device)
    print(f'FAST SCREEN: {len(AU_LIST)*len(SIO2_LIST)} candidates; N={ORDER}, pixels={NX}x{NY}, step={STEP_NM:g} nm, V=+/-{vmax:g} V; target phase={TARGET_PHASE}')
    checks(d,device,ORDER);cache={};rows=[]
    def calc(au,sp,phase,v):
        global H_AU,H_SIO2_PATCH
        key=(float(au),float(sp),phase,float(v));H_AU=float(au);H_SIO2_PATCH=float(sp)
        if key not in cache:
            with torch.no_grad():q=np.asarray([solve_one(d,j,phase,float(v),device,ORDER,mask) for j in range(len(d['w']))])
            cache[key]=(metrics(d,q),q)
        return cache[key]
    for sp in SIO2_LIST:
        for au in AU_LIST:
            row=dict(ID=len(rows)+1,Au_nm=au,SiO2_nm=sp,NoAuControl=au==0);ts=[]
            for phase in ('I','M'):
                a,_=calc(au,sp,phase,-vmax);b,_=calc(au,sp,phase,vmax);ts.extend([a['Tlum'],b['Tlum']])
                plus=min(a['DA'],-b['DA']);minus=min(-a['DA'],b['DA']);margin=max(plus,minus)
                row.update({f'DAminus_{phase}':a['DA'],f'DAplus_{phase}':b['DA'],f'Margin_{phase}':margin,f'Reversal_{phase}':margin>0,f'Orientation_{phase}':'+ to -' if plus>=minus else '- to +',f'ElectricalSpan_{phase}':abs(b['DA']-a['DA'])})
            # Thermal phase reversal at SAME tested bias, either orientation.
            row['PhaseMargin_endpoints']=max(max(min(row[f'DA{x}_I'],-row[f'DA{x}_M']),min(-row[f'DA{x}_I'],row[f'DA{x}_M'])) for x in ('minus','plus'))
            row['Tlum_endpoint_min']=min(ts);row['Qualified']=row[f'Margin_{TARGET_PHASE}']>=MARGIN and min(ts)>=TLUM_FLOOR
            rows.append(row);print(f"{len(rows)}/{len(AU_LIST)*len(SIO2_LIST)} Au={au:g}, spacer={sp:g}: margin {TARGET_PHASE}={row[f'Margin_{TARGET_PHASE}']:+.5f}, T_min={min(ts):.3f}; elapsed={time.perf_counter()-t0:.1f}s",flush=True)
    table=pd.DataFrame(rows);physical=table[~table.NoAuControl];pool=physical if len(physical) else table
    good=pool[pool.Qualified]
    if len(good):selected=good.sort_values(['Tlum_endpoint_min',f'Margin_{TARGET_PHASE}'],ascending=False).iloc[0];reason='largest endpoint transparency among qualified Au-containing sampled candidates'
    else:
        bright=pool[pool.Tlum_endpoint_min>=TLUM_FLOOR];pool=bright if len(bright) else pool
        selected=pool.sort_values([f'Margin_{TARGET_PHASE}','Tlum_endpoint_min'],ascending=False).iloc[0];reason='NO qualified Au-containing candidate: diagnostic best margin (not a successful reversal design)'
    au=float(selected.Au_nm);sp=float(selected.SiO2_nm);detail=[];spectra={}
    for phase in ('I','M'):
        for v in np.unique(np.r_[np.linspace(-vmax,vmax,DETAIL_POINTS),0.]):
            met,q=calc(au,sp,phase,v);detail.append(dict(Phase=phase,Voltage=v,**met));spectra[(phase,float(v))]=q
    detail=pd.DataFrame(detail);base=[]
    for phase in ('I','M'):
        for v in (-vmax,0.,vmax):
            met,_=calc(BASE_AU,BASE_SIO2,phase,v);base.append(dict(Design='Baseline',Phase=phase,Voltage=v,Au_nm=BASE_AU,SiO2_nm=BASE_SIO2,**met))
            met,_=calc(au,sp,phase,v);base.append(dict(Design='Selected',Phase=phase,Voltage=v,Au_nm=au,SiO2_nm=sp,**met))
    comparison=pd.DataFrame(base);print('\nSCAN RESULTS\n'+table.to_string(index=False));print('\nSELECTION: '+reason);print(selected.to_string());print('\nBASELINE / SELECTED\n'+comparison.to_string(index=False))
    H_AU=au;H_SIO2_PATCH=sp
    fig,axs=plt.subplots(2,3,figsize=(16,9),constrained_layout=True);axes=axs.ravel()
    for ax,key,title in zip(axes[:3],['Margin_I','Margin_M','Tlum_endpoint_min'],['I phase: endpoint reversal margin','M phase: endpoint reversal margin','Minimum endpoint luminous T (both phases)']):
        z=table.pivot(index='SiO2_nm',columns='Au_nm',values=key).reindex(index=SIO2_LIST,columns=AU_LIST).to_numpy()
        im=ax.imshow(z,origin='lower',aspect='auto',cmap='viridis');fig.colorbar(im,ax=ax)
        ax.set(xticks=range(len(AU_LIST)),xticklabels=AU_LIST,yticks=range(len(SIO2_LIST)),yticklabels=SIO2_LIST,xlabel='Au thickness (nm); 0 = control',ylabel='SiO2 patch thickness (nm)',title=title)
        ax.plot(AU_LIST.index(au),SIO2_LIST.index(sp),'r*',ms=14)
        for iy in range(len(SIO2_LIST)):
            for ix in range(len(AU_LIST)):ax.text(ix,iy,f'{z[iy,ix]:.3f}',ha='center',va='center',color='white',fontsize=8)
    for phase,c in [('I','tab:blue'),('M','tab:red')]:
        x=detail[detail.Phase==phase];axes[3].plot(x.Voltage,x.DA,'-o',color=c,label=phase);axes[4].plot(x.Voltage,x.Tlum,'-o',color=c,label=phase)
        for v,ls in [(-vmax,'--'),(vmax,'-')]:
            q=spectra[(phase,v)].mean(axis=2);axes[5].plot(d['w'],q[:,0,2]-q[:,1,2],ls,color=c,label=f'{phase} {v:+g}V')
    axes[3].axhline(0,color='.5');axes[3].axhline(MARGIN,color='.6',ls=':');axes[3].axhline(-MARGIN,color='.6',ls=':');axes[4].axhline(TLUM_FLOOR,color='.5',ls=':')
    axes[3].set(title='Selected: sampled voltage response',xlabel='Voltage (V)',ylabel='Solar-weighted A_F - A_B');axes[4].set(title='Selected: min(F,B) luminous T',xlabel='Voltage (V)',ylabel='D65 weighted T');axes[5].set(title='Selected: directional spectra',xlabel='Wavelength (nm)',ylabel='A_F - A_B')
    for ax in axes[3:]:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle(f'FAST exploratory scan: Au={au:g} nm, SiO2 patch={sp:g} nm; no convergence claim')
    print(f'\nElapsed={time.perf_counter()-t0:.1f}s. ONE figure; no exports. Au=0 control excluded from selection when Au-containing candidates exist.')
    print('Margins permit either reversal orientation. Target phase controls selection; transparency floor applies to both phases at both voltage endpoints.')
    print('Endpoint reversal is not proof of a minimum switching voltage. Interior samples are not a continuous-voltage proof. No convergence, thermal feedback or thickness-dependent material model is included.')
    if SHOW_FIGURE:plt.show()
    return dict(ScanTable=table,Selected=selected,DetailTable=detail,ComparisonTable=comparison,Spectra=spectra,Figure=fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path);p.add_argument('--no-show',action='store_true');args=p.parse_args()
    if args.project is not None:PROJECT=args.project.resolve()
    if args.no_show:SHOW_FIGURE=0
    Results=main();globals().update(Results)
