
#!/usr/bin/env python3
"""2D TORCWA: PMMA-embedded periodic VO2/Au contacting squares above planar ITO/DAST/ITO/SiO2.
Fast preliminary design, not converged, not an experimental performance claim.
Place in D:/09_PCloud/14_DAST_VO2 beside optical_constants, Atmospheric data, VO2_temperature_data.
Run directly from VS Code with xhh31132 Python. No MATLAB; one 2x4 figure, no exports.
"""
from pathlib import Path
import time
T0=time.perf_counter()
import numpy as np, pandas as pd, torch, torcwa
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
# ===== EDIT THESE SETTINGS ONLY =====
flag=0; QUICK=0                 # 0 CPU, 1 CUDA; QUICK=1 functional check (~40 nm)
STEP_NM=20.; ORDER=2; NX=NY=40 # preliminary 2D design, NOT convergence settings
P=400.; WIDTH=180.              # 2D square array period and square side, nm
H_CAP=20.; H_VO2=12.; H_AU=12.        # two contacting aligned squares embedded in PMMA; NO SiO2 spacer
H_ITO_TOP=5.; H_DAST=95.; H_ITO_BOTTOM=5.  # continuous electro-optic stack
VOLTAGES=(-8.5,0.,8.5); PHASES=('I','M')
RUN_CONTROL=1; SHOW_FIGURE=1    # 1 adds no-VO2 & no-squares 0-V control spectra
ITO_FILE='4TIO.xlsx'            # Verify it truly describes ITO, not titanium oxide!
DAST_SOURCE='hybrid'            # 'hybrid' existing exploratory n/k, or 'xlsx'
EO_SIGN=1.; E_LIMIT=100.       # scalar EO benchmark, field limit in V/um
PROJECT=Path(__file__).resolve().parent; DTYPE=torch.complex128
TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
# ===== NUMERIC DATA AND MATERIAL MODELS =====
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
def solve_one(d,j,phase,voltage,device,order,mask,control='full'):
    w=float(d['w'][j]);lam=w+1e-5 if abs(w-P)<1e-7 else w
    nsub=float(d['SiO2'][j].real);pm=complex(d['PMMA'][j]**2)
    nd=d['DASTn'][j]+d['beta'][j]*voltage/H_DAST+1j*d['DASTk'][j]
    if nd.real<=0:raise ValueError('Nonphysical DAST n: adjust voltage/electro-optic model')
    s=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    s.add_output_layer(eps=nsub*nsub);s.set_incident_angle(inc_ang=0.,azi_ang=0.)
    s.add_layer(thickness=H_CAP,eps=pm)
    for key,h in [('VO2_'+phase,H_VO2),('Au',H_AU)]:
        n=d[key][j];active=(control!='no_squares' and not(control=='no_vo2' and key.startswith('VO2_')))
        ep=pm+mask*(complex(n*n)-pm) if active else pm
        s.add_layer(thickness=h,eps=ep)
    for n,h in [(d['ITO'][j],H_ITO_TOP),(nd,H_DAST),(d['ITO'][j],H_ITO_BOTTOM)]:
        if h>0:s.add_layer(thickness=h,eps=complex(n*n))
    s.solve_global_smatrix();q=np.zeros((2,4)) # forward/backward × (R,T,A,T0)
    for i,direction in enumerate(('forward','backward')):
        q[i,0],_=port_power(s,lam,nsub,order,direction,'reflection',device)
        q[i,1],q[i,3]=port_power(s,lam,nsub,order,direction,'transmission',device)
        q[i,2]=1-q[i,0]-q[i,1]
    if not np.isfinite(q).all() or q.min()<-3e-4 or q.max()>1+3e-4:raise FloatingPointError(f'Passive device failure at {w:g} nm, {phase}, {voltage:+g} V, {control}: {q}; do not clip!')
    return q
def run_spectrum(d,phase,v,device,order,mask,control='full'):
    t=time.perf_counter();result=[]
    with torch.no_grad():
        for j in range(len(d['w'])):
            result.append(solve_one(d,j,phase,v,device,order,mask,control))
            if (j+1)%max(10,len(d['w'])//4)==0:
                if device.type=='cuda':torch.cuda.synchronize()
                print(f' {phase} {v:+g} V {control}: {j+1}/{len(d["w"])}; {time.perf_counter()-t:.1f}s',flush=True)
    if device.type=='cuda':torch.cuda.synchronize()
    q=np.asarray(result);print(f' DONE {phase} {v:+g} V {control}: {time.perf_counter()-t:.1f}s',flush=True)
    return q
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
# ===== SINGLE 2×4 MULTIPANEL FIGURE =====
def draw(d,spec,table,controls):
    w=d['w'];fig,axs=plt.subplots(2,4,figsize=(20,9.5),layout='constrained');a=axs.ravel();colors={'I':'tab:blue','M':'tab:red'}
    hlist=[('PMMA cap',H_CAP,'#dcece8',False),('VO2',H_VO2,'#c99472',True),('Au',H_AU,'#efce75',True),('ITO',H_ITO_TOP,'#8dc8d7',False),('DAST',H_DAST,'#d6afe5',False),('ITO',H_ITO_BOTTOM,'#8dc8d7',False)]
    z=0.;left=(P-WIDTH)/2
    for label,h,color,patch in hlist:
        if patch:a[0].add_patch(Rectangle((0,z),P,h,color='#dcece8',ec='k',lw=.4))
        a[0].add_patch(Rectangle((left if patch else 0,z),WIDTH if patch else P,h,color=color,ec='k',lw=.5))
        a[0].text(P/2,z+h/2,label,fontsize=8,ha='center',va='center');z+=h
    a[0].text(P/2,-11,'F: air ↓',ha='center');a[0].text(P/2,z+12,'B: SiO2 substrate ↑',ha='center')
    a[0].set(xlim=(0,P),ylim=(z+24,-21),xlabel='x (nm)',ylabel='z (nm)',title='(a) VO2/Au contacting squares in PMMA')
    for ph in PHASES:
        for v in (VOLTAGES[0],VOLTAGES[-1]):
            q=spec[(ph,v)];a[1].plot(w,q[:,0,2]-q[:,1,2],color=colors[ph],ls='--' if v<0 else '-',label=f'{ph}, {v:+g} V')
        zq=spec[(ph,0.)]
        for side in range(2):
            a[2].plot(w,zq[:,side,2],color=colors[ph],ls='-' if side==0 else '--',label=f'{ph} {"F" if side==0 else "B"}')
            a[3].plot(w,zq[:,side,1],color=colors[ph],ls='-' if side==0 else '--',label=f'{ph} {"F" if side==0 else "B"}')
        df=table[table.Phase==ph].sort_values('Voltage');a[4].plot(df.Voltage,df.DA,'-o',label=ph,color=colors[ph]);a[5].plot(df.Voltage,df.Tlum,'-o',label=ph,color=colors[ph])
    a[1].axhline(0,color='.5',lw=.8);a[1].set(title='(b) Spectral Janus absorption contrast',xlabel='Wavelength (nm)',ylabel='A_F - A_B');a[1].legend(fontsize=8)
    a[2].set(title='(c) Front/back absorption at 0 V',xlabel='Wavelength (nm)',ylabel='Absorptance',ylim=(0,1));a[2].legend(fontsize=8)
    a[3].set(title='(d) Front/back total T at 0 V',xlabel='Wavelength (nm)',ylabel='Transmittance',ylim=(0,1));a[3].legend(fontsize=8)
    a[4].axhline(0,color='.5',lw=.8);a[4].set(title='(e) Electrical control',xlabel='DAST voltage (V)',ylabel='Signed DA_sol');a[4].legend()
    a[5].axhline(.6,color='.5',ls=':',lw=.8);a[5].set(title='(f) D65 luminous transmission',xlabel='DAST voltage (V)',ylabel='Min(TlumF,TlumB)',ylim=(0,1));a[5].legend()
    di=table[table.Phase=='I'].set_index('Voltage');dm=table[table.Phase=='M'].set_index('Voltage');vs=np.asarray(VOLTAGES)
    a[6].plot(vs,[dm.loc[v,'DA']-di.loc[v,'DA'] for v in vs],'-o',color='darkred');a[6].axhline(0,color='.5',lw=.8)
    a[6].set(title='(g) VO2 thermal phase effect',xlabel='DAST voltage (V)',ylabel='DA_sol(M)-DA_sol(I)')
    labs=['I/0V','M/0V'];values=[di.loc[0.,'DA'],dm.loc[0.,'DA']];barcolors=['tab:blue','tab:red']
    for k,ctrl in controls.items():labs.append(k.replace(' ','\n'));values.append(ctrl['DA']);barcolors.append('gray' if 'VO2' in k else '#c3a268')
    a[7].bar(labs,values,color=barcolors);a[7].tick_params(axis='x',labelrotation=13)
    a[7].set(title='(h) Solar Janus contrast, controls',ylabel='Signed DA_sol')
    for ax in a:ax.grid(alpha=.14)
    fig.suptitle('Exploratory PMMA / [VO2-Au square array] / ITO / DAST / ITO / SiO2 Janus',fontsize=13)
    return fig
def main():
    if flag not in (0,1) or ORDER<1 or P<=WIDTH or WIDTH<=0:raise ValueError('Invalid flag, Fourier order or geometry')
    if flag and not torch.cuda.is_available():raise RuntimeError('CUDA unavailable: set flag=0')
    if max(abs(v) for v in VOLTAGES)*1000/H_DAST>E_LIMIT:raise ValueError('DAST field exceeds E_LIMIT V/um')
    if 0. not in VOLTAGES or len(PHASES)!=2 or set(PHASES)!={'I','M'}:raise ValueError('This figure expects 0 V and VO2 phases I/M')
    order,nx,ny,step=(1,24,24,40.) if QUICK else (ORDER,NX,NY,STEP_NM)
    if min(nx,ny)<4*order+4 or min(H_CAP,H_VO2,H_AU,H_DAST,H_ITO_TOP,H_ITO_BOTTOM)<=0:raise ValueError('Too coarse spatial grid or nonpositive layer thickness')
    torch.set_num_threads(1);device=torch.device('cuda' if flag else 'cpu')
    if flag:torch.backends.cuda.matmul.allow_tf32=False
    d=load_data(step);mask=square_mask(WIDTH,nx,ny,device)
    print(f'\nDEVICE={device}, PyTorch={torch.__version__}, CUDA available={torch.cuda.is_available()}, complex128; Nx×Ny={nx}×{ny}; N={order}; modes={(2*order+1)**2}; wavelengths={len(d["w"])}')
    print(f'PERIOD={P:g}, square={WIDTH:g}, cover={(WIDTH/P)**2:.2%}; VO2/Au contacting square heights={H_VO2:g}/{H_AU:g} nm (NO internal SiO2 spacer)')
    print(f'ITO/DAST/ITO={H_ITO_TOP:g}/{H_DAST:g}/{H_ITO_BOTTOM:g} nm; ITO source={ITO_FILE} UNVERIFIED; DAST={DAST_SOURCE} scalar EO; maximum field={max(abs(v) for v in VOLTAGES)*1000/H_DAST:.2f} V/um')
    validate(d,device,order,mask);begin=time.perf_counter();spec={};rows=[]
    for ph in PHASES:
        for v in VOLTAGES:
            q=run_spectrum(d,ph,v,device,order,mask);spec[(ph,float(v))]=q;rows.append(dict(Phase=ph,Voltage=v,**metrics(d,q)))
    table=pd.DataFrame(rows);controls={}
    if RUN_CONTROL:
        for label,mode in [('No VO2','no_vo2'),('No squares','no_squares')]:
            q=run_spectrum(d,'I',0.,device,order,mask,mode);controls[label]=metrics(d,q)
    print('\nRESULT TABLE\n'+table.to_string(index=False,float_format=lambda x:f'{x:.6f}'))
    for ph in PHASES:
        row=table[table.Phase==ph].set_index('Voltage');da=np.array([row.loc[v,'DA'] for v in VOLTAGES])
        print(f'{ph}: signed DA range=[{da.min():+.5f},{da.max():+.5f}], EO sampled swing={np.ptp(da):.5f}, endpoints opposite sign={bool(da[0]*da[-1]<0)}')
    for label,r in controls.items():print(f'{label} @0 V: DA={r["DA"]:+.5f}; Tlum={r["Tlum"]:.4f}')
    print('Thermal phase DA contrast at voltages:',{v:round(table[(table.Phase=='M')&(table.Voltage==v)].DA.iloc[0]-table[(table.Phase=='I')&(table.Voltage==v)].DA.iloc[0],6) for v in VOLTAGES})
    print('Zero-order matched transmission reciprocity max=',max(np.max(abs(q[:,0,3]-q[:,1,3])) for q in spec.values()))
    print(f'RCWA simulation={time.perf_counter()-begin:.1f}s; ALL startup + simulation={time.perf_counter()-T0:.1f}s; one figure only, no exports')
    fig=draw(d,spec,table,controls)
    if SHOW_FIGURE:plt.show()
    return table,spec,fig
if __name__=='__main__':ResultTable,Spectra,Figure=main()
