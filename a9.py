#!/usr/bin/env python3
"""Constant-field DAST/ITO thickness NO-Au Janus RCWA search; pure Python, original material folders, no exports.
Air/PMMA cap/VO2 squares/SiO2 (aligned squares OR continuous)/ITO/DAST/ITO/SiO2 half-space.
Same-height spacer comparison, NOT equal-volume comparison. Trial effective scalar materials.
Python 3.9+; numpy scipy pandas openpyxl matplotlib torch torcwa (API 0.1.4.2).
Run beside optical_constants, Atmospheric data, VO2_temperature_data. No other .py required.
"""
from pathlib import Path
import argparse,time
import numpy as np,pandas as pd,torch,torcwa
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
# ===== All lengths nm; only change these for a different exploratory run =====
flag=0; STEP_NM=40.; ORDER=2; NX=NY=32; SHOW_FIGURE=True
P=400.; WIDTH=220.; H_VO2=8.; H_PMMA_TOP=15.
H_ITO1=5.; H_DAST=95.; H_ITO2=5.; VEND=8.5; E_LIMIT=100.
SPACERS=[55.,60.,65.]; WIDTHS=[200.,220.,240.]; HEIGHTS=[6.,8.,10.]
LOCAL_ROUNDS=2; DETAIL_POINTS=7; TLUM_FLOOR=.50; MARGIN=.005; FONT_SIZE=11
EO_SIGN=1.; DAST_SOURCE='hybrid'; ITO_FILE='4TIO.xlsx'; LOSSLESS_PMMA=True
PROJECT=Path(__file__).resolve().parent; DTYPE=torch.complex128
TRAPZ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
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
                     'PMMA':'9PMMA.xlsx'}
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


def solve_one(d,j,phase,bias,device,order,mask,geometry):
    mode,sp,width,height,hd,hi1,hi2=geometry;w=float(d['w'][j]);ns=float(d['SiO2'][j].real);lam=safe_lambda(w,ns,order)
    pm=complex(d['PMMA'][j]**2);nv=d['VI'][j] if phase=='I' else d['VM'][j]
    nd=d['DASTn'][j]+d['beta'][j]*bias/hd+1j*d['DASTk'][j]
    if nd.real<=0:raise ValueError('Nonpositive DAST real index.')
    sim=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    sim.add_output_layer(eps=ns*ns);sim.set_incident_angle(inc_ang=0.,azi_ang=0.)
    def layer(h,eps):
        if h>0:sim.add_layer(thickness=h,eps=eps)
    layer(H_PMMA_TOP,pm);layer(height,pm+mask*(complex(nv*nv)-pm))
    eps=complex(d['SiO2'][j]**2);layer(sp,pm+mask*(eps-pm) if mode=='patch' else eps)
    for h,n in [(hi1,d['ITO'][j]),(hd,nd),(hi2,d['ITO'][j])]:layer(h,complex(n*n))
    sim.solve_global_smatrix();return powers(sim,lam,ns,order,device)

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




# Constant-field comparison preserves the original 95 nm / 8.5 V reference.
DAST_VALUES=[75.,95.,115.]; ITO_VALUES=[5.,8.,12.]
E_SCAN=8500./95.; V_LIMIT=12.; DETAIL_POINTS=7

def margin(a,b):return max(min(a,-b),min(-a,b))
def endpoint_voltage(g):return E_SCAN*g[4]/1000.
def objectives(ep,v):
    out={}
    for ph in ('I','M'):
        a,b=ep[ph,-v],ep[ph,v];span=abs(b['DA']-a['DA'])
        out.update({f'DAminus_{ph}':a['DA'],f'DAplus_{ph}':b['DA'],f'Margin_{ph}':margin(a['DA'],b['DA']),f'SpanDA_{ph}':span,f'SecantDA_perV_{ph}':span/(2*v),f'SpanTsol_{ph}':abs(b['TsolF']-a['TsolF'])})
    out.update(Tmin=min(x['Tlum'] for x in ep.values()),WorstMargin=min(out['Margin_I'],out['Margin_M']),PhaseTsolAbs=max(abs(ep['I',x]['TsolF']-ep['M',x]['TsolF']) for x in (-v,v)),PhaseAbsMargin=max(margin(ep['I',x]['DA'],ep['M',x]['DA']) for x in (-v,v)))
    out['Transparent']=out['Tmin']>=TLUM_FLOOR;out['Qualified']=out['Transparent'] and out['WorstMargin']>=MARGIN
    return out

def choose(t):
    q=t[t.Transparent];cols=['WorstMargin','Tmin','PhaseTsolAbs'] if len(q) else ['Tmin','WorstMargin','PhaseTsolAbs']
    return (q if len(q) else t).sort_values(cols,ascending=False,kind='stable').iloc[0]

def main():
    started=time.perf_counter();torch.set_num_threads(1)
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise ValueError('flag=1 requires CUDA; otherwise use flag=0.')
    if not 0<E_SCAN<=E_LIMIT or min(DAST_VALUES)<=0 or min(ITO_VALUES)<=0:raise ValueError('Invalid field or thickness.')
    if max(DAST_VALUES)*E_SCAN/1000>V_LIMIT:raise ValueError('Increase V_LIMIT or reduce DAST_VALUES/E_SCAN; no silent field clipping.')
    if ORDER<1 or min(NX,NY)<=4*ORDER or STEP_NM<=0 or DETAIL_POINTS<3:raise ValueError('Invalid numerical settings.')
    device=torch.device('cuda' if flag else 'cpu');torch.backends.cuda.matmul.allow_tf32=False
    d=load_data(STEP_NM);checks(d,device,ORDER);cache={};masks={};records={}
    baseline={m:(m,60.,200. if m=='patch' else 220.,8.,95.,5.,5.) for m in ('patch','film')}
    print(f'Constant field +/-{E_SCAN:.6g} V/um; N={ORDER}, grid={NX}x{NY}, wavelength step={STEP_NM:g} nm.')
    print('ITO1=ITO2 in this search. Ranking: transparency, weaker-phase margin, transparency tie-break, phase T contrast.')
    def state(g,ph,v):
        key=(g,ph,round(float(v),12))
        if key not in cache:
            if g[2] not in masks:masks[g[2]]=mask_square(g[2],NX,NY,device)
            with torch.no_grad():q=np.asarray([solve_one(d,j,ph,float(v),device,ORDER,masks[g[2]],g) for j in range(len(d['w']))])
            cache[key]=metrics(d,q)
        return cache[key]
    def evaluate(g):
        if g not in records:
            v=endpoint_voltage(g);ep={(ph,x):state(g,ph,x) for ph in ('I','M') for x in (-v,v)}
            row=dict(ID=len(records)+1,Mode=g[0],Width_nm=g[2],DAST_nm=g[4],ITO1_nm=g[5],ITO2_nm=g[6],Vend=v,E_Vum=E_SCAN,**objectives(ep,v));records[g]=row
            print(f"ID {row['ID']:02d} {g[0]} DAST/ITO={g[4]:g}/{g[5]:g} nm, V=+/-{v:.3f}: Tmin={row['Tmin']:.5f}, margin I/M={row['Margin_I']:+.6f}/{row['Margin_M']:+.6f}; {time.perf_counter()-started:.1f}s",flush=True)
        return records[g]
    for mode,b in baseline.items():
        evaluate(b)
        for hd in DAST_VALUES:
            for hi in ITO_VALUES:evaluate((*b[:4],float(hd),float(hi),float(hi)))
    table=pd.DataFrame(records.values());winners={};comparison=[];detail=[];summary=[]
    for mode,b in baseline.items():
        best=choose(table[table.Mode==mode]);g=next(g for g,r in records.items() if r['ID']==best.ID);winners[mode]=g
        for role,h in [('Baseline',b),('Best',g)]:
            comparison.append(dict(Role=role,**records[h]));v=endpoint_voltage(h)
            for ph in ('I','M'):
                for x in np.unique(np.r_[np.linspace(-v,v,DETAIL_POINTS),0.]):detail.append(dict(Mode=mode,Role=role,Phase=ph,Voltage=x,Field=x*1000/h[4],**state(h,ph,float(x))))
    detail=pd.DataFrame(detail);comparison=pd.DataFrame(comparison)
    for mode,g in winners.items():
        t=detail[(detail.Mode==mode)&(detail.Role=='Best')];i=t[t.Phase=='I'].set_index('Voltage');m=t[t.Phase=='M'].set_index('Voltage')
        z=np.maximum(np.minimum(i.DA,-m.DA),np.minimum(-i.DA,m.DA));j=int(np.argmax(z.to_numpy()));r=records[g]
        summary.append(dict(Mode=mode,DAST_nm=g[4],ITO_nm=g[5],Vend=endpoint_voltage(g),MinTlum_sample=t.Tlum.min(),WorstMargin=r['WorstMargin'],Qualified_sample=bool(t.Tlum.min()>=TLUM_FLOOR and r['WorstMargin']>=MARGIN),PhaseMargin_sample=z.iloc[j],PhaseBias_sample=z.index[j],PhaseTsol_sample=max(abs(i.TsolF-m.TsolF)),SpanDA_I_sample=i.DA.max()-i.DA.min(),SpanDA_M_sample=m.DA.max()-m.DA.min(),SpanTsol_I_sample=i.TsolF.max()-i.TsolF.min(),SpanTsol_M_sample=m.TsolF.max()-m.TsolF.min()))
    summary=pd.DataFrame(summary)
    for title,t in [('ALL ENDPOINT TRIALS',table),('BASELINE / BEST',comparison),('VOLTAGE SAMPLES',detail),('SAMPLED PERFORMANCE',summary)]:print('\n'+title+'\n'+t.to_string(index=False,float_format=lambda x:f'{x:.6g}'))
    best=choose(table);print('\nPREFERRED ENDPOINT CANDIDATE\n'+best.to_string())
    if not best.Qualified:print('Preferred candidate is diagnostic: endpoint targets not met.')
    if not summary.loc[summary.Mode==best.Mode,'Qualified_sample'].iloc[0]:print('Preferred candidate fails sampled qualification; inspect SampleTable before selecting.')
    # Equal-voltage control: compare all four baseline/best structures over the shared safe voltage range.
    vc=min(endpoint_voltage(g) for g in [*baseline.values(),*winners.values()]);control=[]
    for mode in baseline:
        for role,g in [('Baseline',baseline[mode]),('Best',winners[mode])]:
            ep={(ph,v):state(g,ph,v) for ph in ('I','M') for v in (-vc,vc)}
            control.append(dict(Mode=mode,Role=role,Vend=vc,E_Vum=vc*1000/g[4],DAST_nm=g[4],ITO_nm=g[5],**objectives(ep,vc)))
    control=pd.DataFrame(control);print('\nCOMMON-VOLTAGE CONTROL (same voltage, different field)\n'+control.to_string(index=False,float_format=lambda x:f'{x:.6g}'))
    plt.rcParams.update({'font.family':'Times New Roman','font.size':FONT_SIZE,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in'})
    fig,aa=plt.subplots(2,4,figsize=(19,9),constrained_layout=True);ax=aa.ravel();colors={'patch':'tab:orange','film':'tab:blue'}
    for mode in baseline:
        t=table[table.Mode==mode];c=colors[mode];ax[0].scatter(t.Tmin,t.WorstMargin,c=c,label=mode)
        for r in t.itertuples():ax[0].annotate(str(r.ID),(r.Tmin,r.WorstMargin),fontsize=7)
        for role,alpha in [('Baseline',.35),('Best',1.)]:
            for ph,ls in [('I','-'),('M','--')]:
                t=detail[(detail.Mode==mode)&(detail.Role==role)&(detail.Phase==ph)];label=f'{mode} {role} {ph}'
                for k,col in [(2,'DA'),(3,'Tlum'),(4,'TsolF')]:ax[k].plot(t.Voltage,t[col],ls+'o',c=c,alpha=alpha,label=label,ms=3)
        t=detail[(detail.Mode==mode)&(detail.Role=='Best')];i=t[t.Phase=='I'].set_index('Voltage');m=t[t.Phase=='M'].set_index('Voltage')
        z=np.maximum(np.minimum(i.DA,-m.DA),np.minimum(-i.DA,m.DA));ax[5].plot(z.index,z,'-o',c=c,label=mode);ax[6].plot(i.index,i.TsolF-m.TsolF,'-o',c=c,label=mode)
    for k,t in [(1,comparison),(7,control)]:
        x=np.arange(len(t));ax[k].bar(x-.18,t.Margin_I,.36,label='I');ax[k].bar(x+.18,t.Margin_M,.36,label='M');ax[k].set_xticks(x);ax[k].set_xticklabels([f'{r.Mode}\n{r.Role}' for r in t.itertuples()]);ax[k].axhline(MARGIN,c='.5',ls=':')
    settings=[('Thickness trials','Min endpoint luminous T','Weaker-phase margin'),('Equal-field margin','','Reversal margin'),('Electrical directional response','Voltage (V)','Solar A_F - A_B'),('Luminous transmittance','Voltage (V)','Min directional luminous T'),('Forward solar transmittance','Voltage (V)','Solar T_F'),('Phase directional switching','Voltage (V)','Phase reversal margin'),('Phase solar contrast','Voltage (V)','Solar T_F,I - T_F,M'),('Common-voltage margin','','Reversal margin')]
    for k,(title,x,y) in enumerate(settings):
        ax[k].set(xlabel=x,ylabel=y);ax[k].set_title(f'({chr(97+k)}) {title}',loc='left',fontsize=FONT_SIZE+2);ax[k].grid(alpha=.15);ax[k].legend(fontsize=7)
    ax[0].axvline(TLUM_FLOOR,c='.5',ls=':');ax[0].axhline(MARGIN,c='.5',ls=':');ax[3].axhline(TLUM_FLOOR,c='.5',ls=':')
    for k in (2,5,6):ax[k].axhline(0,c='.5',lw=.8)
    print(f'\nElapsed={time.perf_counter()-started:.1f}s; {len(records)} structures; {len(cache)} spectra. ONE figure; no exports.')
    print('Screening only: fixed effective nk and EO model, no thickness-dependent electrode continuity/resistance or thermal dynamics. F=air, B=SiO2 half-space; total directional T can differ due to diffraction.')
    print('Equal-field improvement can require higher voltage. Common-voltage control isolates that comparison. SecantDA_perV is an endpoint average, not a local derivative. Sampled checks are not continuous-voltage proofs.')
    if SHOW_FIGURE:plt.show()
    return dict(AllTable=table,ComparisonTable=comparison,DetailTable=detail,SampleTable=summary,CommonVoltageTable=control,Selected=best,Figure=fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path);p.add_argument('--no-show',action='store_true');a=p.parse_args()
    if a.project is not None:PROJECT=a.project.resolve()
    if a.no_show:SHOW_FIGURE=False
    Results=main();globals().update(Results)
