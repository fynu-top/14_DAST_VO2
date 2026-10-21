#!/usr/bin/env python3
"""Two-round bounded local NO-Au Janus RCWA search; pure Python, original material folders, no exports.
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
    mode,sp,width,height=geometry;w=float(d['w'][j]);ns=float(d['SiO2'][j].real);lam=safe_lambda(w,ns,order)
    pm=complex(d['PMMA'][j]**2);nv=d['VI'][j] if phase=='I' else d['VM'][j]
    nd=d['DASTn'][j]+d['beta'][j]*bias/H_DAST+1j*d['DASTk'][j]
    if nd.real<=0:raise ValueError('Nonpositive DAST real index.')
    sim=torcwa.rcwa(freq=1/lam,order=[order,order],L=[P,P],dtype=DTYPE,device=device,stable_eig_grad=False)
    sim.add_output_layer(eps=ns*ns);sim.set_incident_angle(inc_ang=0.,azi_ang=0.)
    def layer(h,eps):
        if h>0:sim.add_layer(thickness=h,eps=eps)
    layer(H_PMMA_TOP,pm);layer(height,pm+mask*(complex(nv*nv)-pm))
    eps=complex(d['SiO2'][j]**2);layer(sp,pm+mask*(eps-pm) if mode=='patch' else eps)
    for h,n in [(H_ITO1,d['ITO'][j]),(H_DAST,nd),(H_ITO2,d['ITO'][j])]:layer(h,complex(n*n))
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




def objectives(endpoints):
    out={};ts=[]
    for phase in ('I','M'):
        a,b=endpoints[(phase,-VEND)],endpoints[(phase,VEND)];pm=min(a['DA'],-b['DA']);mp=min(-a['DA'],b['DA'])
        out.update({f'DAminus_{phase}':a['DA'],f'DAplus_{phase}':b['DA'],f'Margin_{phase}':max(pm,mp),f'Orientation_{phase}':'+ to -' if pm>=mp else '- to +',f'VoltageSpanDA_{phase}':abs(b['DA']-a['DA'])})
        ts.extend([a['Tlum'],b['Tlum']])
    contrast=[endpoints[('I',v)]['TsolF']-endpoints[('M',v)]['TsolF'] for v in (-VEND,VEND)];ix=int(np.argmax(np.abs(contrast)))
    out.update(Tmin=min(ts),WorstMargin=min(out['Margin_I'],out['Margin_M']),PhaseTsolAbs=abs(contrast[ix]),PhaseTsolSigned=contrast[ix],PhaseTsolBias=(-VEND,VEND)[ix])
    out['PhaseAbsMargin']=max(max(min(endpoints[('I',v)]['DA'],-endpoints[('M',v)]['DA']),min(-endpoints[('I',v)]['DA'],endpoints[('M',v)]['DA'])) for v in (-VEND,VEND))
    out['Transparent']=out['Tmin']>=TLUM_FLOOR;out['BothReverse']=out['WorstMargin']>0
    out['Qualified']=out['Transparent'] and out['WorstMargin']>=MARGIN
    return out

def pareto_mask(table):
    x=table[['WorstMargin','PhaseTsolAbs']].to_numpy();ok=table.Transparent.to_numpy(bool);answer=np.zeros(len(x),bool)
    for i in np.flatnonzero(ok):answer[i]=not np.any(ok & np.all(x>=x[i],axis=1) & np.any(x>x[i],axis=1))
    return answer

def choose(table):
    good=table[table.Transparent]
    if len(good):return good.sort_values(['WorstMargin','Tmin','PhaseTsolAbs'],ascending=False).iloc[0]
    return table.sort_values(['Tmin','WorstMargin','PhaseTsolAbs'],ascending=False).iloc[0]

def neighbors(g):
    out=[g]
    for index,values in [(1,SPACERS),(2,WIDTHS),(3,HEIGHTS)]:
        for value in values:
            h=list(g);h[index]=float(value);out.append(tuple(h))
    return list(dict.fromkeys(out))

def main():
    started=time.perf_counter();torch.set_num_threads(1)
    if flag not in (0,1) or (flag and not torch.cuda.is_available()):raise ValueError('Use CPU flag=0 or available CUDA flag=1.')
    if not 0<VEND<=10 or H_DAST<=0 or VEND*1000/H_DAST>E_LIMIT:raise ValueError('Invalid voltage/field.')
    if STEP_NM<=0 or ORDER<1 or min(NX,NY)<=4*ORDER or DETAIL_POINTS<3 or LOCAL_ROUNDS<1:raise ValueError('Invalid numerical settings.')
    if min(SPACERS)<0 or min(HEIGHTS)<=0 or not all(0<w<P for w in WIDTHS):raise ValueError('Invalid search bounds.')
    device=torch.device('cuda' if flag else 'cpu');torch.backends.cuda.matmul.allow_tf32=False;d=load_data(STEP_NM);checks(d,device,ORDER)
    print(f'LOCAL SEARCH: {LOCAL_ROUNDS} rounds, N={ORDER}, pixels={NX}x{NY}, step={STEP_NM:g} nm. No Au; phase I/M; V=+/-{VEND:g} V.')
    print('Ranking: endpoint transparency >=50% (editable), then largest weaker-phase margin; ties: transparency, then phase T contrast.')
    cache={};records={};masks={};history=[];baseline={mode:(mode,60.,220.,8.) for mode in ('film','patch')}
    def state(g,ph,v):
        key=(g,ph,float(v))
        if key not in cache:
            if g[2] not in masks:masks[g[2]]=mask_square(g[2],NX,NY,device)
            with torch.no_grad():q=np.asarray([solve_one(d,j,ph,float(v),device,ORDER,masks[g[2]],g) for j in range(len(d['w']))])
            cache[key]=(metrics(d,q),q)
        return cache[key]
    def evaluate(g,r):
        if g not in records:
            met={(ph,v):state(g,ph,v)[0] for ph in ('I','M') for v in (-VEND,VEND)}
            row=dict(ID=len(records)+1,Round=r,Mode=g[0],Spacer_nm=g[1],Width_nm=g[2],VO2_nm=g[3],**objectives(met));records[g]=row
            print(f"ID {row['ID']:02d} {g}: Tmin={row['Tmin']:.4f}, margins I/M={row['Margin_I']:+.6f}/{row['Margin_M']:+.6f}, elapsed={time.perf_counter()-started:.1f}s",flush=True)
        return records[g]
    winners={}
    for mode in ('film','patch'):
        center=baseline[mode];evaluate(center,0)
        for r in range(1,LOCAL_ROUNDS+1):
            for g in neighbors(center):evaluate(g,r)
            pool=pd.DataFrame([v for k,v in records.items() if k[0]==mode]);best=choose(pool)
            new=(mode,float(best.Spacer_nm),float(best.Width_nm),float(best.VO2_nm));history.append(dict(Mode=mode,Round=r,**{k:best[k] for k in ['ID','Spacer_nm','Width_nm','VO2_nm','Tmin','WorstMargin','Qualified']}))
            print(f'Round {r} best {mode}: {new}; worst margin={best.WorstMargin:+.6f}',flush=True)
            if new==center:break
            center=new
        winners[mode]=center
    table=pd.DataFrame(records.values());best=choose(table);selected=(best.Mode,float(best.Spacer_nm),float(best.Width_nm),float(best.VO2_nm))
    detail=[];comparison=[];bounds=[]
    for mode,g in winners.items():
        for ph in ('I','M'):
            for v in np.unique(np.r_[np.linspace(-VEND,VEND,DETAIL_POINTS),0.]):
                met,_=state(g,ph,float(v));detail.append(dict(Mode=mode,Phase=ph,Voltage=v,**met))
        b=records[baseline[mode]];a=records[g]
        for label,row in [('Baseline',b),('LocalBest',a)]:comparison.append(dict(Role=label,**row))
        for parameter,value,values in [('SiO2',g[1],SPACERS),('Width',g[2],WIDTHS),('VO2Height',g[3],HEIGHTS)]:bounds.append(dict(Mode=mode,Parameter=parameter,Value_nm=value,AtLower=value==min(values),AtUpper=value==max(values)))
    detail=pd.DataFrame(detail);comparison=pd.DataFrame(comparison);bounds=pd.DataFrame(bounds);history=pd.DataFrame(history);sample=[]
    for mode in winners:
        t=detail[detail.Mode==mode];i=t[t.Phase=='I'].set_index('Voltage');m=t[t.Phase=='M'].set_index('Voltage');vs=i.index.to_numpy();am=np.maximum(np.minimum(i.DA.to_numpy(),-m.DA.to_numpy()),np.minimum(-i.DA.to_numpy(),m.DA.to_numpy()));idx=int(np.argmax(am))
        sample.append(dict(Mode=mode,MinTlum_sample=t.Tlum.min(),Transparent_sample=bool(t.Tlum.min()>=TLUM_FLOOR),BestPhaseAbsMargin_sample=am[idx],BestPhaseBias_sample=vs[idx],PhaseTsolContrast_sample=np.max(abs(i.TsolF.to_numpy()-m.TsolF.to_numpy()))))
    sample=pd.DataFrame(sample)
    for title,t in [('ALL TRIALS',table),('ROUND HISTORY',history),('BASELINE / LOCAL BEST',comparison),('VOLTAGE SAMPLES',detail),('INTERIOR SAMPLE CHECK',sample),('SEARCH BOUNDARIES',bounds)]:print('\n'+title+'\n'+t.to_string(index=False,float_format=lambda x:f'{x:.6g}'))
    print('\nPREFERRED ENDPOINT CANDIDATE\n'+best.to_string());print('Geometry:',selected)
    if not best.Qualified:print('No preferred candidate meets BOTH phase margins and endpoint transparency; this is a diagnostic candidate.')
    if not bool(sample[sample.Mode==best.Mode].Transparent_sample.iloc[0]):print('Preferred endpoint candidate violates transparency at an interior sample; it is NOT full-sample qualified.')
    plt.rcParams.update({'font.size':FONT_SIZE,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in'})
    fig,axs=plt.subplots(2,3,figsize=(16,9),constrained_layout=True);ax=axs.ravel();colors={'film':'tab:blue','patch':'tab:orange'}
    for mode in winners:
        t=table[table.Mode==mode];ax[0].scatter(t.Tmin,t.WorstMargin,color=colors[mode],label=mode)
        b=records[baseline[mode]];a=records[winners[mode]];ax[0].plot([b['Tmin'],a['Tmin']],[b['WorstMargin'],a['WorstMargin']],'-*',color=colors[mode],ms=11)
        for ph,ls in [('I','-'),('M','--')]:
            t=detail[(detail.Mode==mode)&(detail.Phase==ph)];lab=mode+' '+ph
            ax[2].plot(t.Voltage,t.DA,ls+'o',color=colors[mode],label=lab);ax[3].plot(t.Voltage,t.Tlum,ls+'o',color=colors[mode],label=lab);ax[4].plot(t.Voltage,t.TsolF,ls+'o',color=colors[mode],label=lab)
        t=detail[detail.Mode==mode];i=t[t.Phase=='I'].set_index('Voltage');m=t[t.Phase=='M'].set_index('Voltage');z=np.maximum(np.minimum(i.DA,-m.DA),np.minimum(-i.DA,m.DA));ax[5].plot(z.index,z.values,'-o',color=colors[mode],label=mode)
    x=np.arange(4);sub=comparison.reset_index(drop=True);ax[1].bar(x-.17,sub.Margin_I,.34,label='I');ax[1].bar(x+.17,sub.Margin_M,.34,label='M');ax[1].set(xticks=x,xticklabels=[r.Mode+'\n'+r.Role for r in sub.itertuples()])
    ax[0].axvline(TLUM_FLOOR,color='.5',ls=':');ax[0].axhline(MARGIN,color='.5',ls=':');ax[0].set(xlabel='Min endpoint luminous T',ylabel='Weaker-phase reversal margin');ax[0].set_title('(a) Local trials and baseline-to-best',loc='left')
    ax[1].axhline(MARGIN,color='.5',ls=':');ax[1].set_ylabel('Reversal margin');ax[1].set_title('(b) Baselines versus local candidates',loc='left')
    for j,title,y in [(2,'Electrical directional response','Solar-weighted A_F - A_B'),(3,'Minimum directional luminous T','D65 weighted T'),(4,'Forward solar transmittance','Solar weighted T'),(5,'Fixed-bias phase reversal','Either-orientation phase margin')]:ax[j].set(xlabel='Voltage (V)',ylabel=y);ax[j].set_title(f'({chr(97+j)}) '+title,loc='left')
    for j in (2,5):ax[j].axhline(0,color='.5');ax[j].axhline(MARGIN,color='.5',ls=':')
    ax[2].axhline(-MARGIN,color='.5',ls=':');ax[3].axhline(TLUM_FLOOR,color='.5',ls=':')
    for a in ax:a.legend(fontsize=8);a.grid(alpha=.15)
    print(f'\nElapsed={time.perf_counter()-started:.1f}s; unique geometries={len(records)}; ONE figure; no exports.')
    print('Bounded two-round coordinate search, NOT a full factorial/global optimum. N=2, coarse spectra: screening only. Margins allow either orientation.')
    print('Selection uses four endpoints; interior samples separately reported. Fixed-bias phase reversal now checks intermediate sampled voltages too; no continuous-voltage proof.')
    print('No Au. Other layers: PMMA15/ITO5/DAST95/ITO5; P400. Reverse source is SiO2 half-space. No thermal dynamics or thickness-dependent nk.')
    if SHOW_FIGURE:plt.show()
    return dict(AllTable=table,HistoryTable=history,ComparisonTable=comparison,DetailTable=detail,SampleTable=sample,BoundaryTable=bounds,Selected=best,Figure=fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path);p.add_argument('--no-show',action='store_true');p.add_argument('--one-round',action='store_true');a=p.parse_args()
    if a.project is not None:PROJECT=a.project.resolve()
    if a.no_show:SHOW_FIGURE=False
    if a.one_round:LOCAL_ROUNDS=1
    Results=main();globals().update(Results)
