%% build_COMSOL56_Air_HfO2_SiO2.m
% COMSOL Multiphysics 5.6 + Wave Optics + LiveLink for MATLAB
% Baseline validation model:
%   Air -> HfO2(45 nm) -> SiO2(semi-infinite, represented by output port)
% Normal incidence, 2D, TE (Ez), periodic side boundaries.
%
% Put this .m file in the SAME folder as:
%   2SiO2.xlsx
%   8HfO2.xlsx
%
% Recommended launch method:
%   Windows Start Menu -> COMSOL Multiphysics 5.6 with MATLAB
% Then in MATLAB:
%   cd('folder_containing_this_file')
%   build_COMSOL56_Air_HfO2_SiO2
%
% Output:
%   COMSOL56_Air_HfO2_45nm_SiO2.mph
%   COMSOL56_Air_HfO2_45nm_SiO2_results.csv
%   one MATLAB figure comparing COMSOL port-S R/T/A with analytic TMM

clc; close all;
import com.comsol.model.*
import com.comsol.model.util.*

%% ---------------- User settings ----------------
root = fileparts(mfilename('fullpath'));
fileSiO2 = fullfile(root,'2SiO2.xlsx');
fileHfO2 = fullfile(root,'8HfO2.xlsx');
outMph = fullfile(root,'COMSOL56_Air_HfO2_45nm_SiO2.mph');
outCsv = fullfile(root,'COMSOL56_Air_HfO2_45nm_SiO2_results.csv');

lambda_nm = (300:5:2500)';       % change 5 -> 2 nm later for a denser final sweep
dHf = 45e-9;                    % HfO2 thickness
W = 200e-9;                     % periodic cell width; planar result is independent of W
hAir = 400e-9;                  % distance from film to input port
hSub = 400e-9;                  % SiO2 section before output port
meshMax = 2e-9;                 % conservative baseline mesh

assert(exist(fileSiO2,'file')==2,'Missing file: %s',fileSiO2);
assert(exist(fileHfO2,'file')==2,'Missing file: %s',fileHfO2);

%% ---------------- Read n,k data ----------------
S = readmatrix(fileSiO2,'Sheet','Sheet1');
H = readmatrix(fileHfO2,'Sheet','Sheet1');
S = S(:,1:3); H = H(:,1:3);
S = S(all(isfinite(S),2),:); H = H(all(isfinite(H),2),:);
S = sortrows(S,1); H = sortrows(H,1);
assert(min(S(:,1))<=min(lambda_nm)/1000 && max(S(:,1))>=max(lambda_nm)/1000,'SiO2 data do not cover sweep.');
assert(min(H(:,1))<=min(lambda_nm)/1000 && max(H(:,1))>=max(lambda_nm)/1000,'HfO2 data do not cover sweep.');

%% ---------------- New COMSOL model ----------------
ModelUtil.clear;
model = ModelUtil.create('Model');
model.modelPath(root);
model.label('Air_HfO2_45nm_SiO2_baseline.mph');
comp = model.component.create('comp1',true);
geom = comp.geom.create('geom1',2);
geom.lengthUnit('m');

% Three stacked full-width rectangles. y increases upward.
r1 = geom.create('air','Rectangle');
r1.set('size',[W hAir]); r1.set('pos',[-W/2 0]);
r2 = geom.create('hfo2','Rectangle');
r2.set('size',[W dHf]); r2.set('pos',[-W/2 -dHf]);
r3 = geom.create('sio2','Rectangle');
r3.set('size',[W hSub]); r3.set('pos',[-W/2 -(dHf+hSub)]);
geom.run;

%% ---------------- Robust coordinate selections ----------------
tol = 1e-12;
domAir  = mphselectbox(model,'geom1',[-W/2-tol W/2+tol; -tol hAir+tol],'domain');
domHf   = mphselectbox(model,'geom1',[-W/2-tol W/2+tol; -dHf-tol tol],'domain');
domSiO2 = mphselectbox(model,'geom1',[-W/2-tol W/2+tol; -(dHf+hSub)-tol -dHf+tol],'domain');
bTop    = mphselectbox(model,'geom1',[-W/2-tol W/2+tol; hAir-tol hAir+tol],'boundary');
bBot    = mphselectbox(model,'geom1',[-W/2-tol W/2+tol; -(dHf+hSub)-tol -(dHf+hSub)+tol],'boundary');
bLeft   = mphselectbox(model,'geom1',[-W/2-tol -W/2+tol; -(dHf+hSub)-tol hAir+tol],'boundary');
bRight  = mphselectbox(model,'geom1',[ W/2-tol  W/2+tol; -(dHf+hSub)-tol hAir+tol],'boundary');
assert(numel(domAir)==1 && numel(domHf)==1 && numel(domSiO2)==1,'Domain selection failed.');
assert(numel(bTop)==1 && numel(bBot)==1,'Top/bottom boundary selection failed.');
assert(~isempty(bLeft) && ~isempty(bRight),'Periodic-side boundary selection failed.');
fprintf('Domains: Air=%d, HfO2=%d, SiO2=%d\n',domAir,domHf,domSiO2);
fprintf('Ports: top=%d, bottom=%d | side boundaries=%d\n',bTop,bBot,numel([bLeft(:);bRight(:)]));

%% ---------------- Interpolation functions ----------------
add_interp(model,'nSi','n_SiO2',S(:,1),S(:,2));
add_interp(model,'kSi','k_SiO2',S(:,1),S(:,3));
add_interp(model,'nHf','n_HfO2',H(:,1),H(:,2));
add_interp(model,'kHf','k_HfO2',H(:,1),H(:,3));

%% ---------------- Materials ----------------
% COMSOL harmonic convention: passive nk data are represented as N=n-i*k.
% Therefore epsilon_r=(n-i*k)^2.
matAir = comp.material.create('matAir','Common');
matAir.label('Air'); matAir.selection.set(domAir);
matAir.propertyGroup('def').set('relpermittivity',{'1'});
matAir.propertyGroup('def').set('relpermeability',{'1'});
matAir.propertyGroup('def').set('electricconductivity',{'0[S/m]'});

matHf = comp.material.create('matHf','Common');
matHf.label('HfO2 from 8HfO2.xlsx'); matHf.selection.set(domHf);
matHf.propertyGroup('def').set('relpermittivity',{'(nHf(c_const/freq)-i*kHf(c_const/freq))^2'});
matHf.propertyGroup('def').set('relpermeability',{'1'});
matHf.propertyGroup('def').set('electricconductivity',{'0[S/m]'});

matSi = comp.material.create('matSi','Common');
matSi.label('SiO2 from 2SiO2.xlsx'); matSi.selection.set(domSiO2);
matSi.propertyGroup('def').set('relpermittivity',{'(nSi(c_const/freq)-i*kSi(c_const/freq))^2'});
matSi.propertyGroup('def').set('relpermeability',{'1'});
matSi.propertyGroup('def').set('electricconductivity',{'0[S/m]'});

%% ---------------- Wave Optics physics ----------------
ewfd = comp.physics.create('ewfd','ElectromagneticWavesFrequencyDomain','geom1');
ewfd.label('Electromagnetic Waves, Frequency Domain (TE)');
% Force true 2D TE: Ez only.
try
    ewfd.prop('components').set('components','outofplane');
catch ME
    warning('Could not force out-of-plane mode via API: %s\nModel will still be excited with E0=[0 0 1].',ME.message);
end

% Periodic input/output ports. Top port is excited; bottom port is matched listener.
p1 = ewfd.create('port1','Port',1); p1.selection.set(bTop);
p1.set('PortType','Periodic');
try, p1.set('InputType','E'); catch, end
p1.set('E0',{'0','0','1'});
try, p1.set('PortExcitation','on'); catch, end

p2 = ewfd.create('port2','Port',1); p2.selection.set(bBot);
p2.set('PortType','Periodic');
try, p2.set('InputType','E'); catch, end
p2.set('E0',{'0','0','1'});
try, p2.set('PortExcitation','off'); catch, end

% Floquet periodic side boundaries, normal incidence -> tangential k = 0.
pc = ewfd.create('pc1','PeriodicCondition',1);
pc.selection.set(unique([bLeft(:);bRight(:)])');
try
    pc.set('PeriodicType','Floquet');
    pc.set('kFloquet',{'0','0','0'});
catch ME
    warning('Floquet explicit settings were not accepted; using PeriodicCondition defaults: %s',ME.message);
end

%% ---------------- Mesh ----------------
mesh = comp.mesh.create('mesh1','geom1');
sz = mesh.create('size1','Size');
sz.set('custom','on');
sz.set('hmax',sprintf('%.12g[m]',meshMax));
sz.set('hmin','0.25[nm]');
sz.set('hgrad','1.2');
% Mapped quads are ideal for a planar layered benchmark and give matching side meshes.
map = mesh.create('map1','Map');
map.selection.geom('geom1',2); map.selection.set([domAir domHf domSiO2]);
mesh.run;

%% ---------------- Frequency sweep ----------------
c0 = 299792458;
freq_Hz = c0./(lambda_nm*1e-9);
freqList = sprintf('%.15g ',freq_Hz);
std = model.study.create('std1');
fq = std.create('freq','Frequency');
fq.set('plist',strtrim(freqList));

% Save a runnable pre-solve model first.
mphsave(model,outMph);
fprintf('\nModel created and saved:\n  %s\n',outMph);
fprintf('Solving %d wavelengths (%.0f-%.0f nm)...\n',numel(lambda_nm),min(lambda_nm),max(lambda_nm));

%% ---------------- Solve ----------------
tStart = tic;
std.run;
fprintf('COMSOL solve finished in %.1f s.\n',toc(tStart));
mphsave(model,outMph);

%% ---------------- Extract Port-S results ----------------
fsol = mphglobal(model,'freq');
Rcom = real(mphglobal(model,'abs(ewfd.S11)^2'));
Tcom = real(mphglobal(model,'abs(ewfd.S21)^2'));
lambdaCom_nm = c0./fsol(:)*1e9;
Rcom = Rcom(:); Tcom = Tcom(:);
[lambdaCom_nm,ord] = sort(lambdaCom_nm); Rcom=Rcom(ord); Tcom=Tcom(ord);
Acom = 1-Rcom-Tcom;

%% ---------------- Independent analytic TMM ----------------
lambda_um = lambdaCom_nm/1000;
nH = interp1(H(:,1),H(:,2),lambda_um,'linear'); kH = interp1(H(:,1),H(:,3),lambda_um,'linear');
nS = interp1(S(:,1),S(:,2),lambda_um,'linear'); kS = interp1(S(:,1),S(:,3),lambda_um,'linear');
N0 = ones(size(lambda_um)); N1 = nH-1i*kH; Ns = nS-1i*kS;
r01=(N0-N1)./(N0+N1); r12=(N1-Ns)./(N1+Ns);
beta=2*pi.*N1*dH./(lambdaCom_nm*1e-9);
r=(r01+r12.*exp(-2i*beta))./(1+r01.*r12.*exp(-2i*beta));
t=(2*N0./(N0+N1)).*(2*N1./(N1+Ns)).*exp(-1i*beta)./(1+r01.*r12.*exp(-2i*beta));
Rtmm=abs(r).^2; Ttmm=real(Ns)./real(N0).*abs(t).^2; Atmm=1-Rtmm-Ttmm;

fprintf('\n================ VALIDATION =================\n');
fprintf('max |R_COMSOL-R_TMM| = %.4g\n',max(abs(Rcom-Rtmm)));
fprintf('max |T_COMSOL-T_TMM| = %.4g\n',max(abs(Tcom-Ttmm)));
fprintf('COMSOL max |R+T+A-1| = %.3e\n',max(abs(Rcom+Tcom+Acom-1)));
fprintf('=============================================\n');

%% ---------------- Save numeric results ----------------
Result = table(lambdaCom_nm,Rcom,Tcom,Acom,Rtmm,Ttmm,Atmm, ...
    'VariableNames',{'lambda_nm','R_COMSOL','T_COMSOL','A_COMSOL','R_TMM','T_TMM','A_TMM'});
writetable(Result,outCsv);
fprintf('Results saved:\n  %s\n',outCsv);

%% ---------------- One figure only ----------------
figure('Color','w','Position',[100 100 920 700]);
subplot(2,1,1);
plot(lambdaCom_nm,Rcom,'LineWidth',1.7); hold on;
plot(lambdaCom_nm,Rtmm,'--','LineWidth',1.5);
plot(lambdaCom_nm,Tcom,'LineWidth',1.4);
plot(lambdaCom_nm,Ttmm,'--','LineWidth',1.3);
grid on; box on; xlim([min(lambdaCom_nm) max(lambdaCom_nm)]); ylim([0 1]);
xlabel('Wavelength (nm)'); ylabel('Power coefficient');
legend('R COMSOL','R TMM','T COMSOL','T TMM','Location','best');
title('Air / HfO_2(45 nm) / SiO_2: COMSOL 5.6 vs analytic TMM');

subplot(2,1,2);
plot(lambdaCom_nm,Rcom-Rtmm,'LineWidth',1.4); hold on;
plot(lambdaCom_nm,Tcom-Ttmm,'LineWidth',1.4);
yline(0,'k:'); grid on; box on; xlim([min(lambdaCom_nm) max(lambdaCom_nm)]);
xlabel('Wavelength (nm)'); ylabel('Difference');
legend('\DeltaR','\DeltaT','Location','best');
title('Numerical difference');

%% ================= Local helper =================
function add_interp(model,tag,fname,x,y)
    f = model.func.create(tag,'Interpolation');
    f.label(fname); f.set('funcname',tag);
    f.set('interp','linear');
    try, f.set('argunit','um'); catch, end
    try, f.set('fununit','1'); catch, end
    tab = cell(numel(x),2);
    for ii=1:numel(x)
        tab{ii,1}=sprintf('%.15g',x(ii));
        tab{ii,2}=sprintf('%.15g',y(ii));
    end
    f.set('table',tab);
end
