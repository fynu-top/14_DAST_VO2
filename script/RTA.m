% jreftran_rt.m
% Reflection and transmission coefficient calculator of lossy multilayer 
% media for optical scientists

% Shawn Divitt ETH Zurich Photonics Group Ver. 2.1, 15 Nov 2016 -iwt time dependence

%This program produces the complex reflection and transmission coefficients
%of a multilayer stack given the angle of incidence, polarization,
%wavelength, complex-refractive index of each layer, and thickness of each
%layer. The program assumes lossless dielectric incident and exit media. 
%The program also assumes that all layers are non-magnetic. Magnetic media
%can be handled by more general theory (see technical report cited below). 

%This program was inspired by the technical report written by K. Pascoe, 
%"Reflectivity and Transmissivity through Layered, Lossy Media: A User-Friendly Approach," 2001. See the following link for the technical
%paper:
%http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix=html&identifier=ADA389099
% function [r,t,R,T,A]=RTA(lambda, dM, nk, IncidentAngle, polarization)
function [R T A] = RTA(lambda, nk, dM, IncidentAngle, polarization)
% [R1(m,1) T1(m,1) A1(m,1)] = RTA( wave(m), nk_M(m,:),  h  ,IncidentAngle, 0);
nk=conj(nk);
nk = [1 nk 1];
dM = [0 dM 0];

%l = free space wavelength, nm
%d = layer thickness vector, nm
%n = layer complex refractive index vector
%t0= angle of incidence
%polarization should be 0 for TE (s-polarized), otherwise TM (p-polarized)

Z0 = 120*pi; %impedance of free space, Ohms

%the line below had mistakenly been a Z0/n instead of a n/Z0 in version 1!
Y=nk./Z0; %admittance in terms of impedance of free space and refractive index, assuming non-magnetic media
g=1i*2*pi*nk/lambda; %propagation constant in terms of free space wavelength and refractive index

%all of the calculations rely on cosine of the complex angle, but we can
%only find the sine of the complex angle from snells law. So we use the
%fact that cos(asin(x))=sqrt(1-x^2)
%t=asin(n(1)./n*sin(t0)), complex theta for each layer
ct=sqrt(1-(nk(1)./nk*sin(IncidentAngle)).^2); %cosine theta

if polarization==0
    eta=Y.*ct; %tilted admittance, TE case
else
    eta=Y./ct; %tilted admittance, TM case
end

delta=1i*g.*dM.*ct;

M = zeros(2,2,length(dM));
for j=1:length(dM)
    M(1,1,j)=cos(delta(j));
    M(1,2,j)=1i./eta(j).*sin(delta(j));
    M(2,1,j)=1i*eta(j).*sin(delta(j));
    M(2,2,j)=cos(delta(j));
end    
%%
%     M1 = zeros(2,2,length(dM));
%     M1(1,1,:) = cos(delta);
%     M1(1,2,:) = 1i./eta.*sin(delta);
%     M1(2,1,:) = 1i*eta.*sin(delta);
%     M1(2,2,:) = cos(delta);
%     M = M1;
    
M_t=[1,0;0,1]; %M total
for j=2:(length(dM)-1)
    M_t=M_t*M(:,:,j);
end

r=(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))-(M_t(2,1)+M_t(2,2)*eta(end)))/(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))+(M_t(2,1)+M_t(2,2)*eta(end)));
t=2*eta(1)/(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))+(M_t(2,1)+M_t(2,2)*eta(end)));

R = abs(r)^2;
T = real(eta(end)/eta(1))*abs(t)^2;

A = (4*eta(1)*real((M_t(1,1)+M_t(1,2)*eta(end))*conj(M_t(2,1)+M_t(2,2)*eta(end))-eta(end)))/abs(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))+(M_t(2,1)+M_t(2,2)*eta(end)))^2;




