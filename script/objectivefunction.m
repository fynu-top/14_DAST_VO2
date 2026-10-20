
function Y = objectivefunction(X)

Visible(:,1) = [100	300	400	405	410	415	420	425	430	435	440	445	450	455	460	465	470	475	480	485	490	495	500	505	510	515	520	525	530	535	540	545	550	555	560	565	570	575	580	585	590	595	600	605	610	615	620	625	630	635	640	645	650	655	660	665	670	675	680	685	690	695	700	705	710	715	720	725	730	755	790	815	830	10000];
Visible(:,2) = [0 0 0 0.001 0.001 0.002 0.004 0.007 0.012 0.017 0.023 0.03 0.038 0.048 0.06 0.074 0.091 0.113 0.139 0.169 0.208 0.259 0.323 0.407 0.503 0.608 0.71 0.793 0.862 0.915 0.954 0.98 0.995 1 0.995 0.979 0.952 0.915 0.87 0.816 0.757 0.695 0.631 0.567 0.503 0.441 0.381 0.321 0.265 0.217 0.175 0.138 0.107 0.082 0.061 0.045 0.032 0.023 0.017 0.012 0.008 0.006 0.004 0.003 0.002 0.001 0.001 0.001 0.001 0 0 0 0 0];

% objectivefunction([69 50 9 51])
% objectivefunction([73 50 10 71])
% objectivefunction([73 50 10])
% objectivefunction([10000])
% X = [5 200 26];
h = X;
% h = [73 50 10 71];
shoulian = load('shoulian.mat');shoulian = shoulian.shoulian;
nk_M = load('nk_M.mat') ;nk_M = nk_M.nk_M;
% load nk_M.mat;
% ceng_layer = length(h)+1;
% h(ceng_layer) = 500;
h = [50 h 500]; % [ 上层PMMA 中间层 底层基底];
%%
clear Absorption xz
% wave_range = (300:20:10000)'; 
wave = load('wave.mat');wave = wave.wave';
%% case 1
%{
for m = 1:length(wave_range)
    n0 = 1;
    lambda = wave_range(m)*1e-9;
    nsubs = 1; % 基底材质为空气
    n_layers = nk_M(m,:);
    d_layers = h*1e-9;
    theta_in = 0;
    [R(m,1) T(m,1)] = reflectance_TE(n0,nsubs,n_layers,d_layers,theta_in,lambda);
end
Absorption=1-R;
plot(wave_range,R);hold on;plot(wave_range,T);
%}

%% case 4
%{1
for m = 1:length(wave)
    n0 = 1;
    lambda = wave(m)*1e-9;
    nsubs = 1; % 基底材质为空气
    n_layers = nk_M(m,:);
    d_layers = h*1e-9;
    theta_in = 0;
%     [R(m,1) T(m,1)] = reflectance_TE(n0,nsubs,n_layers,d_layers,theta_in,lambda);
    [R1(m,1) T1(m,1) A1(m,1)] = tm1( wave(m), nk_M(m,:),  h  );
    [R2(m,1) T2(m,1) A2(m,1)] = tm1( wave(m), rot90(rot90(nk_M(m,:))),  rot90(rot90(h))  );
end

% figure('color',[1 1 1]); 
% plot(wave_range,R1,'b',  wave_range,T1,'r',  wave_range,A1,'g',  'LineWidth',2);hold on;
% plot(wave_range,R2,'b--',wave_range,T2,'r--',wave_range,A2,'g--','LineWidth',2);hold on;
% legend('\fontsize{20}R','\fontsize{20}T','\fontsize{20}A'); xlabel('wavelength [nm]');ylabel('intensity');
%% calculation
V = interp1(Visible(:,1),Visible(:,2),wave)';
Visible_I = V*T1/sum(V);
Target=[540 660 640 960];
for ix = 1:length(Target)
    [ipx(ix) ipy(ix)] = find(wave == Target(ix));
end
Visible_intensity = sum(T1(ipx(1):ipx(2)))/length(T1(ipx(1):ipx(2)));
A_difference =abs(sum(A1(ipx(3):ipx(4)))/length(A1(ipx(3):ipx(4)))-sum(A2(ipx(3):ipx(4)))/length(A2(ipx(3):ipx(4))) );

% if shoulian == 1
%     RTA =[R1 T1 A1 R2 T2 A2]; save('RTA.mat','RTA');
% end
%}
%% case 2
%{
for m = 1:length(wave_range)
    n0 = 1;
    lambda = wave_range(m);
%     nsubs = nk_M(m,end); % 基底材质
    nsubs = 1; % 基底材质为空气
    n_layers = nk_M(m,:);
    d_layers = h;
    theta_in = 0;
%     R_TE(m) = reflectance_TE(n0,nsubs,n_layers,d_layers,theta_in,lambda);
    [R(m,1) T(m,1)] = tmm(lambda,n_layers,d_layers);%          n0,nsubs,n_layers,d_layers,theta_in,lambda);
end
Absorption=1-R;

%}
%% case 3
%{
for xz=1:length(wave_range) %%%%%%%%%%%%%%%对波长的扫参
%             xz=wave_range(xz); 
% wave_range = (300:20:10000)';  
%             wa=wa+1;
            wa = xz;
            n0=1;a0=0;
            d=wave_range(xz);
            n2=nk_M(xz,1);%%%%%%%%%%%%%%%%%%%%%%%%
            a1=0;   a2=0;a3=0;
            
            u0=4*pi*1e-7;e0=1e-9/(36*pi);b=sqrt(e0/u0);c1=b*n2/cos(a1);%波阻抗
            c0=b*n0/cos(a0);c00=b*n0/cos(a3);
            p1=2*pi*n2*h(1)*cos(a1)/d;%相位厚度
            A1=cos(p1);B1=-1i*sin(p1)/c1;C1=-1i*c1*sin(p1);D1=cos(p1);
            M=[A1 B1;C1 D1];%介质1的传输矩阵
                for ceng=2:1:ceng_layer   %%%%%%%%%%%%%结构的层数
                    d2=h(ceng);
                    n2=nk_M(xz,ceng);
                    p2=2*pi*n2*d2*cos(a2)/d;c2=b*n2/cos(a2);
                    A2=cos(p2);B2=-1i*sin(p2)/c2;C2=-1i*c2*sin(p2);D2=cos(p2);              
                    M2=[A2 B2;C2 D2];M=M*M2;
                end
            A=M(1,1);B=M(1,2);C=M(2,1);D=M(2,2);
            r=(A*c0+B*c0*c00-C-D*c00)/(A*c0+B*c0*c00+C+D*c00);%反射系数
            t=(2*c0)/(A*c0+B*c0*c00+C+D*c00);%透射系数
            R(xz)=r*conj(r);T(xz)=t*conj(t);
%             a(wa) = 1-R(xz);
%             a(wa)=1-r*conj(r);
            AA(xz)=1-R(xz);
end
Absorption=1-R';
%}
%{
%% 求目标函数
%% 输入加进黑体辐射后的谱（与温度，太阳浓度有关）
solar=[0.001 0.2053 0.5018 0.5982 0.7008 1.1141 1.1232 1.3499 1.5291 1.6181 1.5451 1.5236 1.4825 1.474 1.502 1.4753 1.4739 1.434 1.3992 1.3969 1.2823 0.9855 1.2195 0.266 1.1636 1.0725 0.8619 1.0157 0.9882 0.9396 0.7426 0.7441 0.4718 0.4207 0.6047 0.7353 0.699 0.6717 0.6359 0.5972 0.4858 0.1419 0.256 0.2865 0.4407 0.4482 0.4581 0.4608 0.4311 0.422 0.3531 0.2587 0.1683 0 0.0001 0 0.0083 0.0396 0.0854 0.0606 0.2506 0.2645 0.2649 0.2657 0.2446 0.2381 0.2345 0.2151 0.2233 0.2056 0.1998 0.187 0.1682 0.16 0.1005 0.0318 0.001 0 0 0.0001 0 0.0005 0.0033 0.0219 0.0755 0.0382 0.045 0.0898 0.0692 0.0868 0.0861 0.0876 0.0908 0.0842 0.0818 0.0712 0.0777 0.0731 0.0669 0.0663 0.0588 0.052 0.0458 0.0502 0.0426 0.0441 0.0266 0.0432 0.0334 0.008 0.0071 0.0004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0002 0.0008 0.0029 0.0016 0.0046 0.0013 0.0078 0.0006 0.002 0.0063 0.0036 0.0044 0.0098 0.0033 0.0092 0.0106 0.0004 0.0016 0.0037 0.0012 0.0029 0.0018 0.0001 0.0035 0.0052 0.0051 0.0125 0.0132 0.008 0.0125 0.0112 0.0119 0.0121 0.009 0.0108 0.0102 0.0103 0.0116 0.0115 0.0109 0.0083 0.0109 0.0104 0.0088 0.0089 0.0096 0.0099 0.0097 0.009 0.008 0.0065 0.0079 0.0069 0.0074 0.0077 0.0074 0.0071 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
solar = interp1(300:20:10000,solar,wave_range)';
% x_wavelength = 0.3:0.02:10;
x_wavelength = wave_range'*1e-3;
weight(1) = 1;
for mm = 2:length(wave_range)
    weight(mm) = (wave_range(mm) - wave_range(mm-1))/(wave_range(2) - wave_range(1));
end
weight = weight/(sum(weight)/length(weight));
clear mm
T_t = 100 +273.15;% K
% T_t = 1000 ;% K
T_spec = 3.742./ ((x_wavelength.^5).*(exp(1.439e4./(x_wavelength*T_t))-1)) ; % gives the result in W/m2/um x 1e8
%plot(x_wavelength,T_spec);xlabel('wavelength (\mum)');ylabel('spectral intensity (W/m^2/\mum) x 10^8');title('Blackbody Radiation');legend(sprintf('T = %.0f K',T_t));
% T_spectrum=T_spectrum';
T_spec=T_spec*5.67*0.00000001*(T_t^4-300^4)/1000/1000;

    A_return =  solar*Absorption/sum(solar)-0/50*T_spec*Absorption/sum(T_spec)*5.67*1e-8*(T_t^4-300^4)/1000     ;   %% 光热转化效率
%     A_return = (solar*Absorption/sum(solar)-1/50*T_spec*Absorption/sum(T_spec)*5.67*1e-8*(T_t^4-300^4)/1000)*(1-300/T_t)    ; %% 卡诺基效率
%}



if shoulian == 1     Y = -(  Visible_intensity+ A_difference*3                                  )  ; end
if shoulian == 2     Y = [R1 T1 A1 R2 T2 A2];Y(1,7)=Visible_I;  end

% objectivefunction([61 32 9.2 46.3])







