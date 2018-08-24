function FiveSteph2
clc
options.MaxFunEvals = 9000000;
x0=[608,0.3366,12,25,40,55,.33,.33,.33,.33,.57,-18,109,1.44,0,-19,73,0,102,154,1.0104,0,136,187,1.266,1.66,5,-3,7.2,1.000825,-152,-3.3,1.8];
fsolve(@fan,x0)
% tau=x(1) a5=x(2) n61=x(3) n62=x(4) n63=x(5) n64=x(6) a1=x(7) a2=x(8)
% a3=x(9) a4=x(10) psi=x(11) tD=x(12) tB=x(13) theta=x(14) psi2=x(15)
% tD2=x(16) tB2=x(17) psi3=x(18) tD2A=x(19) tB2A=x(20) theta2A=x(21) 
% psi3A=x(22) tD3A=x(23) tB3A=x(24) theta3A=x(25) theta2B=x(26)
% psi3B=x(27) tD3B=x(28) tB3B=x(29) theta3B=x(30) tD4B=x(31) tB4B=x(32)
% theta4B=x(33)	
 
function y=fan(x)
options.MaxFunEvals = 9000000;
 
k10=1.04E+06; % Pre-exponential Factor
k20=7.30E+15;
k30=8.62E+11;
E1=3.43E+04; % Activation Energy
E2=1.11E+05;
E3=8.47E+04;
%---ADJUSTABLE PARAMETERS-----------------SPECIFY------------%
cat=3*10^-5; % concentration of catalyst, mol/L
p=800; % rxr total pressure, psig
ncyc=150; % mol/s cyclohexane
t=80; % rxr temperature, degC
%-----------------------------------------------------%
T=t+273.15; % temperature, kelvin
 
% vapor pressure calculations
pcyc=0.019337*(exp(15.7527-2766.63/(T+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p2=(-0.0573)*(t^2)+21.523*t+461.88; %vapor P,ethylene,psi (depends on t)
p6=0.019337*exp(15.8089-2654.81/(T-47.3)); %vapor P, hexene, psi
p8=0.019337*exp(15.963-3116.52/(T-60.39)); %vapor P, octene
p10=0.019337*exp(16.0129-3448.18/(T-76.09)); %vapor P, decene
poct=0.019337*exp(15.7428-3017.81/(T-137.1));%vapor P, octanol, psi 
peth=0.019337*exp(15.6637-1511.42/(T-17.16)); %vapor P, ethane 
pmeth=0.019337*exp(15.2243-597.84/(T-7.16)); %vapor P, methane
 
k1=cat*k10*exp(-E1/(8.3145*(t+273.15)));
k2=cat*k20*exp(-E2/(8.3145*(t+273.15)));
k3=cat*k30*exp(-E3/(8.3145*(t+273.15)));
 
n20=ncyc*((pcyc-p)/(p-p2)); % calculates inlet ethylene feed (mol/s) to rxr 1 based 
                            % on set pressure and set cyclohex feed
% conversion per reactor                          
a1=x(7); 
a2=x(8);
a3=x(9);
a4=x(10);
a5=x(2);
 
% RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR RXR
% 1 REACTOR------------------------------------------------------------
n21=n20*(1-a1); 
y(1)=n21-((n20)/x(1))/(3*k1+4*k2+2*k1*k3*x(1)/(1+k3*x(1))+1/x(1)); %rxr 1, ethylene tau
y(2)=x(3)-k1*x(1)*n21/(1+k3*x(1));                                 %rxr 1, hexene tau
n81=x(1)*k2*n21; % octene outlet from rxr1, mol/s
n101=x(1)*k3*((k1*x(1)*n21)/(1+k3*x(1))); % decene outlet rxr1, mol/s
% 1 NEW CODE BELOW
n21_=((p6*x(3)+p10*n101+p8*n81+pcyc*ncyc)-(p*(x(3)+n101+n81+ncyc)))/(p-p2); % ethylene inlet to rxr 2
n21m=n21_-n21; % makeup ethylene feed to rxr 2
sum1=n21+x(3)+n101+n81+ncyc;
f21=n21/sum1; f61=x(3)/sum1; f101=n101/sum1; f81=n81/sum1; fcyc1=ncyc/sum1; f1=f21+f61+f101+f81+fcyc1;
P1before=p2*f21+p6*f61+p10*f101+p8*f81+pcyc*fcyc1;
sum1=n21_+x(3)+n101+n81+ncyc;
f21=n21_/sum1; f61=x(3)/sum1; f101=n101/sum1; f81=n81/sum1; fcyc1=ncyc/sum1; f1=f21+f61+f101+f81+fcyc1;
P1after=p2*f21+p6*f61+p10*f101+p8*f81+pcyc*fcyc1;
% 1 END OF NEW CODE
 
% 2 REACTOR------------------------------------------------------------
n22=n21_*(1-a2); 
y(3)=n22-(n21_/x(1)-2*k3*x(3)/(1+k3*x(1)))/(3*k1+4*k2+2*k1*k3*x(1)/(1+k3*x(1))+1/x(1));
y(4)=x(4)-(k1*x(1)*n22+x(3))/(1+k3*x(1));
n82=x(1)*k2*n22+n81; % octene out of rxr2, mol/s
n102=x(1)*k3*((k1*x(1)*n22+x(3))/(1+k3*x(1)))+n101; % decene outlet rxr2, mol/s
% 2 NEW CODE BELOW- calculates makeup ethylene feed for next reactor
n22_=((p6*x(4)+p10*n102+p8*n82+pcyc*ncyc)-(p*(x(4)+n102+n82+ncyc)))/(p-p2); % ethylene inlet to rxr 3
n22m=n22_-n22; % makeup ethylene feed to rxr 3
sum2=n22+x(4)+n102+n82+ncyc;
f22=n22/sum2; f62=x(4)/sum2; f102=n102/sum2; f82=n82/sum2; fcyc2=ncyc/sum2; f2=f22+f62+f102+f82+fcyc2;
P2before=p2*f22+p6*f62+p10*f102+p8*f82+pcyc*fcyc2;
sum2=n22_+x(4)+n102+n82+ncyc;
f22=n22_/sum2; f62=x(4)/sum2; f102=n102/sum2; f82=n82/sum2; fcyc2=ncyc/sum2; f2=f22+f62+f102+f82+fcyc2;
P2after=p2*f22+p6*f62+p10*f102+p8*f82+pcyc*fcyc2;
% 2 END OF NEW CODE
 
% 3 REACTOR------------------------------------------------------------
n23=n22_*(1-a3);
y(5)=n23-(n22_/x(1)-2*k3*x(4)/(1+k3*x(1)))/(3*k1+4*k2+2*k1*k3*x(1)/(1+k3*x(1))+1/x(1));
y(6)=x(5)-(k1*x(1)*n23+x(4))/(1+k3*x(1));
n83=x(1)*k2*n23+n82;
n103=x(1)*k3*((k1*x(1)*n23+x(4))/(1+k3*x(1)))+n102;
% 3 NEW CODE BELOW
n23_=((p6*x(5)+p10*n103+p8*n83+pcyc*ncyc)-(p*(x(5)+n103+n83+ncyc)))/(p-p2); % ethylene inlet to rxr 4
n23m=n23_-n23; % makeup ethylene feed to rxr 4
sum3=n23+x(5)+n103+n83+ncyc;
f23=n23/sum3; f63=x(5)/sum3; f103=n103/sum3; f83=n83/sum3; fcyc3=ncyc/sum3; f3=f23+f63+f103+f83+fcyc3;
P3before=p2*f23+p6*f63+p10*f103+p8*f83+pcyc*fcyc3;
sum3=n23_+x(5)+n103+n83+ncyc;
f23=n23_/sum3; f63=x(5)/sum3; f103=n103/sum3; f83=n83/sum3; fcyc3=ncyc/sum3; f3=f23+f63+f103+f83+fcyc3;
P3after=p2*f23+p6*f63+p10*f103+p8*f83+pcyc*fcyc3;
% 3 END OF NEW CODE
 
% 4 REACTOR------------------------------------------------------------
n24=n23_*(1-a4); 
y(7)=n24-(n23_/x(1)-2*k3*x(5)/(1+k3*x(1)))/(3*k1+4*k2+2*k1*k3*x(1)/(1+k3*x(1))+1/x(1));
y(8)=x(6)-(k1*x(1)*n24+x(5))/(1+k3*x(1));
n84=x(1)*k2*n24+n83;
n104=x(1)*k3*((k1*x(1)*n24+x(5))/(1+k3*x(1)))+n103;
% 4 NEW CODE BELOW
n24_=((p6*x(6)+p10*n104+p8*n84+pcyc*ncyc)-(p*(x(6)+n104+n84+ncyc)))/(p-p2); % ethylene inlet to rxr 5
n24m=n24_-n24; % makeup ethylene feed to rxr 4
sum4=n24+x(6)+n104+n84+ncyc;
f24=n24/sum4; f64=x(6)/sum4; f104=n104/sum4; f84=n84/sum4; fcyc4=ncyc/sum4; f4=f24+f64+f104+f84+fcyc4;
P4before=p2*f24+p6*f64+p10*f104+p8*f84+pcyc*fcyc4;
sum4=n24_+x(6)+n104+n84+ncyc;
f24=n24_/sum4; f64=x(6)/sum4; f104=n104/sum4; f84=n84/sum4; fcyc4=ncyc/sum4; f4=f24+f64+f104+f84+fcyc4;
P4after=p2*f24+p6*f64+p10*f104+p8*f84+pcyc*fcyc4;
% 4 END OF NEW CODE
 
% 5 REACTOR------------------------------------------------------------
n25=n24_*(1-a5);
y(9)=n25-(n24_/x(1)-2*k3*x(6)/(1+k3*x(1)))/(3*k1+4*k2+2*k1*k3*x(1)/(1+k3*x(1))+1/x(1));
y(10)=72-(k1*x(1)*n25+x(6))/(1+k3*x(1));
n85=x(1)*k2*n25+n84;
n105=x(1)*k3*((k1*x(1)*n25+x(6))/(1+k3*x(1)))+n104;
% 5 NEW CODE BELOW
%n25_=((p6*72+p10*n105+p8*n85+pcyc*ncyc)-(p*(72+n105+n85+ncyc)))/(p-p2); % ethylene inlet to DISTILLATION SEQUENCE
%n25m=n25_-n25; % makeup ethylene feed to DISTILLATION, not necessary!
sum5=n25+72+n105+n85+ncyc;
f25=n25/sum5; f65=72/sum5;   f105=n105/sum5; f85=n85/sum5; fcyc5=ncyc/sum5; f5=f25+f65+f105+f85+fcyc5;
P5before=p2*f25+p6*f65+p10*f105+p8*f85+pcyc*fcyc5;
%sum5=n25_+72+n105+n85+ncyc;
%f25=n25_/sum5; f65=72/sum5;   f105=n105/sum5; f85=n85/sum5; fcyc5=ncyc/sum5; f5=f25+f65+f105+f85+fcyc5;
%P5after=p2*f25+p6*f65+p10*f105+p8*f85+pcyc*fcyc5;
% 5 END OF NEW CODE
%-------------------------------------------------------------------------
 
n2total=(n20+n21m+n22m+n23m+n24m); % total amount of ethylene (mol/s) inlet to rxr train (fresh+recycle)
fresh=(n2total)-n25;               % fresh c2 reqd, assumes 100% of outlet ethylene (pure) is recycled
recycle=n25;                       % c2 recycle, assumes 100% of outlet ethylene (pure) is recycled
 
A=((n2total)-n25)/((n2total)); % overall conversion of ethylene
s6=3*72/(n2total-n25);         % overall hexene selectivity
s10=(5*n105)/(n2total-n25);    % selectivity of decene
s8=1-s10-s6;                   % selectivity octene
 
y6=s6*A;         % yield of hexene (mol percent)
y8=s8*A;         % yield octene
y10=s10*A;       % yield decene
hours=x(1)/3600; % tau, hours
 
%---------convert from molar flow to mass flow rates-----------------------
dcyc=779000; % liquid density cyclohexane, g/m3
d6=673000; 
d8=715000;
d10=740000;
d2=567920; % this is LIQUID density of ethylene
mwcyc=84.16; % molar mass cyclohexane, g/mol
mw6=84.1608; 
mw8=112.24;
mw10=140.27;
mw2=28.05;
mweth=30.07;
mwmeth=16.04;
mwoct=130.2279;
 
% lb/h                             m3/s
mcyc=ncyc*mwcyc*3600*0.00220462;   vcyc=(ncyc*mwcyc)/dcyc;  % cyclohexane
m61=x(3)*mw6*3600*0.00220462;      v61=(x(3)*mw6)/d6; % hexene
m62=x(4)*mw6*3600*0.00220462;      v62=(x(4)*mw6)/d6;
m63=x(5)*mw6*3600*0.00220462;      v63=(x(5)*mw6)/d6;
m64=x(6)*mw6*3600*0.00220462;      v64=(x(6)*mw6)/d6;
m65=72*mw6*3600*0.00220462;        v65=(72*mw6)/d6;
m81=n81*mw8*3600*0.00220462;       v81=(n81*mw8)/d8; % octene
m82=n82*mw8*3600*0.00220462;       v82=(n82*mw8)/d8;
m83=n83*mw8*3600*0.00220462;       v83=(n83*mw8)/d8;
m84=n84*mw8*3600*0.00220462;       v84=(n84*mw8)/d8;
m85=n85*mw8*3600*0.00220462;       v85=(n85*mw8)/d8;
m101=n101*mw10*3600*0.00220462;    v101=(n101*mw10)/d10; % decene 
m102=n102*mw10*3600*0.00220462;    v102=(n102*mw10)/d10;
m103=n103*mw10*3600*0.00220462;    v103=(n103*mw10)/d10;
m104=n104*mw10*3600*0.00220462;    v104=(n104*mw10)/d10;
m105=n105*mw10*3600*0.00220462;    v105=(n105*mw10)/d10;
m20=n20*mw2*3600*0.00220462;       v20=(n20*mw2)/d2; % ethylene
m21=n21*mw2*3600*0.00220462;       v21=(n21*mw2)/d2; 
m22=n22*mw2*3600*0.00220462;       v22=(n22*mw2)/d2;
m23=n23*mw2*3600*0.00220462;       v23=(n23*mw2)/d2;
m24=n24*mw2*3600*0.00220462;       v24=(n24*mw2)/d2;
m25=n25*mw2*3600*0.00220462;       v25=(n25*mw2)/d2;
m21m=n21m*mw2*3600*0.00220462;       v21m=(n21m*mw2)/d2; 
m22m=n22m*mw2*3600*0.00220462;       v22m=(n22m*mw2)/d2;
m23m=n23m*mw2*3600*0.00220462;       v23m=(n23m*mw2)/d2;
m24m=n24m*mw2*3600*0.00220462;       v24m=(n24m*mw2)/d2;
m21_=n21_*mw2*3600*0.00220462;       v21_=(n21_*mw2)/d2; 
m22_=n22_*mw2*3600*0.00220462;       v22_=(n22_*mw2)/d2;
m23_=n23_*mw2*3600*0.00220462;       v23_=(n23_*mw2)/d2;
m24_=n24_*mw2*3600*0.00220462;       v24_=(n24_*mw2)/d2;
m2total=n2total*mw2*3600*0.00220462; v2total=(n2total*mw2)/d2;
 
 
fprintf('\n-------------------MOL/S BASIS----------------------------------------------------------\n')
fprintf('Stream    Ethylene    Hexene      Decene     Octene     Cyclohexane\n')
fprintf('   0      %f    %f    %f   %f   %f   mol/s\n',0,0,0,0,ncyc) % stream 0, mol/s
fprintf('   1      %f   %f   %f   %f   %f   mol/s\n',n21,x(3),n101,n81,ncyc) % stream 1
fprintf('   2      %f   %f   %f   %f   %f   mol/s\n',n22,x(4),n102,n82,ncyc) % stream 2
fprintf('   3      %f   %f   %f   %f   %f   mol/s\n',n23,x(5),n103,n83,ncyc) % stream 3
fprintf('   4      %f   %f   %f   %f   %f   mol/s\n',n24,x(6),n104,n84,ncyc) % stream 4
fprintf('   5      %f  %f   %f   %f   %f   mol/s\n',n25,72,n105,n85,ncyc) % stream 5
fprintf('   0m     %f  mol/s makeup\n',n20)
fprintf('   1m     %f   mol/s makeup\n',n21m)
fprintf('   2m     %f   mol/s makeup\n',n22m)
fprintf('   3m     %f   mol/s makeup\n',n23m)
fprintf('   4m     %f   mol/s makeup\n\n',n24m)
fprintf('Reactor   TotalC2Inlet    TotalC2Outlet    Conversion\n')
fprintf('   1      %f      %f        %f\n',n20,n21,x(7))
fprintf('   2      %f      %f        %f\n',n21_,n22,x(8))
fprintf('   3      %f      %f        %f\n',n22_,n23,x(9))
fprintf('   4      %f      %f        %f\n',n23_,n24,x(10))
fprintf('   5      %f      %f       %f\n',n24_,n25,x(2))
fprintf('overall                                    %f\n\n',A)
fprintf('Total C2 into   RXR train:   %f mol/s\n',n2total)
fprintf('Total C2 out of RXR train:   %f mol/s\n',n25)
fprintf('Minimum Fresh C2 to RXRs:    %f mol/s\n\n',fresh)
fprintf('------------------------LB/H BASIS-----------------------------------------------------\n')
fprintf('Stream    Ethylene       Hexene        Decene       Octene     Cyclohexane\n')
fprintf('   0      %f       %f       %f    %f      %f   lb/h\n',0,0,0,0,mcyc) % stream 0, mol/s
fprintf('   1      %f   %f    %f   %f   %f   lb/h\n',m21,m61,m101,m81,mcyc) % stream 1
fprintf('   2      %f   %f   %f   %f    %f   lb/h\n',m22,m62,m102,m82,mcyc) % stream 2
fprintf('   3      %f   %f   %f   %f    %f   lb/h\n',m23,m63,m103,m83,mcyc) % stream 3
fprintf('   4      %f   %f   %f   %f    %f   lb/h\n',m24,m64,m104,m84,mcyc) % stream 4
fprintf('   5      %f   %f   %f   %f   %f   lb/h\n',m25,m65,m105,m85,mcyc) % stream 5
fprintf('   0m     %f   lb/h makeup\n',m20)
fprintf('   1m     %f   lb/h makeup\n',m21m)
fprintf('   2m     %f   lb/h makeup\n',m22m)
fprintf('   3m     %f   lb/h makeup\n',m23m)
fprintf('   4m     %f   lb/h makeup\n\n',m24m)
fprintf('Reactor   TotalC2Inlet    TotalC2Outlet    Conversion\n')
fprintf('   1      %f      %f    %f\n',m20,m21,x(7))
fprintf('   2      %f      %f    %f\n',m21_,m22,x(8))
fprintf('   3      %f      %f    %f\n',m22_,m23,x(9))
fprintf('   4      %f      %f    %f\n',m23_,m24,x(10))
fprintf('   5      %f      %f    %f\n',m24_,m25,x(2))
fprintf('overall                                     %f\n\n',A)
fprintf('Total C2 into   RXR train:   %f lb/h\n',m2total)
fprintf('Total C2 out of RXR train:   %f lb/h\n',m25)
fprintf('Minimum Fresh C2 to RXRs:    %f lb/h\n\n',(m2total-m25))
 
%--MASS BALANCE
stream0=mcyc;                    % total mass of stream, lb/h
stream1=m21+m61+m101+m81+mcyc;
stream2=m22+m62+m102+m82+mcyc;
stream3=m23+m63+m103+m83+mcyc;
stream4=m24+m64+m104+m84+mcyc;
stream5=m25+m65+m105+m85+mcyc;
stream0m=m20;
stream1m=m21m;
stream2m=m22m;
stream3m=m23m;
stream4m=m24m;
inletstream=stream0+stream0m+stream1m+stream2m+stream3m+stream4m;
outletstream=stream5;
%----VOLUME OF REACTOR
% use following eqn to calculate stream density:
% rho_mix = 1 / sum(x_i/rho_i) where x_i=mass fraction, rho_i=density
% ASSUMES IT IS AN IDEAL LIQUID
rho1=1/((1000/stream1)*(m21/d2+m61/d6+m101/d10+m81/d8+mcyc/dcyc)); % density of stream1, kg/m3
rho2=1/((1000/stream2)*(m22/d2+m62/d6+m102/d10+m82/d8+mcyc/dcyc));
rho3=1/((1000/stream3)*(m23/d2+m63/d6+m103/d10+m83/d8+mcyc/dcyc));
rho4=1/((1000/stream4)*(m24/d2+m64/d6+m104/d10+m84/d8+mcyc/dcyc));
rho5=1/((1000/stream5)*(m25/d2+m65/d6+m105/d10+m85/d8+mcyc/dcyc));
q1=stream1*0.453592*(1/rho1)*(1/3600); % vol flow rate, stream1, m3/s
q2=stream2*0.453592*(1/rho2)*(1/3600);
q3=stream3*0.453592*(1/rho3)*(1/3600);
q4=stream4*0.453592*(1/rho4)*(1/3600);
q5=stream5*0.453592*(1/rho5)*(1/3600);
V1=q1*x(1); % volume of rxr1 (m3), calculated using stream 1 and tau
V2=q2*x(1);
V3=q3*x(1);
V4=q4*x(1);
V5=q5*x(1);
 
fprintf('---------MASS BALANCE-----RXR VOLUME--------------\n')
fprintf('RXR   in               out               delta            RXRvolume\n')
fprintf(' 1    %f    %f     %f lb/h    %f m3\n',stream0+stream0m,stream1,-stream0-stream0m+stream1,V1)
fprintf(' 2    %f    %f     %f lb/h    %f m3\n',stream1+stream1m,stream2,-stream1-stream1m+stream2,V2)
fprintf(' 3    %f    %f     %f lb/h    %f m3\n',stream2+stream2m,stream3,-stream2-stream2m+stream3,V3)
fprintf(' 4    %f    %f     %f lb/h    %f m3\n',stream3+stream3m,stream4,-stream3-stream3m+stream4,V4)
fprintf(' 5    %f    %f     %f lb/h    %f m3\n',stream4+stream4m,stream5,-stream4-stream4m+stream5,V5)
fprintf('all   %f    %f     %f lb/h\n\n',inletstream,outletstream,-inletstream+outletstream)
fprintf('------------------OVERALL-----------------------------------\n')
fprintf('                Hexene         Decene        Octene\n')
fprintf('Yield           %f       %f      %f\n',y6,y10,y8)
fprintf('Selectivity     %f       %f      %f\n\n',s6,s10,s8)
fprintf('Tau   %f seconds\n',x(1))
fprintf('Tau   %f    hours\n\n',hours)
 
% DISTILLATION DISTILLATION DISTILLATION DISTILLATION DISTILLATION
%----Calculate bubble point of outlet stream from rxr train (stream5)-------
noct=2*q5*1000*cat; % mol/s octanol
% m2total-m25 = MIN FRESH C2 TO RXR TRAIN
 
% Iteration # 1
% ASSUMPTION = 100% recycle ethylene
mfresh2=m2total-m25; % minimum fresh C2 to rxr train, lb/h
m_eth=(mfresh2/0.999)*0.001; % mass ethane from fresh ethylene stock, lb/h
m_meth=(mfresh2/0.999)*0.0001; % mass methane from fresh ethylene stock, lb/h
n_eth=m_eth/(mweth*3600*0.00220462); % ethane, mol/s
n_meth=m_meth/(mwmeth*3600*0.00220462); % methane, mol/s
 
ztotal=n25+72+n105+n85+ncyc+noct+n_eth+n_meth;
z2=n25/(ztotal); % mol fraction in feed (any phase)
z6=72/(ztotal);
z10=n105/(ztotal);
z8=n85/(ztotal);
zcyc=ncyc/(ztotal);
zoct=noct/(ztotal);
zeth=n_eth/(ztotal);
zmeth=n_meth/(ztotal);
Pbubble=z2*p2+z6*p6+z10*p10+z8*p8+zcyc*pcyc+zoct*poct+zeth*peth+zmeth*pmeth; % bubble P of stream5
 
% - - - - - PRE-COLUMN 1 THROTTLE - - - - - - - - - 
pthrottle=37; % pressure we throttle to prior to column 1, still at rxr temp. SPECIFY
K2=p2/pthrottle;
K6=p6/pthrottle;
K10=p10/pthrottle;
K8=p8/pthrottle;
Kcyc=pcyc/pthrottle;
Koct=poct/pthrottle;
Keth=peth/pthrottle;
Kmeth=pmeth/pthrottle;
 
alphaF6=K6/Kcyc; % HK=cyclohex is reference
alphaF2=K2/Kcyc;
alphaF8=K8/Kcyc;
alphaF10=K10/Kcyc;
alphaFcyc=Kcyc/Kcyc;
alphaFoct=Koct/Kcyc;
alphaFeth=Keth/Kcyc;
alphaFmeth=Kmeth/Kcyc;
 
psi2=(z2*(1-K2))/(1+x(11)*(K2-1));
psi6=(z6*(1-K6))/(1+x(11)*(K6-1));
psi10=(z10*(1-K10))/(1+x(11)*(K10-1));
psi8=(z8*(1-K8))/(1+x(11)*(K8-1));
psicyc=(zcyc*(1-Kcyc))/(1+x(11)*(Kcyc-1));
psioct=(zoct*(1-Koct))/(1+x(11)*(Koct-1));
psieth=(zeth*(1-Keth))/(1+x(11)*(Keth-1));
psimeth=(zmeth*(1-Kmeth))/(1+x(11)*(Kmeth-1));
y(11)=psi2+psi6+psi10+psi8+psicyc+psioct+psieth+psimeth; % solve for psi, gives us VLE
 
x2=z2/(1+x(11)*(K2-1));
y2=x2*K2;
x6=z6/(1+x(11)*(K6-1));
y6=x6*K6;
V=(ztotal)*x(11);
L=(ztotal)-V;
q=1-(V/(V+L)); % q used for column 1 feed conditions
 
% - - - - - - - - - - - - -COLUMN 1 - - - - - - - - - -
% feed composition: ztotal=n25+72+n105+n85+ncyc+noct+n_eth+n_meth;
% feed temp               = rxr temp
% feed pressure           = pthrottle
 
% specify split of 2 keys     SPECIFY
B6=.03; % hexene lost in first column
D6=72-B6; % LK = hexene
Dcyc=0.2;   % HK = cyclohexane
Bcyc=ncyc-Dcyc;
 
% estimate split of nonkey components
B2e=0;                      % lights
D2e=n25-B2e;
Bethe=0;
Dethe=n_eth-Bethe;
Bmethe=0;
Dmethe=n_meth-Bmethe;
D8e=0;                      % heavies
B8e=n85-D8e;
D10e=0;
B10e=n105-D10e;
Docte=0;
Bocte=noct-Docte;
 
Dtotale=D6+Dcyc+D2e+D8e+D10e+Docte+Dethe+Dmethe;
Btotale=B6+Bcyc+B2e+B8e+B10e+Bocte+Bethe+Bmethe;
 
% estimated Distillate mol fractions
zD6e=D6/Dtotale;
zDcyce=Dcyc/Dtotale;
zD2e=D2e/Dtotale;
zD8e=D8e/Dtotale;
zD10e=D10e/Dtotale;
zDocte=Docte/Dtotale;
zDethe=Dethe/Dtotale;
zDmethe=Dmethe/Dtotale;
% estimated Bottoms mol fractions
zB6e=B6/Btotale;
zBcyce=Bcyc/Btotale;
zB2e=B2e/Btotale;
zB8e=B8e/Btotale;
zB10e=B10e/Btotale;
zBocte=Bocte/Btotale;
zBethe=Bethe/Btotale;
zBmethe=Bmethe/Btotale;
 
pB=32; % bottoms pressure      SPECIFY
pD=32; % distillate pressure   SPECIFY
 
% distillate calc
tD=x(12); % distillate temp
TD=tD+273.15;
p6D=0.019337*exp(15.8089-2654.81/(TD-47.3)); %vapor P, hexene, psi
pcycD=0.019337*(exp(15.7527-2766.63/(TD+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p2D=(-0.0573)*(tD^2)+21.523*tD+461.88; %vapor P,ethylene,psi (depends on t)
p8D=0.019337*exp(15.963-3116.52/(TD-60.39)); %vapor P, octene
p10D=0.019337*exp(16.0129-3448.18/(TD-76.09)); %vapor P, decene
poctD=0.019337*exp(15.7428-3017.81/(TD-137.1));%vapor P, octanol, psi 
pethD=0.019337*exp(15.6637-1511.42/(TD-17.16)); %vapor P, ethane 
pmethD=0.019337*exp(15.2243-597.84/(TD-7.16)); %vapor P, methane 
 
K2D=p2D/pD;
K6D=p6D/pD;
K10D=p10D/pD;
K8D=p8D/pD;
KcycD=pcycD/pD;
KoctD=poctD/pD;
KethD=pethD/pD;
KmethD=pmethD/pD;
y(12)=1-(zD2e*K2D+zD6e*K6D+zD10e*K10D+zD8e*K8D+zDcyce*KcycD+zDocte*KoctD+zDethe*KethD+zDmethe*KmethD); %calculate bubble point temp, distillation, degC, DISTILLATE
alphaD6=K6D/KcycD; % HK=cyclohex is reference
alphaD2=K2D/KcycD;
alphaD8=K8D/KcycD;
alphaD10=K10D/KcycD;
alphaDoct=KoctD/KcycD;
alphaDeth=KethD/KcycD;
alphaDmeth=KmethD/KcycD;
 
% bottoms calc
tB=x(13); % bottoms temp
TB=tB+273.15;
p2B=(-0.0573)*(tB^2)+21.523*tB+461.88; %vapor P,ethylene,psi (depends on t)
p6B=0.019337*exp(15.8089-2654.81/(TB-47.3)); %vapor P, hexene, psi
pcycB=0.019337*(exp(15.7527-2766.63/(TB+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8B=0.019337*exp(15.963-3116.52/(TB-60.39)); %vapor P, octene
p10B=0.019337*exp(16.0129-3448.18/(TB-76.09)); %vapor P, decene
poctB=0.019337*exp(15.7428-3017.81/(TB-137.1));%vapor P, octanol, psi
pethB=0.019337*exp(15.6637-1511.42/(TB-17.16)); %vapor P, ethane 
pmethB=0.019337*exp(15.2243-597.84/(TB-7.16)); %vapor P, methane 
 
K2B=p2B/pB;
K6B=p6B/pB;
K10B=p10B/pB;
K8B=p8B/pB;
KcycB=pcycB/pB;
KoctB=poctB/pB;
KethB=pethB/pB;
KmethB=pmethB/pB;
 
y(13)=1-(zB2e*K2B+zB6e*K6B+zB10e*K10B+zB8e*K8B+zBcyce*KcycB+zBocte*KoctB+zBethe*KethB+zBmethe*KmethB); %calculate bubble point temp, distillation, degC, BOTTOMS
alphaB6=K6B/KcycB; % HK = cyclohex is reference
alphaB2=K2B/KcycB;
alphaB8=K8B/KcycB;
alphaB10=K10B/KcycB;
alphaBoct=KoctB/KcycB;
alphaBeth=KethB/KcycB;
alphaBmeth=KmethB/KcycB;
 
% mean relative volatilities:
alphaAVG6=(alphaD6*alphaB6)^.5; % HK = cyclohex is reference
alphaAVG2=(alphaD2*alphaB2)^.5;
alphaAVG8=(alphaD8*alphaB8)^.5;
alphaAVG10=(alphaD10*alphaB10)^.5;
alphaAVGoct=(alphaDoct*alphaBoct)^.5;
alphaAVGeth=(alphaDeth*alphaBeth)^.5;
alphaAVGmeth=(alphaDmeth*alphaBmeth)^.5;
alphaAVGcyc=1;
Nmin=log((D6/Dcyc)*(Bcyc/B6))/log(alphaAVG6); % number of stages for first distillation, LK=hexene, HK=cyclohex
 
% nonkey component calculations: calculate using smaller of bi and di
%lights
B2=n25/(1+(Dcyc/Bcyc)*((alphaAVG2)^Nmin)); % mol/s ethylene in bottoms
D2=n25-B2;                                 %                   distillate 
Beth=n_eth/(1+(Dcyc/Bcyc)*((alphaAVGeth)^Nmin));
Deth=n_eth-Beth; 
Bmeth=n_meth/(1+(Dcyc/Bcyc)*((alphaAVGmeth)^Nmin));
Dmeth=n_meth-Bmeth; 
%heavies
D8=(n85*(Dcyc/Bcyc)*((alphaAVG8)^Nmin))/(1+(Dcyc/Bcyc)*((alphaAVG8)^Nmin));
B8=n85-D8;
D10=(n105*(Dcyc/Bcyc)*((alphaAVG10)^Nmin))/(1+(Dcyc/Bcyc)*((alphaAVG10)^Nmin));
B10=n105-D10;
Doct=(noct*(Dcyc/Bcyc)*((alphaAVGoct)^Nmin))/(1+(Dcyc/Bcyc)*((alphaAVGoct)^Nmin));
Boct=noct-Doct;
%total
Dtotal=D2+D6+D8+D10+Doct+Dcyc+Deth+Dmeth; % molar flowrate of distillate
Btotal=B2+B6+B8+B10+Boct+Bcyc+Beth+Bmeth; %                   bottoms
 
% actual distillate mol fractions
zD6=D6/Dtotal;      % mol fraction hexene in distillate
zDcyc=Dcyc/Dtotal;
zD2=D2/Dtotal;
zD8=D8/Dtotal;
zD10=D10/Dtotal;
zDoct=Doct/Dtotal;
zDeth=Deth/Dtotal;
zDmeth=Dmeth/Dtotal;
% actual bottoms mol fractions
zB6=B6/Btotal;
zBcyc=Bcyc/Btotal;
zB2=B2/Btotal;
zB8=B8/Btotal;
zB10=B10/Btotal;
zBoct=Boct/Btotal;
zBeth=Beth/Btotal;
zBmeth=Bmeth/Btotal;
 
% compare nonkey estimate to calculated
nonkeyerror=abs(B2e-B2)+abs(D2e-D2)+abs(D8e-D8)+abs(B8e-B8)+abs(D10e-D10)+abs(B10e-B10)+abs(Docte-Doct)+abs(Bocte-Boct)+abs(Bethe-Beth)+abs(Bmethe-Bmeth);
 
% UNDERWOOD EQNS
theta2=(alphaF2*z2)/(alphaF2-x(14));
theta6=(alphaF6*z6)/(alphaF6-x(14));
theta8=(alphaF8*z8)/(alphaF8-x(14));
theta10=(alphaF10*z10)/(alphaF10-x(14));
thetaoct=(alphaFoct*zoct)/(alphaFoct-x(14));
thetaeth=(alphaFeth*zeth)/(alphaFeth-x(14));
thetameth=(alphaFmeth*zmeth)/(alphaFmeth-x(14));
thetacyc=(alphaFcyc*zcyc)/(alphaFcyc-x(14));
y(14)=(1-q)-(theta6+theta8+theta10+thetaoct+thetacyc+thetaeth+thetameth); % VERY sensitive to initial conditions
Theta2=(alphaF2*(D2/Dtotal))/(alphaF2-x(14));
Theta6=(alphaF6*(D6/Dtotal))/(alphaF6-x(14));
Theta8=(alphaF8*(D8/Dtotal))/(alphaF8-x(14));
Theta10=(alphaF10*(D10/Dtotal))/(alphaF10-x(14));
Thetaoct=(alphaFoct*(Doct/Dtotal))/(alphaFoct-x(14));
Thetacyc=(alphaFcyc*(Dcyc/Dtotal))/(alphaFcyc-x(14));
Thetaeth=(alphaFeth*(Deth/Dtotal))/(alphaFeth-x(14));
Thetameth=(alphaFmeth*(Dmeth/Dtotal))/(alphaFmeth-x(14));
 
Rmin=(Theta6+Theta8+Theta10+Thetaoct+Thetacyc+Thetaeth+Thetameth)-1;
 
% GILLILAND EQNS
R=1.5*Rmin; % Assume optimal R/Rmin ratio     SPECIFY
XX=(R-Rmin)/(R+1);
YY=1-exp(((1+54.4*XX)/(11+117.2*XX))*((XX-1)/(XX^0.5)));
Nactual=(-Nmin-YY)/(YY-1);
 
% KIRKBRIDE EQN
stageratio=((zcyc/z6)*(((B6/Btotal)/(Dcyc/Dtotal))^2)*(Btotal/Dtotal))^0.206;
Ns=Nactual/(stageratio+1); % # stages in stripping section (below feed stage)
Nr=stageratio*Ns;          % # stages in rectifying section (above feed stage)
 
 
fprintf('----------------------------------DISTILLATION COLUMN 1\n')
fprintf('Column 1     Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',Dmeth,Bmeth)
fprintf('Ethylene     %f   %f  mol/s\n',D2,B2)
fprintf('Ethane       %f     %f\n',Deth,Beth)
fprintf('Hexene       %f    %f\n',D6,B6)
fprintf('Cyclohexane  %f     %f\n',Dcyc,Bcyc)
fprintf('Octene       %f     %f\n',D8,B8)
fprintf('Decene       %f     %f\n',D10,B10)
fprintf('Octanol      %f     %f\n',Doct,Boct)
fprintf('Total        %f   %f\n\n',Dtotal,Btotal)
fprintf('Temp         %f   %f degC\n',tD,tB)
fprintf('Pressure     %f   %f psi\n\n',pD,pB)
fprintf('N  = %f\n',Nactual)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr)
fprintf('R = %f\n\n',R)
fprintf('Nonkey Error = %f\n\n',nonkeyerror)
 
 
 
%------------------- COLUMN 2B -----------------------------------------
% feed = light ends from C1 = Dtotal=D2+D6+D8+D10+Doct+Dcyc+Deth+Dmeth;
% majority = methane, ethylene, ethane, hexene (some cyclohexane)
% temp = tD
% pressure = pD
 
% determine q of feed
pflash=pD; % pressure of stream
K2Dflash=p2D/pflash;
K6Dflash=p6D/pflash;
K10Dflash=p10D/pflash;
K8Dflash=p8D/pflash;
KcycDflash=pcycD/pflash;
KoctDflash=poctD/pflash;
KethDflash=pethD/pflash;
KmethDflash=pmethD/pflash;
 
alphaF2Bmeth=KmethDflash/K6Dflash;
alphaF2B2=K2Dflash/K6Dflash;
alphaF2Beth=KethDflash/K6Dflash; % LK=ethane
alphaF2B6=K6Dflash/K6Dflash; % HK=hexene is reference
alphaF2Bcyc=KcycDflash/K6Dflash;
alphaF2B10=K10Dflash/K6Dflash; 
alphaF2B8=K8Dflash/K6Dflash;  
alphaF2Boct=KoctDflash/K6Dflash;
 
psi2D=(zD2*(1-K2Dflash))/(1+x(15)*(K2Dflash-1));
psi6D=(zD6*(1-K6Dflash))/(1+x(15)*(K6Dflash-1));
psi10D=(zD10*(1-K10Dflash))/(1+x(15)*(K10Dflash-1));
psi8D=(zD8*(1-K8Dflash))/(1+x(15)*(K8Dflash-1));
psicycD=(zDcyc*(1-KcycDflash))/(1+x(15)*(KcycDflash-1));
psioctD=(zDoct*(1-KoctDflash))/(1+x(15)*(KoctDflash-1));
psimethD=(zDmeth*(1-KmethDflash))/(1+x(15)*(KmethDflash-1));
psiethD=(zDeth*(1-KethDflash))/(1+x(15)*(KethDflash-1));
y(15)=psi2D+psi6D+psi10D+psi8D+psicycD+psioctD+psimethD+psiethD; % solve for psi, gives us VLE
xD2=zD2/(1+x(15)*(K2Dflash-1)) ;% mol fraction hexene in liq phase
yD2=xD2*K2Dflash;
xD6=zD6/(1+x(15)*(K6Dflash-1));
yD6=xD6*K6Dflash;
xDcyc=zDcyc/(1+x(15)*(KcycDflash-1));
yDcyc=xDcyc*KcycDflash;
VD=(Dtotal)*x(15);
LD=(Dtotal)-VD;
qD=1-(VD/(VD+LD));
%mxD6=xD6*LD*mw6*3600*0.00220462;
%myD6=yD6*VD*mw6*3600*0.00220462;
%mxDcyc=xDcyc*LD*mwcyc*3600*0.00220462;
%myDcyc=yDcyc*VD*mwcyc*3600*0.00220462;
%mxD2=xD2*LD*mw2*3600*0.00220462;
%myD2=yD2*VD*mw2*3600*0.00220462;
%hexenepurity=mxD6/(mxD6+mxDcyc+mxD2);
%fprintf('        Liq            Vapor\n')
%fprintf('Hexene  %f     %f\n',xD6*LD*mw6*3600*0.00220462,yD6*VD*mw6*3600*0.00220462)
 
% specify split of 2 keys (LK=ethane, HK=hexene)    SPECIFY
B2eth=.00001;           %LK=ethane
D2eth=Deth-B2eth;
D26=0.000001;            %HK=hexene         
B26=D6-D26;
 
%estimatedhexenepurity=(D6*mw6)/(D6*mw6+Dcyc*mwcyc+B2eth*mweth+B22*mw2) % ONLYcontaminated with ethylene,ethane,cyclohex
 
% estimate split of nonkeys
Bmethe2=0;                  % light
Dmethe2=Dmeth-Bmethe2;
B2e2=0.3;
D2e2=D2-B2e2;
Dcyce2=0;                   % heavies
Bcyce2=Dcyc-Dcyce2;
D8e2=0;                     
B8e2=D8-D8e2;
D10e2=0;
B10e2=D10-D10e2;
Docte2=0;
Bocte2=Doct-Docte2;
 
Dtotale2=D2e2+D2eth+Dmethe2+D26+Dcyce2+D8e2+D10e2+Docte2;
Btotale2=B2e2+B2eth+Bmethe2+B26+Bcyce2+B8e2+B10e2+Bocte2;
 
% estimated Distillate mol fractions
zD26=D26/Dtotale2;       %HK
zDcyce2=Dcyce2/Dtotale2;
zD2e2=D2e2/Dtotale2; 
zD8e2=D8e2/Dtotale2;
zD10e2=D10e2/Dtotale2;
zDocte2=Docte2/Dtotale2;
zDethe2=D2eth/Dtotale2; %LK
zDmethe2=Dmethe2/Dtotale2;
% estimated Bottoms mol fractions
zB6e2=B26/Btotale2;
zBcyce2=Bcyce2/Btotale2;
zB2e2=B2e2/Btotale2; 
zB8e2=B8e2/Btotale2;
zB10e2=B10e2/Btotale2;
zBocte2=Bocte2/Btotale2;
zBethe2=B2eth/Btotale2; 
zBmethe2=Bmethe2/Btotale2;
 
pB2=27; % bottoms pressure 2            SPECIFY
pD2=27; % distillate pressure 2         SPECIFY
 
% distillate calc
tD2=x(16); % distillate temp
TD2=tD2+273.15;
pmethD2=0.019337*exp(15.2243-597.84/(TD2-7.16)); %vapor P, methane 
p2D2=(-0.0573)*(tD2^2)+21.523*tD2+461.88; %vapor P,ethylene,psi (depends on t)
pethD2=0.019337*exp(15.6637-1511.42/(TD2-17.16)); %vapor P, ethane 
p6D2=0.019337*exp(15.8089-2654.81/(TD2-47.3)); %vapor P, hexene, psi
pcycD2=0.019337*(exp(15.7527-2766.63/(TD2+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8D2=0.019337*exp(15.963-3116.52/(TD2-60.39)); %vapor P, octene
p10D2=0.019337*exp(16.0129-3448.18/(TD2-76.09)); %vapor P, decene
poctD2=0.019337*exp(15.7428-3017.81/(TD2-137.1));%vapor P, octanol, psi 
K2D2=p2D2/pD2;
K6D2=p6D2/pD2;
K10D2=p10D2/pD2;
K8D2=p8D2/pD2;
KcycD2=pcycD2/pD2;
KoctD2=poctD2/pD2;
KethD2=pethD2/pD2;
KmethD2=pmethD2/pD2;
y(16)=1-(zD2e2/K2D2+zD26/K6D2+zD10e2/K10D2+zD8e2/K8D2+zDcyce2/KcycD2+zDocte2/KoctD2+zDethe2/KethD2+zDmethe2/KmethD2); %calculate dew point temp, PARTIAL CONDENSER, degC, DISTILLATE
alphaD26=K6D2/K6D2; % HK=hexene is reference
alphaD22=K2D2/K6D2;
alphaD28=K8D2/K6D2;
alphaD210=K10D2/K6D2;
alphaD2oct=KoctD2/K6D2;
alphaD2eth=KethD2/K6D2;
alphaD2meth=KmethD2/K6D2;
alphaD2cyc=KcycD2/K6D2;
 
% bottoms calc
tB2=x(17); % bottoms temp
TB2=tB2+273.15;
pmethB2=0.019337*exp(15.2243-597.84/(TB2-7.16)); %vapor P, methane
p2B2=(-0.0573)*(tB2^2)+21.523*tB2+461.88; %vapor P,ethylene,psi (depends on t)
pethB2=0.019337*exp(15.6637-1511.42/(TB2-17.16)); %vapor P, ethane
p6B2=0.019337*exp(15.8089-2654.81/(TB2-47.3)); %vapor P, hexene, psi
pcycB2=0.019337*(exp(15.7527-2766.63/(TB2+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8B2=0.019337*exp(15.963-3116.52/(TB2-60.39)); %vapor P, octene
p10B2=0.019337*exp(16.0129-3448.18/(TB2-76.09)); %vapor P, decene
poctB2=0.019337*exp(15.7428-3017.81/(TB2-137.1));%vapor P, octanol, psi 
K2B2=p2B2/pB2;
K6B2=p6B2/pB2;
K10B2=p10B2/pB2;
K8B2=p8B2/pB2;
KcycB2=pcycB2/pB2;
KoctB2=poctB2/pB2;
KethB2=pethB2/pB2;
KmethB2=pmethB2/pB2;
y(17)=1-(zB2e2*K2B2+zB6e2*K6B2+zB10e2*K10B2+zB8e2*K8B2+zBcyce2*KcycB2+zBocte2*KoctB2+zBethe2*KethB2+zBmethe2*KmethB2); %calculate bubble point temp, distillation, degC, BOTTOMS
alphaB26=K6B2/K6B2; % HK=hexene is reference
alphaB22=K2B2/K6B2;
alphaB28=K8B2/K6B2;
alphaB210=K10B2/K6B2;
alphaB2oct=KoctB2/K6B2;
alphaB2eth=KethB2/K6B2;
alphaB2meth=KmethB2/K6B2;
alphaB2cyc=KcycB2/K6B2;
 
% mean relative volatilities:
alphaAVG26=(alphaD26*alphaB26)^.5; % HK = hexene is reference, LK = ethane
alphaAVG22=(alphaD22*alphaB22)^.5;
alphaAVG28=(alphaD28*alphaB28)^.5;
alphaAVG210=(alphaD210*alphaB210)^.5;
alphaAVG2oct=(alphaD2oct*alphaB2oct)^.5;
alphaAVG2eth=(alphaD2eth*alphaB2eth)^.5;
alphaAVG2meth=(alphaD2meth*alphaB2meth)^.5;
alphaAVG2cyc=(alphaD2cyc*alphaB2cyc)^.5;
 
Nmin2=log((D2eth/D26)*(B26/B2eth))/log(alphaAVG2eth);
 
% nonkey component calculations: calculate using smaller of bi and di
%lights
B2meth=Dmeth/(1+(D26/B26)*((alphaAVG2meth)^Nmin2));
D2meth=Dmeth-B2meth;
B22=D2/(1+(D26/B26)*((alphaAVG22)^Nmin2)); % mol/s ethylene in bottoms HK=hexene is reference
D22=D2-B22;
%heavies
D2cyc=(Dcyc*(D26/B26)*((alphaAVG2cyc)^Nmin2))/(1+(D26/B26)*((alphaAVG2cyc)^Nmin2));
B2cyc=Dcyc-D2cyc;
D28=(D8*(D26/B26)*((alphaAVG28)^Nmin2))/(1+(D26/B26)*((alphaAVG28)^Nmin2));
B28=D8-D28;
D210=(D10*(D26/B26)*((alphaAVG210)^Nmin2))/(1+(D26/B26)*((alphaAVG210)^Nmin2));
B210=D10-D210;
D2oct=(Doct*(D26/B26)*((alphaAVG2oct)^Nmin2))/(1+(D26/B26)*((alphaAVG2oct)^Nmin2));
B2oct=Doct-D2oct;
 
D2total=D2meth+D22+D2eth+D26+D2cyc+D210+D28+D2oct;
B2total=B2meth+B22+B2eth+B26+B2cyc+B210+B28+B2oct;
estimatedhexenepurity=(B26*mw6)/(B2meth*mwmeth+B22*mw2+B2eth*mweth+B26*mw6+B2cyc*mwcyc+B28*mw8+B210*mw10+B2oct*mwoct);
 
% actual distillate mol fractions
zD26=D26/D2total;
zDcyc2=D2cyc/D2total;
zD22=D22/D2total; 
zD82=D28/D2total;
zD102=D210/D2total;
zDoct2=D2oct/D2total;
zDeth2=D2eth/D2total; 
zDmeth2=D2meth/D2total;
 
% Underwood
theta2B10   =   (alphaF2B10*zD10)/(alphaF2B10-x(26));
theta2B2    =   (alphaF2B2*zD2)/(alphaF2B2-x(26));
theta2B6    =   (alphaF2B6*zD6)/(alphaF2B6-x(26));
theta2B8    =   (alphaF2B8*zD8)/(alphaF2B8-x(26));
theta2Bcyc  =   (alphaF2Bcyc*zDcyc)/(alphaF2Bcyc-x(26));
theta2Beth  =   (alphaF2Beth*zDeth)/(alphaF2Beth-x(26));
theta2Bmeth =   (alphaF2Bmeth*zDmeth)/(alphaF2Bmeth-x(26));
theta2Boct  =   (alphaF2Boct*zDoct)/(alphaF2Boct-x(26));
y(26)=(1-qD)-(theta2B10+theta2B2+theta2B6+theta2B8+theta2Bcyc+theta2Beth+theta2Bmeth+theta2Boct);
Theta2B10=(alphaF2B10*(D210/D2total))/(alphaF2B10-x(26));
Theta2B2=(alphaF2B2*(D22/D2total))/(alphaF2B2-x(26));
Theta2B6=(alphaF2B6*(D26/D2total))/(alphaF2B6-x(26));
Theta2B8=(alphaF2B8*(D28/D2total))/(alphaF2B8-x(26));
Theta2Bcyc=(alphaF2Bcyc*(D2cyc/D2total))/(alphaF2Bcyc-x(26));
Theta2Beth=(alphaF2Beth*(D2eth/D2total))/(alphaF2Beth-x(26));
Theta2Bmeth=(alphaF2Bmeth*(D2meth/D2total))/(alphaF2Bmeth-x(26));
Theta2Boct=(alphaF2Boct*(D2oct/D2total))/(alphaF2Boct-x(26));
Rmin2B=(Theta2B10+Theta2B2+Theta2B6+Theta2B8+Theta2Bcyc+Theta2Beth+Theta2Bmeth+Theta2Boct)-1;
 
% GILLILAND EQNS
R2B=1.5*Rmin2B; % Assume optimal R/Rmin ratio       SPECIFY
XX2B=(R2B-Rmin2B)/(R2B+1);
YY2B=1-exp(((1+54.4*XX2B)/(11+117.2*XX2B))*((XX2B-1)/(XX2B^0.5)));
Nactual2B=(-Nmin2-YY2B)/(YY2B-1);
 
% KIRKBRIDE EQN
stageratio2B=((zD6/zDeth)*(((B2eth/B2total)/(D26/D2total))^2)*(B2total/D2total))^0.206;
Ns2B=Nactual2B/(stageratio2B+1); % # stages in stripping section (below feed stage)
Nr2B=stageratio2B*Ns2B;          % # stages in rectifying section (above feed stage)
 
 
fprintf('----------------------------------DISTILLATION COLUMN 2B\n')
fprintf('Column 2B    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D2meth,B2meth)
fprintf('Ethylene     %f   %f  mol/s\n',D22,B22)
fprintf('Ethane       %f     %f\n',D2eth,B2eth)
fprintf('Hexene       %f    %f\n',D26,B26)
fprintf('Cyclohexane  %f     %f\n',D2cyc,B2cyc)
fprintf('Octene       %f     %f\n',D28,B28)
fprintf('Decene       %f     %f\n',D210,B210)
fprintf('Octanol      %f     %f\n',D2oct,B2oct)
fprintf('Total        %f   %f\n\n',D2total,B2total)
fprintf('Temp         %f   %f degC\n',tD2,tB2)
fprintf('Pressure     %f   %f psi\n\n',pD2,pB2)
fprintf('N  = %f\n',Nactual2B)
fprintf('Nmin  = %f\n',Nmin2)
fprintf('Hexene Purity = %f\n',estimatedhexenepurity)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns2B)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr2B)
fprintf('R = %f\n\n',R2B)
 
 
 
%------------------- COLUMN 3B-----------------------------------------
% Condenser and knockout drums prior to inlet!
% feed- compressed lights from C2B: C2total=C2meth+C22+C2eth+C26+C2cyc+C210+C28+C2oct;
% majority = methane, ethylene, ethane (some hexene)
% temp = compressor outlet temp, specify below (TCmp or tCmp)
% pressure = compressor outlet p, specify below (pCmp)
 
%COMPRESSOR CONDITIONS 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
pCmp=400;           % outlet P of compressors, psi  SPECIFY
tCmp=15.3;        % outlet T of compressors, degC    SPECIFY
TCmp=tCmp+273.15;   % Kelvin
 
C2meth=D2meth;      % composition of outlet stream, mol/s  SPECIFY
C22   =D22;
C2eth =D2eth;
C26   =D26;
C2cyc =D2cyc;
C210  =D210;
C28   =D28;
C2oct =D2oct;
 
C2total=C2meth+C22+C2eth+C26+C2cyc+C210+C28+C2oct;
 
zC26=C26/C2total     ;      % mol fractions, feed
zCcyc2=C2cyc/C2total;
zC22=C22/C2total;
zC82=C28/C2total;
zC102=C210/C2total;
zCoct2=C2oct/C2total;
zCeth2=C2eth/C2total;
zCmeth2=C2meth/C2total;
 
%0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 
% 3B determine q of feed
pmethCmp=0.019337*exp(15.2243-597.84/(TCmp-7.16)); %vapor P, methane (TEMP RANGE!)
p2Cmp=(-0.0573)*(tCmp^2)+21.523*tCmp+461.88; %vapor P,ethylene,psi (depends on t)
%p2Cmp=0.019337*exp(15.5368-1347.01/(TCmp-18.15)); %KNOVEL VAPOR P!
pethCmp=0.019337*exp(15.6637-1511.42/(TCmp-17.16)); %vapor P, ethane (TEMP RANGE!)
p6Cmp=0.019337*exp(15.8089-2654.81/(TCmp-47.3)); %vapor P, hexene, psi
pcycCmp=0.019337*(exp(15.7527-2766.63/(TCmp+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8Cmp=0.019337*exp(15.963-3116.52/(TCmp-60.39)); %vapor P, octene
p10Cmp=0.019337*exp(16.0129-3448.18/(TCmp-76.09)); %vapor P, decene
poctCmp=0.019337*exp(15.7428-3017.81/(TCmp-137.1));%vapor P, octanol, psi 
 
pflash3B=pCmp; % pressure of stream
K2C2flash=p2Cmp/pflash3B;
K6C2flash=p6Cmp/pflash3B;
K10C2flash=p10Cmp/pflash3B;
K8C2flash=p8Cmp/pflash3B;
KcycC2flash=pcycCmp/pflash3B;
KoctC2flash=poctCmp/pflash3B;
KethC2flash=pethCmp/pflash3B;
KmethC2flash=pmethCmp/pflash3B;
 
alphaF3Bmeth=KmethC2flash/KethC2flash;
alphaF3B2=K2C2flash/KethC2flash  ;    % LK=ethylene
alphaF3Beth=KethC2flash/KethC2flash;  % HK=ethane is reference
alphaF3B6=K6C2flash/KethC2flash;
alphaF3Bcyc=KcycC2flash/KethC2flash;
alphaF3B10=K10C2flash/KethC2flash; 
alphaF3B8=K8C2flash/KethC2flash;  
alphaF3Boct=KoctC2flash/KethC2flash;
psi3B2=(zC22*(1-K2C2flash))/(1+x(27)*(K2C2flash-1));
psi3B6=(zC26*(1-K6C2flash))/(1+x(27)*(K6C2flash-1));
psi3B10=(zC102*(1-K10C2flash))/(1+x(27)*(K10C2flash-1));
psi3B8=(zC82*(1-K8C2flash))/(1+x(27)*(K8C2flash-1));
psi3Bcyc=(zCcyc2*(1-KcycC2flash))/(1+x(27)*(KcycC2flash-1));
psi3Boct=(zCoct2*(1-KoctC2flash))/(1+x(27)*(KoctC2flash-1));
psi3Bmeth=(zCmeth2*(1-KethC2flash))/(1+x(27)*(KethC2flash-1));
psi3Beth=(zCeth2*(1-KmethC2flash))/(1+x(27)*(KmethC2flash-1));
y(27)=psi3B2+psi3B6+psi3B10+psi3B8+psi3Bcyc+psi3Boct+psi3Bmeth+psi3Beth; % solve for psi, gives us VLE
VF3B=(C2total)*x(27);
LF3B=(C2total)-VF3B;
qF3B=0; %1-(VF3B/(VF3B+LF3B));
 
% specify split of 2 keys (LK=ethylene, HK=ethane)  SPECIFY
B3B2=.01;         %LK
D3B2=C22-B3B2;
D3Beth=.01;       %HK
B3Beth=C2eth-D3Beth;
 
% estimate split of nonkeys
% light
B3Bemeth=0;
D3Bemeth=C2meth-B3Bemeth;
% heavy
D3Be6=0;
B3Be6=C26-D3Be6;
D3Becyc=0;
B3Becyc=C2cyc-D3Becyc;
D3Be8=0;
B3Be8=C28-D3Be8;
D3Be10=0;
B3Be10=C210-D3Be10;
D3Beoct=0;
B3Beoct=C2oct-D3Beoct;
 
B3Betotal=B3Be10+B3Be8+B3B2+B3Be6+B3Becyc+B3Beth+B3Bemeth+B3Beoct;
D3Betotal=D3Be10+D3Be8+D3B2+D3Be6+D3Becyc+D3Beth+D3Bemeth+D3Beoct;
 
% estimated Distillate mol fractions
zD3Bemeth=D3Bemeth/D3Betotal;
zD3Be2=D3B2/D3Betotal;
zD3Beeth=D3Beth/D3Betotal;
zD3Be6=D3Be6/D3Betotal;
zD3Becyc=D3Becyc/D3Betotal;
zD3Be8=D3Be8/D3Betotal;
zD3Be10=D3Be10/D3Betotal;
zD3Beoct=D3Beoct/D3Betotal;
% estimated Bottoms mol fractions
zB3Bemeth=B3Bemeth/B3Betotal;
zB3Be2=B3B2/B3Betotal;
zB3Beeth=B3Beth/B3Betotal;
zB3Be6=B3Be6/B3Betotal;
zB3Becyc=B3Becyc/B3Betotal;
zB3Be8=B3Be8/B3Betotal;
zB3Be10=B3Be10/B3Betotal;
zB3Beoct=B3Beoct/B3Betotal;
 
pB3B=395; % bottoms pressure        SPECIFY
pD3B=395; % distillate pressure     SPECIFY
% 3B distillate calc
tD3B=x(28); % distillate temp
TD3B=tD3B+273.15;
pmethD3B=0.019337*exp(15.2243-597.84/(TD3B-7.16)); %vapor P, methane
p2D3B=(-0.0573)*(tD3B^2)+21.523*tD3B+461.88; %vapor P,ethylene,psi (depends on t)
pethD3B=0.019337*exp(15.6637-1511.42/(TD3B-17.16)); %vapor P, ethane 
p6D3B=0.019337*exp(15.8089-2654.81/(TD3B-47.3)); %vapor P, hexene, psi
pcycD3B=0.019337*(exp(15.7527-2766.63/(TD3B+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8D3B=0.019337*exp(15.963-3116.52/(TD3B-60.39)); %vapor P, octene
p10D3B=0.019337*exp(16.0129-3448.18/(TD3B-76.09)); %vapor P, decene
poctD3B=0.019337*exp(15.7428-3017.81/(TD3B-137.1));%vapor P, octanol, psi 
K2D3B=p2D3B/pD3B;
K6D3B=p6D3B/pD3B;
K10D3B=p10D3B/pD3B;
K8D3B=p8D3B/pD3B;
KcycD3B=pcycD3B/pD3B;
KoctD3B=poctD3B/pD3B;
KethD3B=pethD3B/pD3B;
KmethD3B=pmethD3B/pD3B;
y(28)=1-(zD3Be10/K10D3B+zD3Be2/K2D3B+zD3Be6/K6D3B+zD3Be8/K8D3B+zD3Becyc/KcycD3B+zD3Beeth/KethD3B+zD3Bemeth/KmethD3B+zD3Beoct/KoctD3B); %calculate dew point temp, partial CONDENSER, degC, DISTILLATE
alphaD3B2=K2D3B/KethD3B; % LK=ethylene
alphaD3B6=K6D3B/KethD3B;
alphaD3B10=K10D3B/KethD3B; 
alphaD3B8=K8D3B/KethD3B;
alphaD3Bcyc=KcycD3B/KethD3B;
alphaD3Boct=KoctD3B/KethD3B;
alphaD3Beth=KethD3B/KethD3B; % HK=ethane is reference
alphaD3Bmeth=KmethD3B/KethD3B;
 
% 3B bottoms calc
tB3B=x(29); % bottoms temp
TB3B=tB3B+273.15;
pmethB3B=0.019337*exp(15.2243-597.84/(TB3B-7.16)); %vapor P, methane
p2B3B=(-0.0573)*(tB3B^2)+21.523*tB3B+461.88; %vapor P,ethylene,psi (depends on t)
pethB3B=0.019337*exp(15.6637-1511.42/(TB3B-17.16)); %vapor P, ethane 
p6B3B=0.019337*exp(15.8089-2654.81/(TB3B-47.3)); %vapor P, hexene, psi
pcycB3B=0.019337*(exp(15.7527-2766.63/(TB3B+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8B3B=0.019337*exp(15.963-3116.52/(TB3B-60.39)); %vapor P, octene
p10B3B=0.019337*exp(16.0129-3448.18/(TB3B-76.09)); %vapor P, decene
poctB3B=0.019337*exp(15.7428-3017.81/(TB3B-137.1));%vapor P, octanol, psi 
K2B3B=p2B3B/pB3B; 
K6B3B=p6B3B/pB3B;
K10B3B=p10B3B/pB3B;
K8B3B=p8B3B/pB3B;
KcycB3B=pcycB3B/pB3B;
KoctB3B=poctB3B/pB3B;
KethB3B=pethB3B/pB3B;
KmethB3B=pmethB3B/pB3B;
y(29)=1-(zB3Be10*K10B3B+zB3Be2*K2B3B+zB3Be6*K6B3B+zB3Be8*K8B3B+zB3Becyc*KcycB3B+zB3Beeth*KethB3B+zB3Bemeth*KmethB3B+zB3Beoct*KoctB3B); %calculate bubble point temp, degC, Bottoms
alphaB3B2=K2B3B/KethB3B;
alphaB3B6=K6B3B/KethB3B;
alphaB3B10=K10B3B/KethB3B; 
alphaB3B8=K8B3B/KethB3B;
alphaB3Bcyc=KcycB3B/KethB3B;
alphaB3Boct=KoctB3B/KethB3B;
alphaB3Beth=KethB3B/KethB3B; % HK=ethane is reference
alphaB3Bmeth=KmethB3B/KethB3B;
 
% mean relative volatilities:
alphaAVG3B2 =(alphaD3B2*alphaB3B2)^.5; %LK
alphaAVG3B6 =(alphaD3B6*alphaB3B6)^.5;
alphaAVG3B10=(alphaD3B10*alphaB3B10)^.5;
alphaAVG3B8 =(alphaD3B8*alphaB3B8)^.5;
alphaAVG3Bcyc=(alphaD3Bcyc*alphaB3Bcyc)^.5;
alphaAVG3Boct=(alphaD3Boct*alphaB3Boct)^.5;
alphaAVG3Beth=(alphaD3Beth*alphaB3Beth)^.5; % HK=ethane is reference
alphaAVG3Bmeth=(alphaD3Bmeth*alphaB3Bmeth)^.5;
 
Nmin3B=log((D3B2/D3Beth)*(B3Beth/B3B2))/log(alphaAVG3B2);
 
% nonkey component calculations: calculate using smaller of bi and di
%lights
B3Bmeth=C2meth/(1+(D3Beth/B3Beth)*((alphaAVG3Bmeth)^Nmin3B));
D3Bmeth=C2meth-B3Bmeth;
%heavies
D3B6=(C26*(D3Beth/B3Beth)*((alphaAVG3B6)^Nmin3B))/(1+(D3Beth/B3Beth)*((alphaAVG3B6)^Nmin3B));
B3B6=C26-D3B6;
D3Bcyc=(C2cyc*(D3Beth/B3Beth)*((alphaAVG3Bcyc)^Nmin3B))/(1+(D3Beth/B3Beth)*((alphaAVG3Bcyc)^Nmin3B));
B3Bcyc=C2cyc-D3Bcyc;
D3B8=(C28*(D3Beth/B3Beth)*((alphaAVG3B8)^Nmin3B))/(1+(D3Beth/B3Beth)*((alphaAVG3B8)^Nmin3B));
B3B8=C28-D3B8;
D3B10=(C210*(D3Beth/B3Beth)*((alphaAVG3B10)^Nmin3B))/(1+(D3Beth/B3Beth)*((alphaAVG3B10)^Nmin3B));
B3B10=C210-D3B10;
D3Boct=(C2oct*(D3Beth/B3Beth)*((alphaAVG3Boct)^Nmin3B))/(1+(D3Beth/B3Beth)*((alphaAVG3Boct)^Nmin3B));
B3Boct=C2oct-D3Boct;
 
D3Btotal=D3Bmeth+D3B2+D3Beth+D3B6+D3Bcyc+D3B10+D3B8+D3Boct;
B3Btotal=B3Bmeth+B3B2+B3Beth+B3B6+B3Bcyc+B3B10+B3B8+B3Boct;
 
% Underwood
theta3B10=(alphaAVG3B10*zC102)/(alphaAVG3B10-x(30));
theta3B2=(alphaAVG3B2*zC22)/(alphaAVG3B2-x(30));
theta3B6=(alphaAVG3B6*zC26)/(alphaAVG3B6-x(30));
theta3B8=(alphaAVG3B8*zC82)/(alphaAVG3B8-x(30));
theta3Bcyc=(alphaAVG3Bcyc*zCcyc2)/(alphaAVG3Bcyc-x(30));
theta3Beth=(alphaAVG3Beth*zCeth2)/(alphaAVG3Beth-x(30));
theta3Bmeth=(alphaAVG3Bmeth*zCmeth2)/(alphaAVG3Bmeth-x(30));
theta3Boct=(alphaAVG3Boct*zCoct2)/(alphaAVG3Boct-x(30));
y(30)=(1-qF3B)-(theta3B10+theta3B2+theta3B6+theta3B8+theta3Bcyc+theta3Beth+theta3Bmeth+theta3Boct);
 
Theta3B10=(alphaAVG3B10*(D3B10/D3Btotal))/(alphaAVG3B10-x(25));
Theta3B2=(alphaAVG3B2*(D3B2/D3Btotal))/(alphaAVG3B2-x(25));
Theta3B6=(alphaAVG3B6*(D3B6/D3Btotal))/(alphaAVG3B6-x(25));
Theta3B8=(alphaAVG3B8*(D3B8/D3Btotal))/(alphaAVG3B8-x(25));
Theta3Bcyc=(alphaAVG3Bcyc*(D3Bcyc/D3Btotal))/(alphaAVG3Bcyc-x(25));
Theta3Beth=(alphaAVG3Beth*(D3Beth/D3Btotal))/(alphaAVG3Beth-x(25));
Theta3Bmeth=(alphaAVG3Bmeth*(D3Bmeth/D3Btotal))/(alphaAVG3Bmeth-x(25));
Theta3Boct=(alphaAVG3Boct*(D3Boct/D3Btotal))/(alphaAVG3Boct-x(25));
 
Rmin3B=(Theta3B10+Theta3B2+Theta3B6+Theta3B8+Theta3Bcyc+Theta3Beth+Theta3Bmeth+Theta3Boct)-1;
 
% GILLILAND EQNS
R3B=1.5*Rmin3B; % Assume optimal R/Rmin ratio       SPECIFY
XX3B=(R3B-Rmin3B)/(R3B+1);
YY3B=1-exp(((1+54.4*XX3B)/(11+117.2*XX3B))*((XX3B-1)/(XX3B^0.5)));
Nactual3B=(-Nmin3B-YY3B)/(YY3B-1);
 
% KIRKBRIDE EQN
stageratio3B=((zCeth2/zC22)*(((B3B2/B3Btotal)/(D3Beth/D3Btotal))^2)*(B3Btotal/D3Btotal))^0.206;
Ns3B=Nactual3B/(stageratio3B+1); % # stages in stripping section (below feed stage)
Nr3B=stageratio3B*Ns3B;          % # stages in rectifying section (above feed stage)
 
fprintf('----------------------------------DISTILLATION COLUMN 3B\n')
fprintf('Column 3B    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D3Bmeth,B3Bmeth)
fprintf('Ethylene     %f   %f  mol/s\n',D3B2,B3B2)
fprintf('Ethane       %f     %f\n',D3Beth,B3Beth)
fprintf('Hexene       %f    %f\n',D3B6,B3B6)
fprintf('Cyclohexane  %f     %f\n',D3Bcyc,B3Bcyc)
fprintf('Octene       %f     %f\n',D3B8,B3B8)
fprintf('Decene       %f     %f\n',D3B10,B3B10)
fprintf('Octanol      %f     %f\n',D3Boct,B3Boct)
fprintf('Total        %f   %f\n\n',D3Btotal,B3Btotal)
fprintf('Temp         %f   %f degC\n',tD3B,tB3B)
fprintf('Pressure     %f   %f psi\n\n',pD3B,pB3B)
fprintf('N  = %f\n',Nactual3B)
fprintf('Nmin  = %f\n',Nmin3B)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns3B)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr3B)
fprintf('R = %f\n\n',R3B)
 
 
%------------------- COLUMN 4B ----------------------------------------
% feed = distillate from C3B = D3Btotal=D3Bmeth+D3B2+D3Beth+D3B6+D3Bcyc+D3B10+D3B8+D3Boct;
% majority = methane, ethylene (some ethane)
% temp = tD3B
% pressure = pD3B
 
 
% determine q of feed
pflash4B=pD3B;  % pressure of stream
K2D3flash=p2D3B/pflash4B;
K6D3flash=p6D3B/pflash4B;
K10D3flash=p10D3B/pflash4B;
K8D3flash=p8D3B/pflash4B;
KcycD3flash=pcycD3B/pflash4B;
KoctD3flash=poctD3B/pflash4B;
KethD3flash=pethD3B/pflash4B;
KmethD3flash=pmethD3B/pflash4B;
alphaF4Bmeth=KmethD3flash/K2D3flash; % LK=methane
alphaF4B2=K2D3flash/K2D3flash;      % HK=ethylene is reference
alphaF4Beth=KethD3flash/K2D3flash;
alphaF4B6=K6D3flash/K2D3flash;
alphaF4Bcyc=KcycD3flash/K2D3flash;
alphaF4B10=K10D3flash/K2D3flash;
alphaF4B8=K8D3flash/K2D3flash; 
alphaF4Boct=KoctD3flash/K2D3flash;
qF4B=0;
 
% specify split of 2 keys (LK=methane, HK=ethylene)        SPECIFY
B4Bmeth=.000001;         %LK=methane
D4Bmeth=D3Bmeth-B4Bmeth;
D4B2=.000001;            %HK=ethylene (2)
B4B2=D3B2-D4B2;
 
% estimate split of nonkeys
% heavy
D4Beeth=0;
B4Beeth=D3Beth-D4Beeth;
D4Be6=0;
B4Be6=D3B6-D4Be6;
D4Becyc=0;
B4Becyc=D3Bcyc-D4Becyc;
D4Be8=0;
B4Be8=D3B8-D4Be8;
D4Be10=0;
B4Be10=D3B10-D4Be10;
D4Beoct=0;
B4Beoct=D3Boct-D4Beoct;
 
B4Betotal=B4Be10+B4Be8+B4B2+B4Be6+B4Becyc+B4Beeth+B4Bmeth+B4Beoct;
D4Betotal=D4Be10+D4Be8+D4B2+D4Be6+D4Becyc+D4Beeth+D4Bmeth+D4Beoct;
 
% estimated Distillate mol fractions
zD4Bemeth=D4Bmeth/D4Betotal;
zD4Be2=D4B2/D4Betotal;
zD4Beeth=D4Beeth/D4Betotal;
zD4Be6=D4Be6/D4Betotal;
zD4Becyc=D4Becyc/D4Betotal;
zD4Be8=D4Be8/D4Betotal;
zD4Be10=D4Be10/D4Betotal;
zD4Beoct=D4Beoct/D4Betotal;
% estimated Bottoms mol fractions
zB4Bemeth=B4Bmeth/B4Betotal;
zB4Be2=B4B2/B4Betotal;
zB4Beeth=B4Beeth/B4Betotal;
zB4Be6=B4Be6/B4Betotal;
zB4Becyc=B4Becyc/B4Betotal;
zB4Be8=B4Be8/B4Betotal;
zB4Be10=B4Be10/B4Betotal;
zB4Beoct=B4Beoct/B4Betotal;
 
pB4B=390; % bottoms pressure         SPECIFY        
pD4B=390; % distillate pressure      SPECIFY 
 
% 4B distillate calc
tD4B=x(31); % distillate temp
TD4B=tD4B+273.15;
pmethD4B=0.019337*exp(15.2243-597.84/(TD4B-7.16)); %vapor P, methane 
%p2D4B=(-0.0573)*(tD4B^2)+21.523*tD4B+461.88; %vapor P,ethylene,psi (depends on t)
p2D4B=0.019337*exp(15.5368-1347.01/(TD4B-18.15)); %KNOVEL VAPOR P!
pethD4B=0.019337*exp(15.6637-1511.42/(TD4B-17.16)); %vapor P, ethane 
p6D4B=0.019337*exp(15.8089-2654.81/(TD4B-47.3)); %vapor P, hexene, psi
pcycD4B=0.019337*(exp(15.7527-2766.63/(TD4B+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8D4B=0.019337*exp(15.963-3116.52/(TD4B-60.39)); %vapor P, octene
p10D4B=0.019337*exp(16.0129-3448.18/(TD4B-76.09)); %vapor P, decene
poctD4B=0.019337*exp(15.7428-3017.81/(TD4B-137.1));%vapor P, octanol, psi 
K2D4B=p2D4B/pD4B;
K6D4B=p6D4B/pD4B;
K10D4B=p10D4B/pD4B;
K8D4B=p8D4B/pD4B;
KcycD4B=pcycD4B/pD4B;
KoctD4B=poctD4B/pD4B;
KethD4B=pethD4B/pD4B;
KmethD4B=pmethD4B/pD4B;
y(31)=1-(zD4Be10/K10D4B+zD4Be2/K2D4B+zD4Be6/K6D4B+zD4Be8/K8D4B+zD4Becyc/KcycD4B+zD4Beeth/KethD4B+zD4Bemeth/KmethD4B+zD4Beoct/KoctD4B); %calculate dew point temp, partial CONDENSER, degC, DISTILLATE
alphaD4B2=K2D4B/K2D4B;  %reference is HK=ethylene (2)
alphaD4B6=K6D4B/K2D4B;
alphaD4B10=K10D4B/K2D4B; 
alphaD4B8=K8D4B/K2D4B;
alphaD4Bcyc=KcycD4B/K2D4B;
alphaD4Boct=KoctD4B/K2D4B;
alphaD4Beth=KethD4B/K2D4B;
alphaD4Bmeth=KmethD4B/K2D4B;
 
% bottoms calc
tB4B=x(32); % distillate temp
TB4B=tB4B+273.15;
pmethB4B=0.019337*exp(15.2243-597.84/(TB4B-7.16)); %vapor P, methane 
p2B4B=(-0.0573)*(tB4B^2)+21.523*tB4B+461.88; %vapor P,ethylene,psi (depends on t)
pethB4B=0.019337*exp(15.6637-1511.42/(TB4B-17.16)); %vapor P, ethane 
p6B4B=0.019337*exp(15.8089-2654.81/(TB4B-47.3)); %vapor P, hexene, psi
pcycB4B=0.019337*(exp(15.7527-2766.63/(TB4B+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8B4B=0.019337*exp(15.963-3116.52/(TB4B-60.39)); %vapor P, octene
p10B4B=0.019337*exp(16.0129-3448.18/(TB4B-76.09)); %vapor P, decene
poctB4B=0.019337*exp(15.7428-3017.81/(TB4B-137.1));%vapor P, octanol, psi 
K2B4B=p2B4B/pB4B;
K6B4B=p6B4B/pB4B;
K10B4B=p10B4B/pB4B;
K8B4B=p8B4B/pB4B;
KcycB4B=pcycB4B/pB4B;
KoctB4B=poctB4B/pB4B;
KethB4B=pethB4B/pB4B;
KmethB4B=pmethB4B/pB4B;
y(32)=1-(zB4Be10/K10B4B+zB4Be2/K2B4B+zB4Be6/K6B4B+zB4Be8/K8B4B+zB4Becyc/KcycB4B+zB4Beeth/KethB4B+zB4Bemeth/KmethB4B+zB4Beoct/KoctB4B); %calculate DEW point temp, PARTIAL REBOILER distillation, degC, BOTTOMS
alphaB4B2=K2B4B/K2B4B;  %reference is HK=ethylene (2)
alphaB4B6=K6B4B/K2B4B;
alphaB4B10=K10B4B/K2B4B;
alphaB4B8=K8B4B/K2B4B;
alphaB4Bcyc=KcycB4B/K2B4B;
alphaB4Boct=KoctB4B/K2B4B;
alphaB4Beth=KethB4B/K2B4B;
alphaB4Bmeth=KmethB4B/K2B4B;
 
% mean relative volatilities:
alphaAVG4B2 =(alphaD4B2*alphaB4B2)^.5; %reference is HK=ethylene (2)
alphaAVG4B6 =(alphaD4B6*alphaB4B6)^.5;
alphaAVG4B10=(alphaD4B10*alphaB4B10)^.5;
alphaAVG4B8 =(alphaD4B8*alphaB4B8)^.5;
alphaAVG4Bcyc=(alphaD4Bcyc*alphaB4Bcyc)^.5;
alphaAVG4Boct=(alphaD4Boct*alphaB4Boct)^.5;
alphaAVG4Beth=(alphaD4Beth*alphaB4Beth)^.5;
alphaAVG4Bmeth=(alphaD4Bmeth*alphaB4Bmeth)^.5;
 
Nmin4B=log((D4Bmeth/D4B2)*(B4B2/B4Bmeth))/log(alphaAVG4Bmeth);
 
% nonkey component calculations: calculate using smaller of bi and di
%heavies
D4Beth=(D3Beth*(D4B2/B4B2)*((alphaAVG4Beth)^Nmin4B))/(1+(D4B2/B4B2)*((alphaAVG4Beth)^Nmin4B));
B4Beth=D3Beth-D4Beth;
D4B6=(D3B6*(D4B2/B4B2)*((alphaAVG4B6)^Nmin4B))/(1+(D4B2/B4B2)*((alphaAVG4B6)^Nmin4B));
B4B6=D3B6-D4B6;
D4Bcyc=(D3Bcyc*(D4B2/B4B2)*((alphaAVG4Bcyc)^Nmin4B))/(1+(D4B2/B4B2)*((alphaAVG4Bcyc)^Nmin4B));
B4Bcyc=D3Bcyc-D4Bcyc;
D4B8=(D3B8*(D4B2/B4B2)*((alphaAVG4B8)^Nmin4B))/(1+(D4B2/B4B2)*((alphaAVG4B8)^Nmin4B));
B4B8=D3B8-D4B8;
D4B10=(D3B10*(D4B2/B4B2)*((alphaAVG4B10)^Nmin4B))/(1+(D4B2/B4B2)*((alphaAVG4B10)^Nmin4B));
B4B10=D3B10-D4B10;
D4Boct=(D3Boct*(D4B2/B4B2)*((alphaAVG4Boct)^Nmin4B))/(1+(D4B2/B4B2)*((alphaAVG4Boct)^Nmin4B));
B4Boct=D3Boct-D4Boct;
 
D4Btotal=D4Bmeth+D4B2+D4Beth+D4B6+D4Bcyc+D4B10+D4B8+D4Boct;
B4Btotal=B4Bmeth+B4B2+B4Beth+B4B6+B4Bcyc+B4B10+B4B8+B4Boct;
 
 
% Underwood
z3B10   =(D3B10/D3Btotal);
z3B2    =(D3B2/D3Btotal);
z3B6    =(D3B6/D3Btotal);
z3B8    =(D3B8/D3Btotal);
z3Bcyc  =(D3Bcyc/D3Btotal);
z3Beth  =(D3Beth/D3Btotal);
z3Bmeth =(D3Bmeth/D3Btotal);
z3Boct  =(D3Boct/D3Btotal);
theta4B10=(alphaAVG4B10*z3B10)/(alphaAVG4B10-x(33));
theta4B2=(alphaAVG4B2*z3B2)/(alphaAVG4B2-x(33));
theta4B6=(alphaAVG4B6*z3B6)/(alphaAVG4B6-x(33));
theta4B8=(alphaAVG4B8*z3B8)/(alphaAVG4B8-x(33));
theta4Bcyc=(alphaAVG4Bcyc*z3Bcyc)/(alphaAVG4Bcyc-x(33));
theta4Beth=(alphaAVG4Beth*z3Beth)/(alphaAVG4Beth-x(33));
theta4Bmeth=(alphaAVG4Bmeth*z3Bmeth)/(alphaAVG4Bmeth-x(33));
theta4Boct=(alphaAVG4Boct*z3Boct)/(alphaAVG4Boct-x(33));
y(33)=(1-qF4B)-(theta4B10+theta4B2+theta4B6+theta4B8+theta4Bcyc+theta4Beth+theta4Bmeth+theta4Boct);
Theta4B10=(alphaAVG4B10*(D4B10/D4Btotal))/(alphaAVG4B10-x(33));
Theta4B2=(alphaAVG4B2*(D4B2/D4Btotal))/(alphaAVG4B2-x(33));
Theta4B6=(alphaAVG4B6*(D4B6/D4Btotal))/(alphaAVG4B6-x(33));
Theta4B8=(alphaAVG4B8*(D4B8/D4Btotal))/(alphaAVG4B8-x(33));
Theta4Bcyc=(alphaAVG4Bcyc*(D4Bcyc/D4Btotal))/(alphaAVG4Bcyc-x(33));
Theta4Beth=(alphaAVG4Beth*(D4Beth/D4Btotal))/(alphaAVG4Beth-x(33));
Theta4Bmeth=(alphaAVG4Bmeth*(D4Bmeth/D4Btotal))/(alphaAVG4Bmeth-x(33));
Theta4Boct=(alphaAVG4Boct*(D4Boct/D4Btotal))/(alphaAVG4Boct-x(33));
 
Rmin4B=(Theta4B10+Theta4B2+Theta4B6+Theta4B8+Theta4Bcyc+Theta4Beth+Theta4Bmeth+Theta4Boct)-1;
 
% GILLILAND EQNS
R4B=1.3*Rmin4B; % Assume optimal R/Rmin ratio       SPECIFY
XX4B=(R4B-Rmin4B)/(R4B+1);
YY4B=1-exp(((1+54.4*XX4B)/(11+117.2*XX4B))*((XX4B-1)/(XX4B^0.5)));
Nactual4B=(-Nmin4B-YY4B)/(YY4B-1);
 
% KIRKBRIDE EQN
stageratio4B=((z3B2/z3Bmeth)*(((B4Bmeth/B4Btotal)/(D4B2/D4Btotal))^2)*(B4Btotal/D4Btotal))^0.206;
Ns4B=Nactual4B/(stageratio4B+1); % # stages in stripping section (below feed stage)
Nr4B=stageratio4B*Ns4B;          % # stages in rectifying section (above feed stage)
 
%D3Btotal=D3Bmeth+D3B2+D3Beth+D3B6+D3Bcyc+D3B10+D3B8+D3Boct;
% specify split of 2 keys (LK=methane, HK=ethylene)        SPECIFY
%B4Bmeth=.001;         %LK=methane
%D4Bmeth=D3Bmeth-B4Bmeth;
%D4B2=.001;            %HK=ethylene (2)
%B4B2=D3B2-D4B2;
 
 
fprintf('----------------------------------DISTILLATION COLUMN 4B\n')
fprintf('Column 4B    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D4Bmeth,B4Bmeth)
fprintf('Ethylene     %f   %f  mol/s\n',D4B2,B4B2)
fprintf('Ethane       %f     %f\n',D4Beth,B4Beth)
fprintf('Hexene       %f    %f\n',D4B6,B4B6)
fprintf('Cyclohexane  %f     %f\n',D4Bcyc,B4Bcyc)
fprintf('Octene       %f     %f\n',D4B8,B4B8)
fprintf('Decene       %f     %f\n',D4B10,B4B10)
fprintf('Octanol      %f     %f\n',D4Boct,B4Boct)
fprintf('Total        %f   %f\n\n',D4Btotal,B4Btotal)
fprintf('Temp         %f   %f degC\n',tD4B,tB4B)
fprintf('Pressure     %f   %f psi\n\n',pD4B,pB4B)
fprintf('N  = %f\n',Nactual4B)
fprintf('Nmin  = %f\n',Nmin4B)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns4B)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr4B)
fprintf('R = %f\n\n',R4B)
 
 
%------------------- COLUMN 2A -----------------------------------------
% feed = heavies from C1 = Btotal=B2+B6+B8+B10+Boct+Bcyc+Beth+Bmeth;
% majority = cyclohex, octene, decene, octanol (some hexene)
% temp = tB
% pressure = pB
 
% determine q of feed
pflash2=pB; % pressure of stream
K2Bflash=p2B/pflash2;
K6Bflash=p6B/pflash2;
K10Bflash=p10B/pflash2;
K8Bflash=p8B/pflash2;
KcycBflash=pcycB/pflash2;
KoctBflash=poctB/pflash2;
KethBflash=pethB/pflash2;
KmethBflash=pmethB/pflash2;
alphaF2Ameth=KmethBflash/K8Bflash;
alphaF2A2=K2Bflash/K8Bflash;
alphaF2Aeth=KethBflash/K8Bflash;
alphaF2A6=K6Bflash/K8Bflash;
alphaF2Acyc=KcycBflash/K8Bflash;
alphaF2A10=K10Bflash/K8Bflash;
alphaF2A8=K8Bflash/K8Bflash; % HK=octene is reference
alphaF2Aoct=KoctBflash/K8Bflash;
psi2B=(zB2*(1-K2Bflash))/(1+x(18)*(K2Bflash-1));
psi6B=(zB6*(1-K6Bflash))/(1+x(18)*(K6Bflash-1));
psi10B=(zB10*(1-K10Bflash))/(1+x(18)*(K10Bflash-1));
psi8B=(zB8*(1-K8Bflash))/(1+x(18)*(K8Bflash-1));
psicycB=(zBcyc*(1-KcycBflash))/(1+x(18)*(KcycBflash-1));
psioctB=(zBoct*(1-KoctBflash))/(1+x(18)*(KoctBflash-1));
psimethB=(zBmeth*(1-KmethBflash))/(1+x(18)*(KmethBflash-1));
psiethB=(zBeth*(1-KethBflash))/(1+x(18)*(KethBflash-1));
y(18)=psi2B+psi6B+psi10B+psi8B+psicycB+psioctB+psimethB+psiethB; % solve for psi, gives us VLE
 
xB2=zB2/(1+x(18)*(K2Bflash-1)); % mol fraction hexene in liq phase
yB2=xB2*K2Bflash;
xB6=zB6/(1+x(18)*(K6Bflash-1));
yB6=xB6*K6Bflash;
xBcyc=zBcyc/(1+x(18)*(KcycBflash-1));
yBcyc=xBcyc*KcycBflash;
VB=(Btotal)*x(18);
LB=(Btotal)-VB;
qB=1-(VB/(VB+LB));
 
% specify split of 2 keys (LK=cyclohexane, HK=octene)       SPECIFY
B2Acyc=.01;           %LK=cyclohexane
D2Acyc=Bcyc-B2Acyc;
D2A8=0.01;            %HK=octene         
B2A8=B8-D2A8;
 
% estimate split of nonkeys
Bmethe2A=0;                  % light
Dmethe2A=Bmeth-Bmethe2A;
B2e2A=0;
D2e2A=B2-B2e2A;
Bethe2A=0;                
Dethe2A=Beth-Bethe2A;
B6e2A=0;                  % light
D6e2A=B6-B6e2A;
D10e2A=0;                % heavies
B10e2A=B10-D10e2A;
Docte2A=0;
Bocte2A=Boct-Docte2A;
 
Dtotale2A=D2e2A+Dethe2A+Dmethe2A+D6e2A+D2Acyc+D2A8+D10e2A+Docte2A;
Btotale2A=B2e2A+Bethe2A+Bmethe2A+B6e2A+B2Acyc+B2A8+B10e2A+Bocte2A;
 
% estimated Distillate mol fractions
zD6e2A=D6e2A/Dtotale2A;
zD2Acyc=D2Acyc/Dtotale2A;
zD2e2A=D2e2A/Dtotale2A; %LK
zD2A8=D2A8/Dtotale2A;
zD10e2A=D10e2A/Dtotale2A;
zDocte2A=Docte2A/Dtotale2A;
zDethe2A=Dethe2A/Dtotale2A; %HK
zDmethe2A=Dmethe2A/Dtotale2A;
% estimated Bottoms mol fractions
zB6e2A=B6e2A/Btotale2A;
zB2Acyc=B2Acyc/Btotale2A;
zB2e2A=B2e2A/Btotale2A; %LK
zB2A8=B2A8/Btotale2A;
zB10e2A=B10e2A/Btotale2A;
zBocte2A=Bocte2A/Btotale2A;
zBethe2A=Bethe2A/Btotale2A; %HK
zBmethe2A=Bmethe2A/Btotale2A;
 
pB2A=27; % bottoms pressure 2           SPECIFY
pD2A=27; % distillate pressure 2        SPECIFY
% distillate calc
tD2A=x(19); % distillate temp
TD2A=tD2A+273.15;
pmethD2A=0.019337*exp(15.2243-597.84/(TD2A-7.16)); %vapor P, methane 
p2D2A=(-0.0573)*(tD2A^2)+21.523*tD2A+461.88; %vapor P,ethylene,psi (depends on t)
pethD2A=0.019337*exp(15.6637-1511.42/(TD2A-17.16)); %vapor P, ethane 
p6D2A=0.019337*exp(15.8089-2654.81/(TD2A-47.3)); %vapor P, hexene, psi
pcycD2A=0.019337*(exp(15.7527-2766.63/(TD2A+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8D2A=0.019337*exp(15.963-3116.52/(TD2A-60.39)); %vapor P, octene
p10D2A=0.019337*exp(16.0129-3448.18/(TD2A-76.09)); %vapor P, decene
poctD2A=0.019337*exp(15.7428-3017.81/(TD2A-137.1));%vapor P, octanol, psi 
K2D2A=p2D2A/pD2A;
K6D2A=p6D2A/pD2A;
K10D2A=p10D2A/pD2A;
K8D2A=p8D2A/pD2A;
KcycD2A=pcycD2A/pD2A;
KoctD2A=poctD2A/pD2A;
KethD2A=pethD2A/pD2A;
KmethD2A=pmethD2A/pD2A;
y(19)=1-(zD2e2A*K2D2A+zD6e2A*K6D2A+zD10e2A*K10D2A+zD2A8*K8D2A+zD2Acyc*KcycD2A+zDocte2A*KoctD2A+zDethe2A*KethD2A+zDmethe2A*KmethD2A); %calculate bubble point temp, total CONDENSER, degC, DISTILLATE
alphaD2A6=K6D2A/K8D2A; % HK=octene is reference
alphaD2A2=K2D2A/K8D2A;
alphaD2A8=K8D2A/K8D2A;
alphaD2A10=K10D2A/K8D2A;
alphaD2Aoct=KoctD2A/K8D2A;
alphaD2Aeth=KethD2A/K8D2A;
alphaD2Ameth=KmethD2A/K8D2A;
alphaD2Acyc=KcycD2A/K8D2A;
 
% bottoms calc
tB2A=x(20); % bottoms temp
TB2A=tB2A+273.15;
pmethB2A=0.019337*exp(15.2243-597.84/(TB2A-7.16)); %vapor P, methane 
p2B2A=(-0.0573)*(tB2A^2)+21.523*tB2A+461.88; %vapor P,ethylene,psi (depends on t)
pethB2A=0.019337*exp(15.6637-1511.42/(TB2A-17.16)); %vapor P, ethane 
p6B2A=0.019337*exp(15.8089-2654.81/(TB2A-47.3)); %vapor P, hexene, psi
pcycB2A=0.019337*(exp(15.7527-2766.63/(TB2A+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8B2A=0.019337*exp(15.963-3116.52/(TB2A-60.39)); %vapor P, octene
p10B2A=0.019337*exp(16.0129-3448.18/(TB2A-76.09)); %vapor P, decene
poctB2A=0.019337*exp(15.7428-3017.81/(TB2A-137.1));%vapor P, octanol, psi 
K2B2A=p2B2A/pB2A;
K6B2A=p6B2A/pB2A;
K10B2A=p10B2A/pB2A;
K8B2A=p8B2A/pB2A;
KcycB2A=pcycB2A/pB2A;
KoctB2A=poctB2A/pB2A;
KethB2A=pethB2A/pB2A;
KmethB2A=pmethB2A/pB2A;
y(20)=1-(zB2e2A*K2B2A+zB6e2A*K6B2A+zB10e2A*K10B2A+zB2A8*K8B2A+zB2Acyc*KcycB2A+zBocte2A*KoctB2A+zBethe2A*KethB2A+zBmethe2A*KmethB2A); %calculate bubble point temp, distillation, degC, BOTTOMS
alphaB2A6=K6B2A/K8B2A; % HK=octene is reference
alphaB2A2=K2B2A/K8B2A;
alphaB2A8=K8B2A/K8B2A;
alphaB2A10=K10B2A/K8B2A;
alphaB2Aoct=KoctB2A/K8B2A;
alphaB2Aeth=KethB2A/K8B2A;
alphaB2Ameth=KmethB2A/K8B2A;
alphaB2Acyc=KcycB2A/K8B2A;
 
% mean relative volatilities:
alphaAVG2A6=(alphaD2A6*alphaB2A6)^.5;  % HK=octene is reference
alphaAVG2A2=(alphaD2A2*alphaB2A2)^.5;
alphaAVG2A8=(alphaD2A8*alphaB2A8)^.5;
alphaAVG2A10=(alphaD2A10*alphaB2A10)^.5;
alphaAVG2Aoct=(alphaD2Aoct*alphaB2Aoct)^.5;
alphaAVG2Aeth=(alphaD2Aeth*alphaB2Aeth)^.5;
alphaAVG2Ameth=(alphaD2Ameth*alphaB2Ameth)^.5;
alphaAVG2Acyc=(alphaD2Acyc*alphaB2Acyc)^.5;
 
Nmin3=log((D2Acyc/D2A8)*(B2A8/B2Acyc))/log(alphaAVG2Acyc);
 
% nonkey component calculations: calculate using smaller of bi and di
%lights
B2Ameth=Bmeth/(1+(D2A8/B2A8)*((alphaAVG2Ameth)^Nmin3));
D2Ameth=Bmeth-B2Ameth;
B2A2=B2/(1+(D2A8/B2A8)*((alphaAVG2A2)^Nmin3)); % mol/s ethylene in bottoms HK=octene is reference
D2A2=B2-B2A2;
B2Aeth=Beth/(1+(D2A8/B2A8)*((alphaAVG2Aeth)^Nmin3));
D2Aeth=Beth-B2Aeth;
B2A6=B6/(1+(D2A8/B2A8)*((alphaAVG2A6)^Nmin3));
D2A6=B6-B2A6;
%heavies
D2A10=(B10*(D2A8/B2A8)*((alphaAVG2A10)^Nmin3))/(1+(D2A8/B2A8)*((alphaAVG2A10)^Nmin3));
B2A10=B10-D2A10;
D2Aoct=(Boct*(D2A8/B2A8)*((alphaAVG2Aoct)^Nmin3))/(1+(D2A8/B2A8)*((alphaAVG2Aoct)^Nmin3));
B2Aoct=Boct-D2Aoct;
 
D2Atotal=D2Ameth+D2A2+D2Aeth+D2A6+D2Acyc+D2A10+D2A8+D2Aoct;
B2Atotal=B2Ameth+B2A2+B2Aeth+B2A6+B2Acyc+B2A10+B2A8+B2Aoct;
 
% Bottoms mol fractions
zB62A=B2A6/B2Atotal;
zB2Acyc=B2Acyc/B2Atotal; %LK
zB2A2=B2A2/B2Atotal; 
zB2A8=B2A8/B2Atotal; %HK
zB2A10=B2A10/B2Atotal;
zB2Aoct=B2Aoct/B2Atotal;  
zB2Aeth=B2Aeth/B2Atotal;
zB2Ameth=B2Ameth/B2Atotal;
 
% Underwood
theta2A10=(alphaF2A10*zB10)/(alphaF2A10-x(21));
theta2A2=(alphaF2A2*zB2)/(alphaF2A2-x(21));
theta2A6=(alphaF2A6*zB6)/(alphaF2A6-x(21));
theta2A8=(alphaF2A8*zB8)/(alphaF2A8-x(21));
theta2Acyc=(alphaF2Acyc*zBcyc)/(alphaF2Acyc-x(21));
theta2Aeth=(alphaF2Aeth*zBeth)/(alphaF2Aeth-x(21));
theta2Ameth=(alphaF2Ameth*zBmeth)/(alphaF2Ameth-x(21));
theta2Aoct=(alphaF2Aoct*zBoct)/(alphaF2Aoct-x(21));
y(21)=(1-qB)-(theta2A10+theta2A2+theta2A6+theta2A8+theta2Acyc+theta2Aeth+theta2Ameth+theta2Aoct); % VERY sensitive to initial conditions
Theta2A10=(alphaF2A10*(D2A10/D2Atotal))/(alphaF2A10-x(21));
Theta2A2=(alphaF2A2*(D2A2/D2Atotal))/(alphaF2A2-x(21));
Theta2A6=(alphaF2A6*(D2A6/D2Atotal))/(alphaF2A6-x(21));
Theta2A8=(alphaF2A8*(D2A8/D2Atotal))/(alphaF2A8-x(21));
Theta2Acyc=(alphaF2Acyc*(D2Acyc/D2Atotal))/(alphaF2Acyc-x(21));
Theta2Aeth=(alphaF2Aeth*(D2Aeth/D2Atotal))/(alphaF2Aeth-x(21));
Theta2Ameth=(alphaF2Ameth*(D2Ameth/D2Atotal))/(alphaF2Ameth-x(21));
Theta2Aoct=(alphaF2Aoct*(D2Aoct/D2Atotal))/(alphaF2Aoct-x(21));
 
Rmin3=(Theta2A10+Theta2A2+Theta2A6+Theta2A8+Theta2Acyc+Theta2Aeth+Theta2Ameth+Theta2Aoct)-1;
 
% GILLILAND EQNS
R3=1.5*Rmin3; % Assume optimal R/Rmin ratio         SPECIFY
XX3=(R3-Rmin3)/(R3+1);
YY3=1-exp(((1+54.4*XX3)/(11+117.2*XX3))*((XX3-1)/(XX3^0.5)));
Nactual3=(-Nmin3-YY3)/(YY3-1);
 
% KIRKBRIDE EQN
stageratio3=((zB8/zBcyc)*(((B2Acyc/B2Atotal)/(D2A8/D2Atotal))^2)*(Btotal/Dtotal))^0.206;
Ns3=Nactual3/(stageratio3+1); % # stages in stripping section (below feed stage)
Nr3=stageratio3*Ns3;          % # stages in rectifying section (above feed stage)
 
fprintf('----------------------------------DISTILLATION COLUMN 2A\n')
fprintf('Column 2A    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D2Ameth,B2Ameth)
fprintf('Ethylene     %f   %f  mol/s\n',D2A2,B2A2)
fprintf('Ethane       %f     %f\n',D2Aeth,B2Aeth)
fprintf('Hexene       %f    %f\n',D2A6,B2A6)
fprintf('Cyclohexane  %f     %f\n',D2Acyc,B2Acyc)
fprintf('Octene       %f     %f\n',D2A8,B2A8)
fprintf('Decene       %f     %f\n',D2A10,B2A10)
fprintf('Octanol      %f     %f\n',D2Aoct,B2Aoct)
fprintf('Total        %f   %f\n\n',D2Atotal,B2Atotal)
fprintf('Temp         %f   %f degC\n',tD2A,tB2A)
fprintf('Pressure     %f   %f psi\n\n',pD2A,pB2A)
fprintf('N  = %f\n',Nactual3)
fprintf('Nmin  = %f\n',Nmin3)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns3)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr3)
fprintf('R = %f\n\n',R3)
 
 
%------------------- COLUMN 3A -----------------------------------------
% feed = heavies from C2A = B2Atotal=B2Ameth+B2A2+B2Aeth+B2A6+B2Acyc+B2A10+B2A8+B2Aoct;
% majority = octene, decene, octanol (some cyclohex)
% temp = tB2A
% pressure = pB2A
 
% determine q of feed
pflash3=pB2A; % pressure of stream
K2B2flash=p2B2A/pflash3;
K6B2flash=p6B2A/pflash3;
K10B2flash=p10B2A/pflash3;
K8B2flash=p8B2A/pflash3;
KcycB2flash=pcycB2A/pflash3;
KoctB2flash=poctB2A/pflash3;
KethB2flash=pethB2A/pflash3;
KmethB2flash=pmethB2A/pflash3;
alphaF3Ameth=KmethB2flash/K10B2flash;
alphaF3A2=K2B2flash/K10B2flash;
alphaF3Aeth=KethB2flash/K10B2flash;
alphaF3A6=K6B2flash/K10B2flash;
alphaF3Acyc=KcycB2flash/K10B2flash;
alphaF3A10=K10B2flash/K10B2flash; % HK=decene is reference
alphaF3A8=K8B2flash/K10B2flash;  % LK=octene
alphaF3Aoct=KoctB2flash/K10B2flash;
psi2B3A=(zB2A2*(1-K2B2flash))/(1+x(22)*(K2B2flash-1));
psi6B3A=(zB62A*(1-K6B2flash))/(1+x(22)*(K6B2flash-1));
psi10B3A=(zB2A10*(1-K10B2flash))/(1+x(22)*(K10B2flash-1));
psi8B3A=(zB2A8*(1-K8B2flash))/(1+x(22)*(K8B2flash-1));
psicycB3A=(zB2Acyc*(1-KcycB2flash))/(1+x(22)*(KcycB2flash-1));
psioctB3A=(zB2Aoct*(1-KoctB2flash))/(1+x(22)*(KoctB2flash-1));
psimethB3A=(zB2Ameth*(1-KethB2flash))/(1+x(22)*(KethB2flash-1));
psiethB3A=(zB2Aeth*(1-KmethB2flash))/(1+x(22)*(KmethB2flash-1));
y(22)=psi2B3A+psi6B3A+psi10B3A+psi8B3A+psicycB3A+psioctB3A+psimethB3A+psiethB3A; % solve for psi, gives us VLE
VF3A=(B2Atotal)*x(22);
LF3A=(B2Atotal)-VF3A;
qF3A=1-(VF3A/(VF3A+LF3A));
 
% specify split of 2 keys (LK=octene, HK=decene)        SPECIFY
B3A8=.001;         %LK=octene (8) 
D3A8=B2A8-B3A8;
D3A10=.001;        %HK=decene (10)
B3A10=B2A10-D3A10;
 
% estimate split of nonkeys
% light
B3Aemeth=0;
D3Aemeth=B2Ameth-B3Aemeth;
B3Ae2=0;
D3Ae2=B2A2-B3Ae2;
B3Aeeth=0;
D3Aeeth=B2Aeth-B3Aeeth;
B3Ae6=0;
D3Ae6=B2A6-B3Ae6;
B3Aecyc=0;
D3Aecyc=B2Acyc-B3Aecyc;
% heavy
D3Aeoct=0;
B3Aeoct=B2Aoct-D3Aeoct;
 
B3Aetotal=B3A10+B3A8+B3Ae2+B3Ae6+B3Aecyc+B3Aeeth+B3Aemeth+B3Aeoct;
D3Aetotal=D3A10+D3A8+D3Ae2+D3Ae6+D3Aecyc+D3Aeeth+D3Aemeth+D3Aeoct;
 
% estimated Distillate mol fractions
zD3Aemeth=D3Aemeth/D3Aetotal;
zD3Ae2=D3Ae2/D3Aetotal;
zD3Aeeth=D3Aeeth/D3Aetotal;
zD3Ae6=D3Ae6/D3Aetotal;
zD3Aecyc=D3Aecyc/D3Aetotal;
zD3Ae8=D3A8/D3Aetotal;
zD3Ae10=D3A10/D3Aetotal;
zD3Aeoct=D3Aeoct/D3Aetotal;
% estimated Bottoms mol fractions
zB3Aemeth=B3Aemeth/B3Aetotal;
zB3Ae2=B3Ae2/B3Aetotal;
zB3Aeeth=B3Aeeth/B3Aetotal;
zB3Ae6=B3Ae6/B3Aetotal;
zB3Aecyc=B3Aecyc/B3Aetotal;
zB3Ae8=B3A8/B3Aetotal;
zB3Ae10=B3A10/B3Aetotal;
zB3Aeoct=B3Aeoct/B3Aetotal;
 
pB3A=22; % bottoms pressure         SPECIFY        
pD3A=22; % distillate pressure      SPECIFY    
% 3A distillate calc
tD3A=x(23); % distillate temp
TD3A=tD3A+273.15;
pmethD3A=0.019337*exp(15.2243-597.84/(TD3A-7.16)); %vapor P, methane
p2D3A=(-0.0573)*(tD3A^2)+21.523*tD3A+461.88; %vapor P,ethylene,psi (depends on t)
pethD3A=0.019337*exp(15.6637-1511.42/(TD3A-17.16)); %vapor P, ethane 
p6D3A=0.019337*exp(15.8089-2654.81/(TD3A-47.3)); %vapor P, hexene, psi
pcycD3A=0.019337*(exp(15.7527-2766.63/(TD3A+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8D3A=0.019337*exp(15.963-3116.52/(TD3A-60.39)); %vapor P, octene
p10D3A=0.019337*exp(16.0129-3448.18/(TD3A-76.09)); %vapor P, decene
poctD3A=0.019337*exp(15.7428-3017.81/(TD3A-137.1));%vapor P, octanol, psi
K2D3A=p2D3A/pD3A;
K6D3A=p6D3A/pD3A;
K10D3A=p10D3A/pD3A;
K8D3A=p8D3A/pD3A;
KcycD3A=pcycD3A/pD3A;
KoctD3A=poctD3A/pD3A;
KethD3A=pethD3A/pD3A;
KmethD3A=pmethD3A/pD3A;
y(23)=1-(zD3Ae10*K10D3A+zD3Ae2*K2D3A+zD3Ae6*K6D3A+zD3Ae8*K8D3A+zD3Aecyc*KcycD3A+zD3Aeeth*KethD3A+zD3Aemeth*KmethD3A+zD3Aeoct*KoctD3A); %calculate bubble point temp, total CONDENSER, degC, DISTILLATE
alphaD3A2=K2D3A/K10D3A;
alphaD3A6=K6D3A/K10D3A;
alphaD3A10=K10D3A/K10D3A; % HK=decene is reference
alphaD3A8=K8D3A/K10D3A;
alphaD3Acyc=KcycD3A/K10D3A;
alphaD3Aoct=KoctD3A/K10D3A;
alphaD3Aeth=KethD3A/K10D3A;
alphaD3Ameth=KmethD3A/K10D3A;
 
% bottoms calc
tB3A=x(24); % distillate temp
TB3A=tB3A+273.15;
pmethB3A=0.019337*exp(15.2243-597.84/(TB3A-7.16)); %vapor P, methane
p2B3A=(-0.0573)*(tB3A^2)+21.523*tB3A+461.88; %vapor P,ethylene,psi (depends on t)
pethB3A=0.019337*exp(15.6637-1511.42/(TB3A-17.16)); %vapor P, ethane
p6B3A=0.019337*exp(15.8089-2654.81/(TB3A-47.3)); %vapor P, hexene, psi
pcycB3A=0.019337*(exp(15.7527-2766.63/(TB3A+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
p8B3A=0.019337*exp(15.963-3116.52/(TB3A-60.39)); %vapor P, octene
p10B3A=0.019337*exp(16.0129-3448.18/(TB3A-76.09)); %vapor P, decene
poctB3A=0.019337*exp(15.7428-3017.81/(TB3A-137.1));%vapor P, octanol, psi
K2B3A=p2B3A/pB3A;
K6B3A=p6B3A/pB3A;
K10B3A=p10B3A/pB3A;
K8B3A=p8B3A/pB3A;
KcycB3A=pcycB3A/pB3A;
KoctB3A=poctB3A/pB3A;
KethB3A=pethB3A/pB3A;
KmethB3A=pmethB3A/pB3A;
y(24)=1-(zB3Ae10*K10B3A+zB3Ae2*K2B3A+zB3Ae6*K6B3A+zB3Ae8*K8B3A+zB3Aecyc*KcycB3A+zB3Aeeth*KethB3A+zB3Aemeth*KmethB3A+zB3Aeoct*KoctB3A); %calculate bubble point temp, distillation, degC, BOTTOMS
alphaB3A2=K2B3A/K10B3A;
alphaB3A6=K6B3A/K10B3A;
alphaB3A10=K10B3A/K10B3A; % HK=decene is reference
alphaB3A8=K8B3A/K10B3A;
alphaB3Acyc=KcycB3A/K10B3A;
alphaB3Aoct=KoctB3A/K10B3A;
alphaB3Aeth=KethB3A/K10B3A;
alphaB3Ameth=KmethB3A/K10B3A;
 
% mean relative volatilities:
alphaAVG3A2 =(alphaD3A2*alphaB3A2)^.5;
alphaAVG3A6 =(alphaD3A6*alphaB3A6)^.5;
alphaAVG3A10=(alphaD3A10*alphaB3A10)^.5; % HK=decene is reference
alphaAVG3A8 =(alphaD3A8*alphaB3A8)^.5;
alphaAVG3Acyc=(alphaD3Acyc*alphaB3Acyc)^.5;
alphaAVG3Aoct=(alphaD3Aoct*alphaB3Aoct)^.5;
alphaAVG3Aeth=(alphaD3Aeth*alphaB3Aeth)^.5;
alphaAVG3Ameth=(alphaD3Ameth*alphaB3Ameth)^.5;
 
Nmin3A=log((D3A8/D3A10)*(B3A10/B3A8))/log(alphaAVG3A8);
 
% nonkey component calculations: calculate using smaller of bi and di
%lights
B3Ameth=B2Ameth/(1+(D3A10/B3A10)*((alphaAVG3Ameth)^Nmin3A));
D3Ameth=B2Ameth-B3Ameth;
B3A2=B2A2/(1+(D3A10/B3A10)*((alphaAVG3A2)^Nmin3A));
D3A2=B2A2-B3A2;
B3Aeth=B2Aeth/(1+(D3A10/B3A10)*((alphaAVG3Aeth)^Nmin3A));
D3Aeth=B2Aeth-B3Aeth;
B3A6=B2A6/(1+(D3A10/B3A10)*((alphaAVG3A6)^Nmin3A));
D3A6=B2A6-B3A6;
B3Acyc=B2Acyc/(1+(D3A10/B3A10)*((alphaAVG3Acyc)^Nmin3A));
D3Acyc=B2Acyc-B3Acyc;
%heavies
D3Aoct=(B2Aoct*(D3A10/B3A10)*((alphaAVG3Aoct)^Nmin3A))/(1+(D3A10/B3A10)*((alphaAVG3Aoct)^Nmin3A));
B3Aoct=B2Aoct-D3Aoct;
 
D3Atotal=D3Ameth+D3A2+D3Aeth+D3A6+D3Acyc+D3A10+D3A8+D3Aoct;
B3Atotal=B3Ameth+B3A2+B3Aeth+B3A6+B3Acyc+B3A10+B3A8+B3Aoct;
 
% Underwood
theta3A10=(alphaF3A10*zB2A10)/(alphaF3A10-x(25));
theta3A2=(alphaF3A2*zB2A2)/(alphaF3A2-x(25));
theta3A6=(alphaF3A6*zB62A)/(alphaF3A6-x(25));
theta3A8=(alphaF3A8*zB2A8)/(alphaF3A8-x(25));
theta3Acyc=(alphaF3Acyc*zB2Acyc)/(alphaF3Acyc-x(25));
theta3Aeth=(alphaF3Aeth*zB2Aeth)/(alphaF3Aeth-x(25));
theta3Ameth=(alphaF3Ameth*zB2Ameth)/(alphaF3Ameth-x(25));
theta3Aoct=(alphaF3Aoct*zB2Aoct)/(alphaF3Aoct-x(25));
y(25)=(1-qF3A)-(theta3A10+theta3A2+theta3A6+theta3A8+theta3Acyc+theta3Aeth+theta3Ameth+theta3Aoct);
Theta3A10=(alphaF3A10*(D3A10/D3Atotal))/(alphaF3A10-x(25));
Theta3A2=(alphaF3A2*(D3A2/D3Atotal))/(alphaF3A2-x(25));
Theta3A6=(alphaF3A6*(D3A6/D3Atotal))/(alphaF3A6-x(25));
Theta3A8=(alphaF3A8*(D3A8/D3Atotal))/(alphaF3A8-x(25));
Theta3Acyc=(alphaF3Acyc*(D3Acyc/D3Atotal))/(alphaF3Acyc-x(25));
Theta3Aeth=(alphaF3Aeth*(D3Aeth/D3Atotal))/(alphaF3Aeth-x(25));
Theta3Ameth=(alphaF3Ameth*(D3Ameth/D3Atotal))/(alphaF3Ameth-x(25));
Theta3Aoct=(alphaF3Aoct*(D3Aoct/D3Atotal))/(alphaF3Aoct-x(25));
 
Rmin3A=(Theta3A10+Theta3A2+Theta3A6+Theta3A8+Theta3Acyc+Theta3Aeth+Theta3Ameth+Theta3Aoct)-1;
 
% GILLILAND EQNS
R3A=1.5*Rmin3A; % Assume optimal R/Rmin ratio       SPECIFY
XX3A=(R3A-Rmin3A)/(R3A+1);
YY3A=1-exp(((1+54.4*XX3A)/(11+117.2*XX3A))*((XX3A-1)/(XX3A^0.5)));
Nactual3A=(-Nmin3A-YY3A)/(YY3A-1);
 
% KIRKBRIDE EQN
stageratio3A=((zB2A10/zB2A8)*(((B3A8/B3Atotal)/(D3A10/D3Atotal))^2)*(B3Atotal/D3Atotal))^0.206;
Ns3A=Nactual3A/(stageratio3A+1); % # stages in stripping section (below feed stage)
Nr3A=stageratio3A*Ns3A;          % # stages in rectifying section (above feed stage)
 
octenepurity=(D3A8*mw8)/(D3Ameth*mwmeth+D3A2*mw2+D3Aeth*mweth+D3A6*mw6+D3Acyc*mwcyc+D3A10*mw10+D3A8*mw8+D3Aoct*mwoct);
 
fprintf('----------------------------------DISTILLATION COLUMN 3A\n')
fprintf('Column 3A    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D3Ameth,B3Ameth)
fprintf('Ethylene     %f   %f  mol/s\n',D3A2,B3A2)
fprintf('Ethane       %f     %f\n',D3Aeth,B3Aeth)
fprintf('Hexene       %f    %f\n',D3A6,B3A6)
fprintf('Cyclohexane  %f     %f\n',D3Acyc,B3Acyc)
fprintf('Octene       %f     %f\n',D3A8,B3A8)
fprintf('Decene       %f     %f\n',D3A10,B3A10)
fprintf('Octanol      %f     %f\n',D3Aoct,B3Aoct)
fprintf('Total        %f   %f\n\n',D3Atotal,B3Atotal)
fprintf('Temp         %f   %f degC\n',tD3A,tB3A)
fprintf('Pressure     %f   %f psi\n\n',pD3A,pB3A)
fprintf('N  = %f\n',Nactual3A)
fprintf('Nmin  = %f\n',Nmin3A)
fprintf('Octene Purity = %f\n',octenepurity)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns3A)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr3A)
fprintf('R = %f\n\n',R3A)
fprintf('*** ****** ***** ******** ** ******* **** ***** ****** ***** ****** *****\n')
 
fprintf('-----------lb/h-----------------------DISTILLATION COLUMN 1\n')
fprintf('Column 1     Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',Dmeth*mwmeth*3600*0.00220462,Bmeth*mwmeth*3600*0.00220462)
fprintf('Ethylene     %f   %f  lb/h\n',D2*mw2*3600*0.00220462,B2*mw2*3600*0.00220462)
fprintf('Ethane       %f     %f\n',Deth*mweth*3600*0.00220462,Beth*mweth*3600*0.00220462)
fprintf('Hexene       %f    %f\n',D6*mw6*3600*0.00220462,B6*mw6*3600*0.00220462)
fprintf('Cyclohexane  %f     %f\n',Dcyc*mwcyc*3600*0.00220462,Bcyc*mwcyc*3600*0.00220462)
fprintf('Octene       %f     %f\n',D8*mw8*3600*0.00220462,B8*mw8*3600*0.00220462)
fprintf('Decene       %f     %f\n',D10*mw10*3600*0.00220462,B10*mw10*3600*0.00220462)
fprintf('Octanol      %f     %f\n',Doct*mwoct*3600*0.00220462,Boct*mwoct*3600*0.00220462)
fprintf('Total        %f   %f\n\n',Dtotal,Btotal)
fprintf('Temp         %f   %f \n',tD*1.8+32,tB*1.8+32)
fprintf('Pressure     %f   %f psi\n\n',pD,pB)
fprintf('N  = %f\n',Nactual)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr)
fprintf('R = %f\n\n',R)
fprintf('Nonkey Error = %f\n\n',nonkeyerror)
fprintf('-----------lb/h-----------------------DISTILLATION COLUMN 2B\n')
fprintf('Column 2B    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D2meth*mwmeth*3600*0.00220462,B2meth*mwmeth*3600*0.00220462)
fprintf('Ethylene     %f   %f  \n',D22*mw2*3600*0.00220462,B22*mw2*3600*0.00220462)
fprintf('Ethane       %f     %f\n',D2eth*mweth*3600*0.00220462,B2eth*mweth*3600*0.00220462)
fprintf('Hexene       %f    %f\n',D26*mw6*3600*0.00220462,B26*mw6*3600*0.00220462)
fprintf('Cyclohexane  %f     %f\n',D2cyc*mwcyc*3600*0.00220462,B2cyc*mwcyc*3600*0.00220462)
fprintf('Octene       %f     %f\n',D28*mw8*3600*0.00220462,B28*mw8*3600*0.00220462)
fprintf('Decene       %f     %f\n',D210*mw10*3600*0.00220462,B210*mw10*3600*0.00220462)
fprintf('Octanol      %f     %f\n',D2oct*mwoct*3600*0.00220462,B2oct*mwoct*3600*0.00220462)
fprintf('Total        %f   %f\n\n',D2total,B2total)
fprintf('Temp         %f   %f \n',tD2*1.8+32,tB2*1.8+32)
fprintf('Pressure     %f   %f psi\n\n',pD2,pB2)
fprintf('N  = %f\n',Nactual2B)
fprintf('Nmin  = %f\n',Nmin2)
fprintf('Hexene Purity = %f\n',estimatedhexenepurity)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns2B)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr2B)
fprintf('R = %f\n\n',R2B)
fprintf('---------lb/h-------------------------DISTILLATION COLUMN 3B\n')
fprintf('Column 3B    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D3Bmeth*mwmeth*3600*0.00220462,B3Bmeth*mwmeth*3600*0.00220462)
fprintf('Ethylene     %f   %f  \n',D3B2*mw2*3600*0.00220462,B3B2*mw2*3600*0.00220462)
fprintf('Ethane       %f     %f\n',D3Beth*mweth*3600*0.00220462,B3Beth*mweth*3600*0.00220462)
fprintf('Hexene       %f    %f\n',D3B6*mw6*3600*0.00220462,B3B6*mw6*3600*0.00220462)
fprintf('Cyclohexane  %f     %f\n',D3Bcyc*mwcyc*3600*0.00220462,B3Bcyc*mwcyc*3600*0.00220462)
fprintf('Octene       %f     %f\n',D3B8*mw8*3600*0.00220462,B3B8*mw8*3600*0.00220462)
fprintf('Decene       %f     %f\n',D3B10*mw10*3600*0.00220462,B3B10*mw10*3600*0.00220462)
fprintf('Octanol      %f     %f\n',D3Boct*mwoct*3600*0.00220462,B3Boct*mwoct*3600*0.00220462)
fprintf('Total        %f   %f\n\n',D3Btotal,B3Btotal)
fprintf('Temp         %f   %f \n',tD3B*1.8+32,tB3B*1.8+32)
fprintf('Pressure     %f   %f psi\n\n',pD3B,pB3B)
fprintf('N  = %f\n',Nactual3B)
fprintf('Nmin  = %f\n',Nmin3B)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns3B)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr3B)
fprintf('R = %f\n\n',R3B)
fprintf('---------lb/h-------------------------DISTILLATION COLUMN 4B\n')
fprintf('Column 4B    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D4Bmeth*mwmeth*3600*0.00220462,B4Bmeth*mwmeth*3600*0.00220462)
fprintf('Ethylene     %f   %f  \n',D4B2*mw2*3600*0.00220462,B4B2*mw2*3600*0.00220462)
fprintf('Ethane       %f     %f\n',D4Beth*mweth*3600*0.00220462,B4Beth*mweth*3600*0.00220462)
fprintf('Hexene       %f    %f\n',D4B6*mw6*3600*0.00220462,B4B6*mw6*3600*0.00220462)
fprintf('Cyclohexane  %f     %f\n',D4Bcyc*mwcyc*3600*0.00220462,B4Bcyc*mwcyc*3600*0.00220462)
fprintf('Octene       %f     %f\n',D4B8*mw8*3600*0.00220462,B4B8*mw8*3600*0.00220462)
fprintf('Decene       %f     %f\n',D4B10*mw10*3600*0.00220462,B4B10*mw10*3600*0.00220462)
fprintf('Octanol      %f     %f\n',D4Boct*mwoct*3600*0.00220462,B4Boct*mwoct*3600*0.00220462)
fprintf('Total        %f   %f\n\n',D4Btotal,B4Btotal)
fprintf('Temp         %f   %f \n',tD4B*1.8+32,tB4B*1.8+32)
fprintf('Pressure     %f   %f psi\n\n',pD4B,pB4B)
fprintf('N  = %f\n',Nactual4B)
fprintf('Nmin  = %f\n',Nmin4B)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns4B)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr4B)
fprintf('R = %f\n\n',R4B)
fprintf('------------lb/h----------------------DISTILLATION COLUMN 2A\n')
fprintf('Column 2A    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D2Ameth*mwmeth*3600*0.00220462,B2Ameth*mwmeth*3600*0.00220462)
fprintf('Ethylene     %f   %f  \n',D2A2*mw2*3600*0.00220462,B2A2*mw2*3600*0.00220462)
fprintf('Ethane       %f     %f\n',D2Aeth*mweth*3600*0.00220462,B2Aeth*mweth*3600*0.00220462)
fprintf('Hexene       %f    %f\n',D2A6*mw6*3600*0.00220462,B2A6*mw6*3600*0.00220462)
fprintf('Cyclohexane  %f     %f\n',D2Acyc*mwcyc*3600*0.00220462,B2Acyc*mwcyc*3600*0.00220462)
fprintf('Octene       %f     %f\n',D2A8*mw8*3600*0.00220462,B2A8*mw8*3600*0.00220462)
fprintf('Decene       %f     %f\n',D2A10*mw10*3600*0.00220462,B2A10*mw10*3600*0.00220462)
fprintf('Octanol      %f     %f\n',D2Aoct*mwoct*3600*0.00220462,B2Aoct*mwoct*3600*0.00220462)
fprintf('Total        %f   %f\n\n',D2Atotal,B2Atotal)
fprintf('Temp         %f   %f \n',tD2A*1.8+32,tB2A*1.8+32)
fprintf('Pressure     %f   %f psi\n\n',pD2A,pB2A)
fprintf('N  = %f\n',Nactual3)
fprintf('Nmin  = %f\n',Nmin3)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns3)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr3)
fprintf('R = %f\n\n',R3)
 
fprintf('-------lb/h---------------------------DISTILLATION COLUMN 3A\n')
fprintf('Column 3A    Distillate   Bottoms\n')
fprintf('Methane      %f     %f\n',D3Ameth*mwmeth*3600*0.00220462,B3Ameth*mwmeth*3600*0.00220462)
fprintf('Ethylene     %f   %f  \n',D3A2*mw2*3600*0.00220462,B3A2*mw2*3600*0.00220462)
fprintf('Ethane       %f     %f\n',D3Aeth*mweth*3600*0.00220462,B3Aeth*mweth*3600*0.00220462)
fprintf('Hexene       %f    %f\n',D3A6*mw6*3600*0.00220462,B3A6*mw6*3600*0.00220462)
fprintf('Cyclohexane  %f     %f\n',D3Acyc*mwcyc*3600*0.00220462,B3Acyc*mwcyc*3600*0.00220462)
fprintf('Octene       %f     %f\n',D3A8*mw8*3600*0.00220462,B3A8*mw8*3600*0.00220462)
fprintf('Decene       %f     %f\n',D3A10*mw10*3600*0.00220462,B3A10*mw10*3600*0.00220462)
fprintf('Octanol      %f     %f\n',D3Aoct*mwoct*3600*0.00220462,B3Aoct*mwoct*3600*0.00220462)
fprintf('Total        %f   %f\n\n',D3Atotal,B3Atotal)
fprintf('Temp         %f   %f \n',tD3A*1.8+32,tB3A*1.8+32)
fprintf('Pressure     %f   %f psi\n\n',pD3A,pB3A)
fprintf('N  = %f\n',Nactual3A)
fprintf('Nmin  = %f\n',Nmin3A)
fprintf('Octene Purity = %f\n',octenepurity)
fprintf('Ns = %f = # stages in stripping section (below feed stage)\n',Ns3A)
fprintf('Nr = %f = # stages in rectifying section (above feed stage)\n',Nr3A)
fprintf('R = %f\n\n',R3A)
fprintf('*****************************end************************************************\n')
 
%fprintf('Var         Answer\n')
%fprintf('1           %f\n',x(1))
%fprintf('2           %f\n',x(2))
%fprintf('3           %f\n',x(3))
%fprintf('4           %f\n',x(4))
%fprintf('5           %f\n',x(5))
%fprintf('6           %f\n',x(6))
%fprintf('7           %f\n',x(7))
%fprintf('8           %f\n',x(8))
%fprintf('9           %f\n',x(9))
%fprintf('10           %f\n',x(10))
%fprintf('11           %f\n',x(11))
%fprintf('12           %f\n',x(12))
%fprintf('13          %f\n',x(13))
%fprintf('14           %f\n',x(14))
%fprintf('15           %f\n',x(15))
%fprintf('16           %f\n',x(16))
%fprintf('17           %f\n',x(17))
%fprintf('18           %f\n',x(18))
%fprintf('19           %f\n',x(19))
%fprintf('20           %f\n',x(20))
%fprintf('21           %f\n',x(21))
%fprintf('22           %f\n',x(22))
%fprintf('23           %f\n',x(23))
%fprintf('24           %f\n',x(24))
%fprintf('25           %f\n',x(25))
%fprintf('26           %f\n',x(26))
%fprintf('27           %f\n',x(27))
%fprintf('28           %f\n',x(28))
%fprintf('29           %f\n',x(29))
%fprintf('30          %f\n',x(30))
%fprintf('31          %f\n',x(31))
%fprintf('32          %f\n',x(32))
%fprintf('33           %f\n',x(33))
%-------------------DEAD CODE------------------------------------------
%Compressor 1 Calculations (PRIOR TO COLUMN 3B)
%Rc=0.296; % constant, kJ/(kg*K)
%cp=1.53 ;% constant, kJ/(kg*K)
%cv=1.23 ;% constant, kJ/(kg*K)
%kc=cp/cv;
%Pinlet=pD2*6894.75729; %inlet P to compressor 1, pascals
%pout1=186.5; %outlet P of compressor 1, psi 
%Poutlet=pout1*6894.75729; % convert to pascals
%mflowrate=(D2meth*mwmeth+D22*mw2+D2eth*mweth+D26*mw6+D2cyc*mwcyc+D210*mw10+D28*mw8+D2oct*mwoct)/1000;
                % mass flowrate, kg/s
%Tinlet=TD2; % inlet temp, kelvin
%Had=Rc*Tinlet*(kc/(kc-1))*(((Poutlet/Pinlet)^((kc-1)/kc))-1);
%Pad=mflowrate*Had % power needed, kilowatts
%flash
%-- ----- ------ -------- --------- ----- ---------
 
% PRESSURE ITERATION
% 1. calculate bubble P of distillate at 49 degC (cooling water temp)
%tcw=49;
%Tcw=tcw+273.15;
%p6cw=0.019337*exp(15.8089-2654.81/(Tcw-47.3)); %vapor P, hexene, psi
%pcyccw=0.019337*(exp(15.7527-2766.63/(Tcw+(-50.5)))); %vapor P,cyclohexane,psi (depends on T)
%p2cw=(-0.0573)*(tcw^2)+21.523*tcw+461.88; %vapor P,ethylene,psi (depends on t)
%p8cw=0.019337*exp(15.963-3116.52/(Tcw-60.39)); %vapor P, octene
%p10cw=0.019337*exp(16.0129-3448.18/(Tcw-76.09)); %vapor P, decene
%poctcw=0.019337*exp(15.7428-3017.81/(Tcw-137.1));%vapor P, octanol, psi
%Pcwbubble=zD6e*p6cw+zDcyce*pcyccw+zD2e*p2cw+zD8e*p8cw+zD10e*p10cw+zDocte*poctcw;
% 2. calculate dew P of distillate at 49 degC (cooling water temp)
%Pcwdew=1/(zD6e/p6cw+zDcyce/pcyccw+zD2e/p2cw+zD8e/p8cw+zD10e/p10cw+zDocte/poctcw);
% END PRESSURE ITERATION
 
%psi
%V=n25/.95;
%L=(n25+72+n105+n85+ncyc+noct)-V;
%phi2=2/106.1
%S2=K2*(V/L)
%y(11)=phi2-((S2/(S2-1))*abs(S2^(x(11)+1)-1))/(((S2/(S2-1))*(S2^(x(11)+1)-1))+(K2/(K2-1)))
%y(11)=phi2-((S2/(S2-1))*abs(S2^(x(11)+1)-1))/((S2/(S2-1))*(S2^(x(11)+1)-1)+(K2/(K2-1)))
%x(11)
%phi6
%phi10
%phi8
%phicyc
%phioct
 
%fprintf('P Before      P After \n%f    %f\n%f    %f\n%f    %f\n%f     %f\n%f\n',P1before,P1after,P2before,P2after,P3before,P3after,P4before,P4after,P5before)
%sum0=n20+ncyc; sum1=n21+x(3)+n101+n81+ncyc; sum2=n22+x(4)+n102+n82+ncyc;
%sum3=n23+x(5)+n103+n83+ncyc; sum4=n24+x(6)+n104+n84+ncyc; sum5=n25+72+n105+n85+ncyc; %sum of each stream, mol/s
%f20=n20/sum0; f60=0/sum0; f100=0; f80=0; fcyc0=ncyc/sum0; f0=f20+f60+fcyc0;
%f21=n21/sum1; f61=x(3)/sum1; f101=n101/sum1; f81=n81/sum1; fcyc1=ncyc/sum1; f1=f21+f61+f101+f81+fcyc1; % mol fraction per component
%f22=n22/sum2; f62=x(4)/sum2; f102=n102/sum2; f82=n82/sum2; fcyc2=ncyc/sum2; f2=f22+f62+f102+f82+fcyc2;
%f23=n23/sum3; f63=x(5)/sum3; f103=n103/sum3; f83=n83/sum3; fcyc3=ncyc/sum3; f3=f23+f63+f103+f83+fcyc3;
%f24=n24/sum4; f64=x(6)/sum4; f104=n104/sum4; f84=n84/sum4; fcyc4=ncyc/sum4; f4=f24+f64+f104+f84+fcyc4;
%f25=n25/sum5; f65=72/sum5;   f105=n105/sum5; f85=n85/sum5; fcyc5=ncyc/sum5; f5=f25+f65+f105+f85+fcyc5;
%P0=p2*f20+p6*f60+p10*f100+p8*f80+pcyc*fcyc0;
%P1=p2*f21+p6*f61+p10*f101+p8*f81+pcyc*fcyc1;
%P2=p2*f22+p6*f62+p10*f102+p8*f82+pcyc*fcyc2;
%P3=p2*f23+p6*f63+p10*f103+p8*f83+pcyc*fcyc3;
%P4=p2*f24+p6*f64+p10*f104+p8*f84+pcyc*fcyc4;
%P5=p2*f25+p6*f65+p10*f105+p8*f85+pcyc*fcyc5; % total calculated pressure of stream 5
