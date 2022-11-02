% ex1.m             
%
% Fitting an artificially created frequency response (single element)
%
% -Creating a 3rd order frequency response f(s)
% -Fitting f(s) using vectfit3.m 
%   -Initial poles: 3 logarithmically spaced real poles
%   -1 iteration
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%
clear all

%Frequency samples:
Ns=101;
s=2*pi*i*logspace(0,4,Ns); 

disp('Creating frequency response f(s)...') 
for k=1:Ns
  sk=s(k);
  f(1,k) = 2/(sk+5) + (30+j*40)/(sk-(-100+j*500)) ...
+ (30-j*40)/(sk-(-100-j*500)) + 0.5;
end


%Initial poles for Vector Fitting:
N=3; %order of approximation
poles=-2*pi*logspace(0,4,N); %Initial poles

weight=ones(1,Ns); %All frequency points are given equal weight

opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=3;      %Include both D, E in fitting    
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=0;   %Do NOT skip identification of residues (C,D,E) 
opts.cmplx_ss=1;   %Create complex state space model

opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=1;       %Create magnitude plot for fitting of f(s) 
opts.logx=1;       %Use logarithmic abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis 
opts.errplot=1;    %Include deviation in magnitude plot
opts.phaseplot=1;  %Also produce plot of phase angle (in addition to magnitiude)
opts.legend=1;     %Do include legends in plots


disp('vector fitting...')
[SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts); 
disp('Done.')

disp('Resulting state space model:')
A=full(SER.A)
B=SER.B
C=SER.C
D=SER.D
E=SER.E
rmserr 