% ex2.m             
%
% -Creating an 18th order frequency response f(s) of 2 elements.
% -Fitting f(s) using vectfit3.m 
%   -Initial poles: 9 linearly spaced complex pairs (N=18)
%   -3 iterations
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%

clear all

D=0.2; E=2e-5;
p=[-4500  -41000 (-100+j*5e3) (-100-j*5e3)  (-120+j*15e3) (-120-j*15e3)  (-3e3+j*35e3) (-3e3-j*35e3)];
r=[-3000  -83000     (-5+j*7e3) (-5-j*7e3)  (-20+j*18e3) (-20-j*18e3)   (6e3+j*45e3) (6e3-j*45e3)];
p=[p (-200+j*45e3) (-200-j*45e3)  (-1500+j*45e3) (-1500-j*45e3)    ];
r=[r (40 +j*60e3) (40-j*60e3)     (90 +j*10e3) (90-j*10e3)      ];
p=[p (-5e2+j*70e3) (-5e2-j*70e3) (-1e3+j*73e3) (-1e3-j*73e3)  (-2e3+j*90e3) (-2e3-j*90e3)];
r=[r (5e4+j*80e3) (5e4-j*80e3)  (1e3+j*45e3) (1e3-j*45e3)     (-5e3+j*92e3) (-5e3-j*92e3)];      

w=2*pi*linspace(1,1e5,100); Ns=length(w); s=j*w;

p=2*pi*p;r=2*pi*r;
p1=p(1:10);r1=r(1:10);N1=length(p1);
p2=p(9:18);r2=r(9:18);N2=length(p2);


f=zeros(2,Ns);
for k=1:Ns
  for n=1:N1
    f(1,k)=f(1,k)+r1(n)/(s(k)-p1(n));  
  end
  f(1,k)=f(1,k)+s(k)*E;
end
f(1,:)=f(1,:)+D;


for k=1:Ns
  for n=1:N2
    f(2,k)=f(2,k)+r2(n)/(s(k)-p2(n));  
  end
  f(2,k)=f(2,k)+s(k)*3*E;
end
f(1,:)=f(1,:)+2*D;

%=====================================
% Rational function approximation of f(s):
%=====================================


N=18; %Order of approximation 

%Complex starting poles :
bet=linspace(w(1),w(Ns),N/2);
poles=[];
for n=1:length(bet)
  alf=-bet(n)*1e-2;
  poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ]; 
end

% Real starting poles :
%poles=-linspace(w(1),w(Ns),N); 
 
%Parameters for Vector Fitting : 

weight=ones(1,Ns);

opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=3;      %Include both D, E in fitting    
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=1;   %DO skip identification of residues (C,D,E) 
opts.cmplx_ss=0;   %Create real-only state space model

opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=1;       %Create magnitude plot for fitting of f(s) 
opts.logx=0;       %Use linear abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis 
opts.errplot=1;    %Include deviation in magnitude plot
opts.phaseplot=0;  %Do NOT produce plot of phase angle
opts.legend=1;     %Include legends in plots


disp('vector fitting...')
Niter=3;
for iter=1:Niter
  if iter==Niter, opts.skip_res=0; end
  disp(['   Iter ' num2str(iter)])
  [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
  rms(iter,1)=rmserr;
end
rms



