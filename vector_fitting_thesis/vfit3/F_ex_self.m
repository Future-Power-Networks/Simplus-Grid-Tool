%% data load
clear
load('GROUP3. 15V-5V-7V-10V.mat')
fSelect=1;
Bode_0=Bode_0(:,162*(fSelect-1)+1:162*fSelect);

%load('vfit3/myfile_0.mat')
Ns=length(Bode_0)/2;
nd=(1:Ns)*2-1;
nq=(1:Ns)*2;
dd_abs=Bode_0(3,nd); dd_phase=Bode_0(4,nd);
dq_abs=Bode_0(3,nq); dq_phase=Bode_0(4,nq);
qd_abs=Bode_0(5,nd); qd_phase=Bode_0(6,nd);
qq_abs=Bode_0(5,nq); qq_phase=Bode_0(6,nq);

omega=2*pi*Bode_0(2,nd);

ydd_=dd_abs .* exp(1i*deg2rad(dd_phase));
ydq_=dq_abs .* exp(1i*deg2rad(dq_phase));
yqd_=qd_abs .* exp(1i*deg2rad(qd_phase));
yqq_=qq_abs .* exp(1i*deg2rad(qq_phase));
Ysys_=[ydd_;ydq_;yqd_;yqq_];
%Ysys_=[ydd_;ydd_;yqq_];
%Ysys_=[ydd_;yqq_];

%% vecotr fitting
%Initial poles for Vector Fitting:
N=22; %order of approximation

%Initial poles: Complex conjugate pairs, linearly spaced:
bet=linspace(omega(1),omega(end),N/2);
poles=[];
for k=1:length(bet)
alf=-bet(k)*1e-2;
poles=[poles (alf-1i*bet(k)) (alf+1i*bet(k)) ];
end
%poles=-2*pi*logspace(0,4,N); 
weight=ones(1,Ns); %No weight
%weight=zeros(1,Ns); %Strong inverse weight
%for k=1:Ns
%weight(1,k)=1/sqrt(norm(f_mea(:,k)));
%end
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
opts.legend=1;     %include legends in plots



%% Forming (weighted) column sum:
g=0;
Nc=2;
for n=1:Nc
  %g=g+f(n,:); %unweighted sum     
  g=g+Ysys_(n,:)/norm(Ysys_(n,:));
  %g=g+f(n,:)/sqrt(norm(f(n,:)));     
end
weight_g=1./abs(g);

disp('****Calculating improved initial poles by fitting column sum ...')
Niter1=5;
for iter=1:Niter1
   disp(['   Iter ' num2str(iter)])
   if iter==Niter1,opts.skip_res=0; end
   [SER,poles,rmserr,fit]=vectfit3(g,1i*omega,poles,weight_g,opts);  
end

%% fitting
disp('vector fitting...')
%opts.skip_res=1;
weight=1./sqrt(abs(Ysys_));
for j=1:5
    [SER,poles,rmserr,fit]=vectfit3(Ysys_,1i*omega,poles,weight,opts); 
end

%SER=col2full(SER);
SER=tri2full(SER);
[Ysys_res,Ysys_pole] = ss2pr(SER.A,SER.B,SER.C);

disp('Done.')
Res_dd=reshape(Ysys_res(1,1,:),[length(Ysys_res),1]);
Res_qq=reshape(Ysys_res(2,2,:),[length(Ysys_res),1]);

%[num,den]=residue(Res_dd,Ysys_pole,[]);
%[z,p,k]=tf2zpk(num,den);
%[A,B,C,D]=zp2ss(z,p,k);
%figure(5);
%clf
%bode(A,B,C,D)
% specify frequency grid w or can generate it using bode
%[m,p,w] = bode(sys); % generate freqeuncy grid

%figure(10);
%bode(sysfrdpos,sysfrdneg) % plot both positive and negative frequencies

%figure(10);
%bode(sys_dd);