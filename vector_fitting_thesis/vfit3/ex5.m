% ex5.m             
%
% The program approximates f(s) with rational functions. f(s) is a vector of 5 elements
% represnting one column of the propagation matrix of a transmission line (parallell AC and DC line).
% The elements of the prop. matrix have been backwinded using a common time delay equal
% to the lossless time delay of the line. 
%
% -Reading frequency response f(s) from disk. (contains 5 elements)
% -Fitting f(s) using vectfit3.m 
% -Initial poles: 7 linearly spaced complex pairs (N=14)
% -5 iterations
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%

clear all

Ns=60; f=zeros(5,Ns); w=zeros(Ns,1); 
fid1=fopen('w.txt','r');
fid2=fopen('h.txt','r');
for k=1:Ns
  [w(k)]=fscanf(fid1,'%e',1);
  for n=1:5 
    [a1]=fscanf(fid2,'%e',1); [a2]=fscanf(fid2,'%e',1);
    f(n,k)=a1+i*a2;
  end
end
fclose(fid1);fclose(fid2);
s=i*w;

%=====================================
% Rational function approximation of f(s):
%=====================================


N=14; %Order of approximation 

%Complex starting poles :
bet=linspace(w(1),w(Ns),N/2);
poles=[];
for n=1:length(bet)
  alf=-bet(n)*1e-2;
  poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ]; 
end

%Parameters for Vector Fitting : 
weight=ones(1,Ns);
%weight=1./abs(f);

opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=1;      %Fitting with D=0, E=0 
opts.skip_pole=0;  %Do NOT skip pole identification
opts.skip_res=0;   %Do NOT skip identification of residues (C,D,E) 
opts.cmplx_ss=1;   %Create complex state space model

opts.spy1=0;       %No plotting for first stage of vector fitting
opts.spy2=1;       %Create magnitude plot for fitting of f(s) 
opts.logx=1;       %Use logarithmic abscissa axis
opts.logy=1;       %Use logarithmic ordinate axis 
opts.errplot=1;    %Include deviation in magnitude plot
opts.phaseplot=1;  %Also produce plot of phase angle (in addition to magnitiude)
opts.legend=0;     %Do NOT include legends in plots

Niter=5;
for iter =1:Niter
  if iter==Niter, opts.legend=1;end
  [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
  rms(iter,1)=rmserr;
end




