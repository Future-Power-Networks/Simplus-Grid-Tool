% ex4b.m             
%
% Fitting all columns of the admittance matrix of a six-terminal system 
% (power system distribution network)
%
% -Reading frequency admittance matrix Y(s) from disk.
% -Extracting 1st column: f(s) (contains 6 elements)
% -Fitting g(s)=sum(f(s)) using vectfit3.m --> new initial poles
%   -Initial poles: 25 linearly spaced complex pairs (N=50)
%   -5 iterations
% -Fitting the six-columns one-by-one, stacking state-space model from each
%  column into global state-space model (3 iterations)
% -Plotting of result
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%
clear all


N=50; %order of approximation
Niter1=5; %Fitting column sum: n.o. iterations
Niter2=3; %Fitting column: n.o. iterations

disp('Reading data from file ...') %--> s(1,Ns), bigY(Nc,Nc,Ns)
fid1=fopen('fdne.txt','r');
Nc=fscanf(fid1,'%f',1);
Ns=fscanf(fid1,'%f',1);
bigY=zeros(Nc,Nc,Ns); s=zeros(1,Ns);
for k=1:Ns
  s(k)=fscanf(fid1,'%e',1);
  for row=1:Nc
    for col=1:Nc
      dum1=fscanf(fid1,'%e',1);
      dum2=fscanf(fid1,'%e',1);   
      bigY(row,col,k)=dum1+j*dum2;
    end
  end
end
s=i*s;
fclose(fid1);



tic
disp('-----------------S T A R T--------------------------')

AA=sparse([]); BB=sparse([]); CC=[]; DD=[]; EE=[];
bigf=[]; bigfit=[];

%Complex starting poles :
w=s/i;
bet=linspace(w(1),w(Ns),N/2);
poles=[];
for n=1:length(bet)
  alf=-bet(n)*1e-2;
  poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ]; 
end


for col=1:Nc %Column loop

  %Extracting elements in column:
  for n=1:Nc
    f(n,:)=squeeze(bigY(n,col,:)).';
  end  



  %Fitting options 
  opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint  
  opts.stable=1;     %Enforce stable poles
  opts.asymp=3;       %Fitting includes D and E
  opts.spy1=0; 
  opts.spy2=1; 
  opts.logx=0; 
  opts.logy=1; 
  opts.errplot=1;
  opts.phaseplot=0;

  opts.skip_pole=0; 
  opts.skip_res=1;
  opts.cmplx_ss=0;
  opts.legend=0;

  
  if col==1
    %First we fit the colum sum:
    g=0;
    for n=1:Nc
      %g=g+f(n,:);      
      g=g+f(n,:)/norm(f(n,:));
      %g=g+f(n,:)/sqrt(norm(f(n,:)));     
    end
    weight_g=1./abs(g);

    disp('****Improving initial poles by fitting column sum (1st column)...')
    for iter=1:Niter1
       disp(['   Iter ' num2str(iter)])
       if iter==Niter1,opts.skip_res=0; end
      [SER,poles,rmserr,fit]=vectfit3(g,s,poles,weight_g,opts);  
    end
  end


  disp(['****Fitting column #' num2str(col) ' ...'])
  %weight=ones(1,Ns);
  %weight=1./abs(f);
  weight=1./sqrt(abs(f));
  if col==Nc, opts.legend=1; end
  opts.skip_res=1;
  for iter=1:Niter2
    disp(['   Iter ' num2str(iter)])
    if iter==Niter2, opts.skip_res=0; end
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);  
  end

  %Stacking the column contribution into complete state space model:
  bigf=[bigf;f];
  bigfit=[bigfit;fit];
  AA=blkdiag(AA,SER.A);
  BB=blkdiag(BB,SER.B);
  CC=[CC SER.C];
  DD=[DD SER.D];
  EE=[EE SER.E];
  
end %for col=1:Nc

A=AA; B=BB; C=CC; D=DD; E=EE; %Renaming variables


%Finally, we assess the fitting quality (all matrix elements):
freq=s/(2*pi*i);
figure(3),
h1=semilogy(freq,abs(bigf).','b');hold on
h2=semilogy(freq,abs(bigfit).','r--');
h3=semilogy(freq,abs(bigfit-bigf).','g-');hold off
legend([h1(1) h2(1) h3(1)],'Data','FRVF','Deviation');
xlabel('Frequency [Hz]')
title('All matrix elements');

disp('-------------------E N D----------------------------')
toc
