function [R,a]=ss2pr(A,B,C)
%
% Convert state-space model having COMMON POLE SET into pole-residue model.
%
% Input:
% A,B,C: must have the format produced by vectfit2.m. Both
% formats determined by parameter VF.cmplx_ss are valid
% 
% Output:
% R(Nc,Nc,N) %Residues
% a(N)       %poles 
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%

%Converting real-only state-space model into complex model, if necessary:  
if max(max(abs(A-diag(diag(A)))))~=0
  errflag=0;
  for m=1:length(A)-1
    if A(m,m+1)~=0
      A(m,m)    =A(m,m)+i*A(m,m+1); 
      A(m+1,m+1)=A(m+1,m+1)-i*A(m,m+1);  
      
      B(m,:)  =(B(m,:)+B(m+1,:))/2;
      B(m+1,:)=B(m,:); 
      
      C(:,m)  =C(:,m)+i*C(:,m+1);
      C(:,m+1)=conj(C(:,m));  
    end
  end
end  

%Converting complex state-space model into pole-residue model
Nc=length(C(:,1));
N=length(A)/Nc;
R=zeros(Nc,Nc,N);
for m=1:N
  Rdum=zeros(Nc);
  for n=1:Nc
    ind=(n-1)*N+m;
    Rdum=Rdum +C(:,ind)*B(ind,:);
  end  
  R(:,:,m)=Rdum;
end
a=full(diag(A(1:N,1:N)));