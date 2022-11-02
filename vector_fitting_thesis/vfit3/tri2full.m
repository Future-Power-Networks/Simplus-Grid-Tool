function [SER2]=tri2full(SER);
%
% Convert rational model of lower matrix triangle into state-space model of
% full matrix.
%
% Input:
% SER: Output structure from vectfit2 when fitting lower triangle of a square matrix
% using a common pole set.
% Both formats determined by parameter VF.cmplx_ss are valid
% 
% Output:
% SER: State-space model of full matrix (with common pole set)
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%

A=SER.A; B=SER.B; C=SER.C; D=SER.D; E=SER.E;

tell=0;
for k=1:1e4
  tell=tell+k;
  if tell==length(D), Nc=k; break, end
end

N=length(A);
tell=0;
CC=zeros(Nc,Nc*N);
AA=[]; BB=[];
for col=1:Nc
  AA=blkdiag(AA,A);
  BB=blkdiag(BB,B);  
  for row=col:Nc
    tell=tell+1;
    DD(row,col)=D(tell);
    EE(row,col)=E(tell);
    CC(row,(col-1)*N+1:col*N)=C(tell,:); 
    CC(col,(row-1)*N+1:row*N)=C(tell,:);     
  end
end
DD=DD+(DD-diag(diag(DD))).';
EE=EE+(EE-diag(diag(EE))).';

SER2.A=AA; SER2.B=BB; SER2.C=CC; SER2.D=DD; SER2.E=EE;