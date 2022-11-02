% Author: Yue

function [SER2]=col2full(SER)

A=SER.A; B=SER.B; C=SER.C; D=SER.D; E=SER.E;

%tell=0;
%for k=1:1e4
%  tell=tell+k;
%  if tell==length(D), Nc=k; break, end
%end

Nc=2;
N=length(A);
tell=0;
CC=zeros(Nc,Nc*N);
AA=[]; BB=[];
for col=1:Nc
  AA=blkdiag(AA,A);
  BB=blkdiag(BB,B);  
  for row=1:Nc  
    tell=tell+1;
    DD(row,col)=D(tell);
    EE(row,col)=E(tell);
    CC(row,(col-1)*N+1:col*N)=C(tell,:); 
    CC(col,(row-1)*N+1:row*N)=C(tell,:);
  end
end

SER2.A=AA; SER2.B=BB; SER2.C=CC; SER2.D=DD; SER2.E=EE;
end