function [SER,poles,rmserr,fit,opts]=vectfit3(f,s,poles,weight,opts);
%
%function [SER,poles,rmserr,fit,opts]=vectfit3(f,s,poles,weight,opts)
%function [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts)
%function [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight)
% 
%     ===========================================================
%     =   Fast Relaxed Vector Fitting                           =
%     =   Version 1.0                                           =
%     =   Last revised: 08.08.2008                              = 
%     =   Written by: Bjorn Gustavsen                           =
%     =   SINTEF Energy Research, N-7465 Trondheim, NORWAY      =
%     =   bjorn.gustavsen@sintef.no                             =
%     =   http://www.energy.sintef.no/Produkt/VECTFIT/index.asp =
%     =   Note: RESTRICTED to NON-COMMERCIAL use                =
%     ===========================================================
%
% PURPOSE : Approximate f(s) with a state-space model 
%
%         f(s)=C*(s*I-A)^(-1)*B +D +s*E
%                       
%           where f(s) is a singe element or a vector of elements.  
%           When f(s) is a vector, all elements become fitted with a common
%           pole set.
%
% INPUT :
%
% f(s) : function (vector) to be fitted. 
%        dimension : (Nc,Ns)  
%                     Nc : number of elements in vector
%                     Ns : number of frequency samples 
% 
% s : vector of frequency points [rad/sec] 
%        dimension : (1,Ns)  
%
% poles : vector of initial poles [rad/sec]
%         dimension : (1,N)  
%
% weight: the rows in the system matrix are weighted using this array. Can be used 
%         for achieving higher accuracy at desired frequency samples. 
%         If no weighting is desired, use unitary weights: weight=ones(1,Ns). 
%
%         Two dimensions are allowed:
%           dimension : (1,Ns) --> Common weighting for all vector elements.   
%           dimension : (Nc,Ns)--> Individual weighting for vector elements.  
%
% opts.relax==1 --> Use relaxed nontriviality constraint 
% opts.relax==0 --> Use nontriviality constraint of "standard" vector fitting
%
% opts.stable=0 --> unstable poles are kept unchanged
% opts.stable=1 --> unstable poles are made stable by 'flipping' them
%                   into the left half-plane
% 
%
% opts.asymp=1 --> Fitting with D=0,  E=0 
% opts.asymp=2 --> Fitting with D~=0, E=0 
% opts.asymp=3 --> Fitting with D~=0, E~=0 
%
% opts.spy1=1 --> Plotting, after pole identification (A)
%              figure(3): magnitude functions
%                cyan trace  : (sigma*f)fit              
%                red trace   : (sigma)fit
%                green trace : f*(sigma)fit - (sigma*f)fit
%
% opts.spy2=1 --> Plotting, after residue identification (C,D,E) 
%              figure(1): magnitude functions
%              figure(2): phase angles  
%
% opts.logx=1 --> Plotting using logarithmic absissa axis             
%
% opts.logy=1 --> Plotting using logarithmic ordinate axis
%
% opts.errplot=1   --> Include deviation in magnitude plot
%
% opts.phaseplot=1 -->Show plot also for phase angle
%
%
% opts.skip_pole=1 --> The pole identification part is skipped, i.e (C,D,E) 
%                    are identified using the initial poles (A) as final poles.
%                 
% opts.skip_res =1 --> The residue identification part is skipped, i.e. only the 
%                    poles (A) are identified while C,D,E are returned as zero. 
%
%
% opts.cmplx_ss  =1 -->The returned state-space model has real and complex conjugate 
%                    parameters. Output variable A is diagonal (and sparse). 
%              =0 -->The returned state-space model has real parameters only.
%                    Output variable A is square with 2x2 blocks (and sparse).
%
% OUTPUT :
% 
%     fit(s) = C*(s*I-(A)^(-1)*B +D +s.*E
%
% SER.A(N,N)    : A-matrix (sparse). If cmplx_ss==1: Diagonal and complex. 
%                                Otherwise, square and real with 2x2 blocks. 
%                           
% SER.B(N,1)    : B-matrix. If cmplx_ss=1: Column of 1's. 
%                       If cmplx_ss=0: contains 0's, 1's and 2's)
% SER.C(Nc,N)   : C-matrix. If cmplx_ss=1: complex
%                       If cmplx_ss=0: real-only
% SERD.D(Nc,1)  : constant term (real). Is non-zero if asymp=2 or 3.
% SERE.E(Nc,1)  : proportional term (real). Is non-zero if asymp=3.
%
% poles(1,N)    : new poles 
%
% rmserr(1) : root-mean-square error of approximation for f(s). 
%                   (0 is returned if skip_res==1)  
% fit(Nc,Ns): Rational approximation at samples. (0 is returned if
%             skip_res==1).
%
%
% APPROACH: 
% The identification is done using the pole relocating method known as Vector Fitting [1], 
% with relaxed non-triviality constraint for faster convergence and smaller fitting errors [2], 
% and utilization of matrix structure for fast solution of the pole identifion step [3]. 
%
%******************************************************************************** 
% NOTE: The use of this program is limited to NON-COMMERCIAL usage only.
% If the program code (or a modified version) is used in a scientific work, 
% then reference should be made to the following:  
%
% [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency       
%     domain responses by Vector Fitting", IEEE Trans. Power Delivery,        
%     vol. 14, no. 3, pp. 1052-1061, July 1999.                                
%
% [2] B. Gustavsen, "Improving the pole relocating properties of vector
%     fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,
%     July 2006.   
%
% [3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
%     "Macromodeling of Multiport Systems Using a Fast Implementation of
%     the Vector Fitting Method", IEEE Microwave and Wireless Components 
%     Letters, vol. 18, no. 6, pp. 383-385, June 2008.
%********************************************************************************
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%
% 
def.relax=1;      %Use vector fitting with relaxed non-triviality constraint
def.stable=1;     %Enforce stable poles
def.asymp=2;      %Include only D in fitting (not E)   
def.skip_pole=0;  %Do NOT skip pole identification
def.skip_res=0;   %Do NOT skip identification of residues (C,D,E) 
def.cmplx_ss=1;   %Create complex state space model
def.spy1=0;       %No plotting for first stage of vector fitting
def.spy2=1;       %Create magnitude plot for fitting of f(s) 
def.logx=1;       %Use logarithmic abscissa axis
def.logy=1;       %Use logarithmic ordinate axis 
def.errplot=1;    %Include deviation in magnitude plot
def.phaseplot=0;  %exclude plot of phase angle (in addition to magnitiude)
def.legend=1;     %Do include legends in plots


if nargin<5
  opts=def;
else 
  %Merge default values into opts  
  A=fieldnames(def);    
  for m=1:length(A)
    if ~isfield(opts,A(m))
      dum=char(A(m)); dum2=getfield(def,dum); opts=setfield(opts,dum,dum2);
    end
  end  
end 



%Tolerances used by relaxed version of vector fitting
TOLlow=1e-18; TOLhigh=1e18;

[a,b]=size(poles); 
if s(1)==0 && a==1
  if poles(1)==0 && poles(2)~=0 
    poles(1)=-1;
  elseif poles(2)==0 && poles(1)~=0 
    poles(2)=-1;
  elseif poles(1)==0 && poles(2)==0 
    poles(1)=-1+i*10; poles(2)=-1-i*10;
  end
end  

if (opts.relax~=0) && (opts.relax)~=1
  disp(['    ERROR in vectfit3.m: ==> Illegal value for opts.relax:  ' num2str(opts.asymp)]),return
end
if (opts.asymp~=1) && (opts.asymp)~=2 && (opts.asymp)~=3
  disp(['    ERROR in vectfit3.m: ==> Illegal value for opts.asymp:  ' num2str(opts.asymp)]),return
end
if (opts.stable~=0) && (opts.stable~=1) 
  disp(['    ERROR in vectfit3.m: ==> Illegal value for opts.stable:  ' num2str(opts.stable)]),return
end
if (opts.skip_pole~=0) && (opts.skip_pole)~=1
  disp(['    ERROR in vectfit3.m: ==> Illegal value for opts.skip_pole:  ' num2str(opts.skip_pole)]),return
end
if (opts.skip_res~=0) && (opts.skip_res)~=1
  disp(['    ERROR in vectfit3.m: ==> Illegal value for opts.skip_res:  ' num2str(opts.skip_res)]),return
end
if (opts.cmplx_ss~=0) && (opts.cmplx_ss)~=1
  disp(['    ERROR in vectfit3.m: ==> Illegal value for opts.cmplx_ss:  ' num2str(opts.cmplx_ss)]),return
end

rmserr=[];%SERC=[];
[a,b]=size(s);
if a<b, s=s.'; end 

% Some sanity checks on dimension of input arrays:
if length(s)~=length(f(1,:))
  disp('Error in vectfit3.m!!! ==> Second dimension of f does not match length of s.'); 
  return;
end  
if length(s)~=length(weight(1,:))
  disp('Error in vectfit3.m!!! ==> Second dimension of weight does not match length of s.');
  return;
end
if length(weight(:,1))~=1
  if length(weight(:,1))~=length(f(:,1))
    disp('Error in vectfit3.m!!! ==> First dimension of weight is neither 1 nor matches first dimension of f.');  
    return;
  end
end  
    
  %set(0,'DefaultLineLineWidth',0.5) ; set(0,'DefaultLineMarkerSize',4) ;
  %clear b; clear C; 
  LAMBD=diag(poles);
  Ns=length(s); N=length(LAMBD); Nc=length(f(:,1));
  B=ones(N,1); %I=diag(ones(1,N));                
  SERA=poles;SERC=zeros(Nc,N);SERD=zeros(Nc,1);SERE=zeros(Nc,1);
  roetter=poles;
  fit=zeros(Nc,Ns);

  weight=weight.';
  if length(weight(1,:))==1
    common_weight=1;
  elseif length(weight(1,:))==Nc
    common_weight=0;
  else
    disp('ERROR in vectfit3.m: Invalid size of array weight')
    return
  end

  if opts.asymp==1
    offs=0; 
  elseif opts.asymp==2
    offs=1;  
  else
    offs=2; 
  end
  
  
  
%=========================================================================
%=========================================================================
%  POLE IDENTIFICATION:
%=========================================================================
%=========================================================================
 
  
if opts.skip_pole~=1  
  


Escale=zeros(1,Nc+1);


%=======================================================
% Finding out which starting poles are complex :
%=======================================================
cindex=zeros(1,N);
for m=1:N 
  if imag(LAMBD(m,m))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 || cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end

%=======================================================
% Building system - matrix :
%=======================================================
 %I3=diag(ones(1,Nc));I3(:,Nc)=[]; 
  Dk=zeros(Ns,N);
  for m=1:N
    if cindex(m)==0      %real pole
      Dk(:,m)=1./(s-LAMBD(m,m));
    elseif cindex(m)==1  %complex pole, 1st part
      Dk(:,m)  =1./(s-LAMBD(m,m)) + 1./(s-LAMBD(m,m)');
      Dk(:,m+1)=i./(s-LAMBD(m,m)) - i./(s-LAMBD(m,m)');
    end
  end      
 if opts.asymp==1 || opts.asymp==2    
   Dk(:,N+1)=1; 
 elseif opts.asymp==3   
   Dk(:,N+1)=1; 
   Dk(:,N+2)=s;
 end
    
%Scaling for last row of LS-problem (pole identification)
scale=0;
for m=1:Nc
  if length(weight(1,:))==1
    scale=scale+(norm(weight(:,1).*f(m,:).'))^2;
  else
    scale=scale+(norm(weight(:,m).*f(m,:).'))^2;
  end  
end
scale=sqrt(scale)/Ns; 

 
if opts.relax==1   

  %Escale=zeros(1,Nc*(N+offs)+N+1);
  %scale=norm(f);%/Ns; %Scaling for sigma in LS problem 
  AA=zeros(Nc*(N+1),N+1);
  bb=zeros(Nc*(N+1),1);
  Escale=zeros(1,length(AA(1,:)));
  for n=1:Nc
    A=zeros(Ns,(N+offs) +N+1); %b=zeros(Ns*Nc+1,1);

    if common_weight==1
      weig=weight;
    else
      weig=weight(:,n);
    end 


    for m=1:N+offs %left block
      A(1:Ns,m)=weig.*Dk(1:Ns,m);
    end
    inda=N+offs;
    for m=1:N+1 %right block
      A(1:Ns,inda+m)=-weig.*Dk(1:Ns,m).*f(n,1:Ns).';    
    end
     
    A=[real(A);imag(A)];
    
   
    %Integral criterion for sigma:
    offset=(N+offs);
    if n==Nc
      for mm=1:N+1  
        A(2*Ns+1,offset+mm)=real(scale*sum(Dk(:,mm)));
      end
    end  
    [Q,R]=qr(A,0);
    ind1=N+offs+1;
    ind2=N+offs+N+1;
    R22=R(ind1:ind2,ind1:ind2);
    AA((n-1)*(N+1)+1:n*(N+1),:)=R22;
    if n==Nc
      bb((n-1)*(N+1)+1:n*(N+1),1)=Q(end,N+offs+1:end)'*Ns*scale; 
    end
  end %for n=1:Nc 

  for col=1:length(AA(1,:))
    Escale(col)=1/norm(AA(:,col));  
    AA(:,col)=Escale(col).*AA(:,col);
  end  
  x=AA\bb;
  %size(x),size(Escale)
  x=x.*Escale.';

end %if opts.relax==0   


%Situation: No relaxation, or produced D of sigma extremely small and large. Solve again, without relaxation 
  if opts.relax==0 | abs(x(end))<TOLlow | abs(x(end))>TOLhigh
    AA=zeros(Nc*(N),N);
    bb=zeros(Nc*(N),1);   
    if opts.relax==0
      Dnew=1; 
    else
      if x(end)==0
        Dnew=1;      
      elseif abs(x(end))<TOLlow
        Dnew=sign(x(end))*TOLlow
      elseif abs(x(end))>TOLhigh 
        Dnew=sign(x(end))*TOLhigh   
      end    
    end

   for n=1:Nc
      A=zeros(Ns,(N+offs) +N); %b=zeros(Ns*Nc+1,1);
      Escale=zeros(1,N);
      
      if common_weight==1
        weig=weight;
      else
        weig=weight(:,n);
      end 

      for m=1:N+offs %left block
        A(1:Ns,m)=weig.*Dk(1:Ns,m);
      end
      inda=N+offs;
      for m=1:N %right block
        A(1:Ns,inda+m)=-weig.*Dk(1:Ns,m).*f(n,1:Ns).';    
      end
      b=Dnew*weig.*f(n,1:Ns).';     
      A=[real(A);imag(A)];
      b=[real(b);imag(b)];       
      offset=(N+offs);
      [Q,R]=qr(A,0);
      ind1=N+offs+1;
      ind2=N+offs+N;
      R22=R(ind1:ind2,ind1:ind2);
      AA((n-1)*N+1:n*N,:)=R22;
      bb((n-1)*N+1:n*N,1)=Q(:,ind1:ind2).'*b;
    end %for n=1:Nc 
    for col=1:length(AA(1,:))
      Escale(col)=1./norm(AA(:,col));  
      AA(:,col)=Escale(col).*AA(:,col);
    end  
    
    %if opts.use_normal==1 
    %  x=AA.'*AA\(AA.'*bb);  
    %else     
      x=AA\bb;
    %end 
    x=x.*Escale.';
    x=[x;Dnew];

  end  %if opts.relax==0 | abs(x(end))<TOLlow | abs(x(end))>TOLhigh
 

  %************************************  

C=x(1:end-1);
D=x(end); %NEW!!

%We now change back to make C complex : 
% **************
for m=1:N
  if cindex(m)==1
    for n=1:1%Nc+1
      r1=C(m); r2=C(m+1);
      C(m)=r1+i*r2; C(m+1)=r1-i*r2;
    end
  end
end
% **************  

if opts.spy1==1
  Dk=zeros(Ns,N); 
  for m=1:N
    Dk(:,m)=1./(s-LAMBD(m,m));
  end 
  RES3(:,1)=D+Dk*C; % (sigma)rat
  freq=s./(2*pi*i); 
  if opts.logx==1  
    if opts.logy==1
      figure(3),
      h1=loglog(freq,abs(RES3'),'b'); xlim([freq(1) freq(Ns)]);  %sigma*f
    else %logy=0
      figure(3),
      h1=semilogx(freq,abs(RES3'),'b'); xlim([freq(1) freq(Ns)]);
    end
  else %logx=0
    if opts.logy==1
      figure(3),
      h1=semilogy(freq,abs(RES3'),'b'); xlim([freq(1) freq(Ns)]);
    else %logy=0
      figure(3),
      h1=plot(s./(2*pi*i),abs(RES3'),'b'); xlim([freq(1) freq(Ns)]);
    end
  end
  figure(3),xlabel('Frequency [Hz]'); ylabel('Magnitude'); 
  %figure(3), title('Sigma')
  if opts.legend==1
    legend([h1(1)],'sigma');
  end  
  drawnow;
end

%=============================================================================
% We now calculate the zeros for sigma :
%=============================================================================
%oldLAMBD=LAMBD;oldB=B;oldC=C;
m=0;
for n=1:N
  m=m+1;
  if m<N  
    if( abs(LAMBD(m,m))>abs(real(LAMBD(m,m))) ) %complex number?
      LAMBD(m+1,m)=-imag(LAMBD(m,m)); LAMBD(m,m+1)=imag(LAMBD(m,m));
      LAMBD(m,m)=real(LAMBD(m,m));LAMBD(m+1,m+1)=LAMBD(m,m);
      B(m,1)=2; B(m+1,1)=0; 
      koko=C(m); C(m)=real(koko); C(m+1)=imag(koko);
      m=m+1;
    end
  end
end

ZER=LAMBD-B*C.'/D;
roetter=eig(ZER).';
unstables=real(roetter)>0;  
if opts.stable==1
  roetter(unstables)=roetter(unstables)-2*real(roetter(unstables)); %Forcing unstable poles to be stable...
end
roetter=sort(roetter); N=length(roetter);


%=============================================
%Sorterer polene s.a. de reelle kommer først:
for n=1:N
  for m=n+1:N
    if imag(roetter(m))==0 && imag(roetter(n))~=0
      trans=roetter(n); roetter(n)=roetter(m); roetter(m)=trans;
    end
  end
end
N1=0;
for m=1:N
  if imag(roetter(m))==0, N1=m; end
end
if N1<N, roetter(N1+1:N)=sort(roetter(N1+1:N)); end   % N1: n.o. real poles
%N2=N-N1;                                             % N2: n.o. imag.poles

roetter=roetter-2*i*imag(roetter); %10.11.97 !!!
SERA=roetter.';


end %if skip_pole~=1  

%=========================================================================
%=========================================================================
%  RESIDUE IDENTIFICATION:
%=========================================================================
%=========================================================================

if opts.skip_res~=1

%=============================================================================
% We now calculate SER for f, using the modified zeros of sigma as new poles :
%========================================================================================

%clear LAMBD A A1 xA1 xxA1 A2 xA2 xxA2 b xb xxb C RES1 RES2;

LAMBD=roetter; 

%B=ones(N,1); 

% Finding out which poles are complex :
cindex=zeros(1,N);
for m=1:N 
  if imag(LAMBD(m))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 || cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end


%===============================================================================
% We now calculate the SER for f (new fitting), using the above calculated 
% zeros as known poles :
%=============================================================================== 
if opts.asymp==1
  A=zeros(2*Ns,N); BB=zeros(2*Ns,Nc);
elseif opts.asymp==2
  A=zeros(2*Ns,N+1); BB=zeros(2*Ns,Nc);
else
  A=zeros(2*Ns,N+2); BB=zeros(2*Ns,Nc);
end

 %I3=diag(ones(1,Nc));I3(:,Nc)=[]; 
  Dk=zeros(Ns,N); 
    for m=1:N
      if cindex(m)==0      %real pole
        Dk(:,m)=1./(s-LAMBD(m));
      elseif cindex(m)==1  %complex pole, 1st part
        Dk(:,m)  =1./(s-LAMBD(m)) + 1./(s-LAMBD(m)');
        Dk(:,m+1)=i./(s-LAMBD(m)) - i./(s-LAMBD(m)');
      end
    end 


if common_weight==1

 %I3=diag(ones(1,Nc));I3(:,Nc)=[]; 
  Dk=zeros(Ns,N); 
    for m=1:N
      if cindex(m)==0      %real pole
        Dk(:,m)=weight./(s-LAMBD(m));
      elseif cindex(m)==1  %complex pole, 1st part
        Dk(:,m)  =weight./(s-LAMBD(m)) + weight./(s-LAMBD(m)');
        Dk(:,m+1)=i.*weight./(s-LAMBD(m)) - i.*weight./(s-LAMBD(m)');
      end
    end 

  if opts.asymp==1
    A(1:Ns,1:N)=Dk; 
  elseif opts.asymp==2
    A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=weight;
  else
    A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=weight; A(1:Ns,N+2)=weight.*s;
  end
  for m=1:Nc
    BB(1:Ns,m)=weight.*f(m,:).';
  end
  A(Ns+1:2*Ns,:)=imag(A(1:Ns,:));
  A(1:Ns,:)=real(A(1:Ns,:));
  BB(Ns+1:2*Ns,:)=imag(BB(1:Ns,:));
  BB(1:Ns,:)=real(BB(1:Ns,:));

  if opts.asymp==2
    A(1:Ns,N+1)=A(1:Ns,N+1);             
  elseif opts.asymp==3
    A(1:Ns,N+1)=A(1:Ns,N+1);             
    A(Ns+1:2*Ns,N+2)=A(Ns+1:2*Ns,N+2);  
  end

  %clear Escale;
  Escale=zeros(1,length(A(1,:)));
  for col=1:length(A(1,:));
    Escale(col)=norm(A(:,col),2);
    A(:,col)=A(:,col)./Escale(col);
  end
  X=A\BB;
  for n=1:Nc
    X(:,n)=X(:,n)./Escale.';
  end

  %clear A;
  X=X.';
  C=X(:,1:N); 
  if opts.asymp==2
    SERD=X(:,N+1);
  elseif opts.asymp==3
    SERE=X(:,N+2); 
    SERD=X(:,N+1);
  end

else %if common_weight==1
    
  SERD=zeros(Nc,1); 
  SERE=zeros(Nc,1); 
  C=zeros(Nc,N);
  for n=1:Nc

    if opts.asymp==1
      A(1:Ns,1:N)=Dk; 
    elseif opts.asymp==2
      A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=1;
    else
      A(1:Ns,1:N)=Dk; A(1:Ns,N+1)=1; A(1:Ns,N+2)=s;
    end

    for m=1:length(A(1,:))
      A(1:Ns,m)=weight(:,n).*A(1:Ns,m);
    end
 
    BB=weight(:,n).*f(n,:).';
    A(Ns+1:2*Ns,:)=imag(A(1:Ns,:));
    A(1:Ns,:)=real(A(1:Ns,:));
    BB(Ns+1:2*Ns)=imag(BB(1:Ns));
    BB(1:Ns)=real(BB(1:Ns));

    if opts.asymp==2
      A(1:Ns,N+1)=A(1:Ns,N+1);             
    elseif opts.asymp==3
      A(1:Ns,N+1)=A(1:Ns,N+1);             
      A(Ns+1:2*Ns,N+2)=A(Ns+1:2*Ns,N+2);  
    end

    %clear Escale;
    Escale=zeros(1,length(A(1,:)));
    for col=1:length(A(1,:));
      Escale(col)=norm(A(:,col),2);
      A(:,col)=A(:,col)./Escale(col);
    end
    x=A\BB;
    x=x./Escale.';

    %clear A;
    C(n,1:N)=x(1:N).'; 

    if opts.asymp==2
      SERD(n)=x(N+1);
    elseif opts.asymp==3
      SERE(n)=x(N+2); 
      SERD(n)=x(N+1);
    end


  end %for n=1:Nc

    
end %if common_weight==1    
  

%=========================================================================

%We now change back to make C complex. 
for m=1:N
  if cindex(m)==1
    for n=1:Nc
      r1=C(n,m); r2=C(n,m+1);
      C(n,m)=r1+i*r2; C(n,m+1)=r1-i*r2;
   end
  end
end
% **************  


B=ones(N,1);

%====================================================

SERA  = LAMBD;
SERB  = B;
SERC  = C;



  Dk=zeros(Ns,N); 
  for m=1:N
    Dk(:,m)=1./(s-SERA(m));
  end 
  for n=1:Nc
    fit(n,:)=(Dk*SERC(n,:).').';
    if opts.asymp==2
      fit(n,:)=fit(n,:)+SERD(n);
    elseif opts.asymp==3
      fit(n,:)=fit(n,:)+SERD(n)+s.'.*SERE(n);
    end
  end 
  
  fit=fit.'; f=f.';
  diff=fit-f; rmserr=sqrt(sum(sum(abs(diff.^2))))/sqrt(Nc*Ns);

if opts.spy2==1
  freq=s./(2*pi*i);
  if opts.logx==1     
    if opts.logy==1
      figure(1),
      h1=loglog(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=loglog(freq,abs(fit),'r--'); hold off
      if opts.errplot== 1
        hold on,h3=loglog(freq,abs(f-fit),'g');hold off; 
      end
    else %logy=0 
      figure(1),
      h1=semilogx(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogx(freq,abs(fit),'r--'); hold off
      if opts.errplot== 1
        hold on,h3=semilogx(freq,abs(f-fit),'g');hold off; 
      end
    end
    if opts.phaseplot==1
      figure(2),
      h4=semilogx(freq,180*unwrap(angle(f))/pi,'b'); xlim([freq(1) freq(Ns)]);hold on
      h5=semilogx(freq,180*unwrap(angle(fit))/pi,'r--');hold off
    end  
  else %logx=0
    if opts.logy==1
      figure(1),
      h1=semilogy(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=semilogy(freq,abs(fit),'r--'); hold off
      if opts.errplot== 1
        hold on,h3=semilogy(freq,abs(f-fit),'g');hold off; 
      end
    else %logy=0 
      figure(1), 
      h1=plot(freq,abs(f),'b'); xlim([freq(1) freq(Ns)]);hold on
      h2=plot(freq,abs(fit),'r--'); hold off
      if opts.errplot== 1
        hold on,h3=plot(freq,abs(f-fit),'g');hold off; 
      end
    end
    if opts.phaseplot==1    
      figure(2),
      h4=plot(freq,180*unwrap(angle(f))/pi,'b'); xlim([freq(1) freq(Ns)]);hold on
      h5=plot(freq,180*unwrap(angle(fit))/pi,'r--');hold off
    end  
  end %logy=0
  figure(1),
  xlabel('Frequency [Hz]'); ylabel('Magnitude');
  %title('Approximation of f');
  if opts.legend==1
    if opts.errplot==1
      legend([h1(1) h2(1) h3(1)],'Data','FRVF','Deviation');
    else
      legend([h1(1) h2(1)],'Data','FRVF');       
    end  
  end
  if opts.phaseplot==1
    figure(2),
    xlabel('Frequency [Hz]'); ylabel('Phase angle [deg]');
    %title('Approximation of f');
    if opts.legend==1
      legend([h4(1) h5(1)],'Data','FRVF'); 
    end
  end
  drawnow;
end
fit=fit.';


end %if skip_res~=1

A=SERA;
poles=A;
if opts.skip_res~=1
  B=SERB; C=SERC; D=SERD; E=SERE;
else
  B=ones(N,1); C=zeros(Nc,N); D=zeros(Nc,Nc); E=zeros(Nc,Nc); rmserr=0;  
end;    


%========================================
% Convert into real state-space model
%========================================
if opts.cmplx_ss~=1 

  A=diag(sparse(A));
  
  cindex=zeros(1,N);
  for m=1:N 
    if imag(A(m,m))~=0  
      if m==1 
        cindex(m)=1;
      else
        if cindex(m-1)==0 || cindex(m-1)==2
          cindex(m)=1; cindex(m+1)=2; 
        else
          cindex(m)=2;
        end
      end 
    end
  end
  
  n=0;
  for m=1:N
    n=n+1;
    if cindex(m)==1
      a=A(n,n); a1=real(a); a2=imag(a);
      c=C(:,n); c1=real(c); c2=imag(c);
      b=B(n,:); b1=2*real(b); b2=-2*imag(b);
      Ablock=[a1 a2;-a2 a1];
   
      A(n:n+1,n:n+1)=Ablock;
      C(:,n)=c1;
      C(:,n+1)=c2;
      B(n,:)=b1;
      B(n+1,:)=b2;
    end
  end

else
  A=sparse(diag(A));    % A is complex, make it diagonal
end %if cmplx_ss~=1    
    
SER.A=A; SER.B=B; SER.C=C; SER.D=D; SER.E=E;


