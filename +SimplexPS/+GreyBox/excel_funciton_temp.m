%prefix 'Gb' represents Greybox.
%get all the modes
ModeNum=length(GsysTemp.A);
StateNum=ModeNum;
GbMode=zeros(ModeNum,1);    % eigenvalue
%[GbV,GbLambda,GbW]=eig(GsysDSS.A,GsysDSS.E); % solve generalized eigenvalues
%GbW=inv(GsysDSS.E*GbV); %Left eigenvector matrix.
Gsys=dss2ss(GsysDSS);
[GbV,GbLambda]=eig(Gsys.A);
for i=1:ModeNum
    GbMode(i)=GbLambda(i,i)/(2*pi);
end

temp = GbW*GsysDSS.A*GbV;

ModeSelect = 33; % select the 33rd mode for further analysis.

% for i=1:ModeNum
%     for k=1:StateNum
%         PfState(i,k) = GbW(i,k) * GbV(k,i);
%     end
% end
% 
% temp=0;
% for i=1:StateNum
%     temp = temp + PfState(ModeSelect,i);
% end

%Gbtemp= GbW * (GsysDSS.A) * GbV;
      
% for i=1:N_Device
% xlswrite('GreyBoxConfig.xlsx',DeviceStateStr{i},'state','B');
% end