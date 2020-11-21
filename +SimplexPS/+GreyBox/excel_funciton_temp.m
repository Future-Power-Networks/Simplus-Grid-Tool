%prefix 'Gb' represents Greybox.
%get all the modes
ModeNum=length(GminSS.A);
GbMode=zeros(ModeNum,1);    % eigenvalue
[GbV,GbLambda]=eig(GminSS.A);
GbW=inv(GbV);
for i=1:ModeNum
    GbMode(i)=GbLambda(i,i)/(2*pi);
end

%Gbtemp= GbW * (GsysDSS.A) * GbV;
      
% for i=1:N_Device
% xlswrite('GreyBoxConfig.xlsx',DeviceStateStr{i},'state','B');
% end