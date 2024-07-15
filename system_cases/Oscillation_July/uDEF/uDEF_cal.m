%save('Oscillation_July/uDEF/out_14bus_case1.mat','out');
% DEF_case=2;
Napp=6;
lengendstr = [{'A1'},{'A2'},{'A3'},{'A6'},{'A8'},{'A11'},{'A12'},{'A13'}];
% if DEF_case==2
%     load('Oscillation_July/uDEF/out_14bus_case2.mat');
% elseif DEF_case==1
%     load('Oscillation_July/uDEF/out_14bus_case1.mat');
% elseif DEF_case== 9
%     load('Oscillation_July/uDEF/out_7bus_case_udef.mat');
%     Napp=6;
%     lengendstr = [{'A1'},{'A3'},{'A4'},{'A5'},{'A6'},{'A7'}];
% end

%Fs=2.5e4;
new_rate = 1e3;
figure(11);
figure(12);
for i=1:8
    W=zeros(2,4);
    col=mod(i-1,4)+1;
    if i>=5
        W(2,col)=1;
    else
        W(1,col)=1;
    end
    W
    DL=data_length*Fs;
    km=Napp; % 8 machines  
    uDEF_val=zeros(new_rate*data_length,km);
    
    for kk=1:km
       FieldName = strcat('pmu',num2str(kk));
       idq=timeseries2timetable(out.(FieldName){1}.Values);
       idq = retime(idq, 'regular','mean','TimeStep', seconds(1/new_rate));
       vdq=timeseries2timetable(out.(FieldName){2}.Values);
       vdq = retime(vdq, 'regular','mean','TimeStep', seconds(1/new_rate));
        for k=1:length(vdq.Time)
            if k>2
                uDEF_val(k,kk)=idq.Data(k,:) * W * vdq.Data(k,:).'+uDEF_val(k-1,kk);
            else
                uDEF_val(k,kk)=idq.Data(k,:) * W * vdq.Data(k,:).';
            end
        end
    end
    uDEF_TT = timetable(idq.Time,uDEF_val);
    figure(12);
    subplot(2,4,i);
    plot(uDEF_TT.Time, uDEF_TT.Var1);
    %legend(lengendstr);
    figure(11);
    subplot(2,4,i);
    bar(uDEF_TT.Var1(end,:));
end
