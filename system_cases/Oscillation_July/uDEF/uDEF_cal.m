%save('Oscillation_July/uDEF/out_14bus_case1.mat','out');
uDEF_parameters
DEF_case=9;
if DEF_case==2
    load('Oscillation_July/uDEF/out_14bus_case2.mat');
    lengendstr = {'A1','A2','A3','A6','A8','A11','A12','A13'};
    Napp=8;
elseif DEF_case==1
    load('Oscillation_July/uDEF/out_14bus_case1.mat');
    lengendstr = {'A1','A2','A3','A6','A8','A11','A12','A13'};
    Napp=8;
elseif DEF_case== 9
    load('Oscillation_July/uDEF/out_7bus_case_udef.mat');
    Napp=6;
    lengendstr = {'A1','A3','A4','A5','A6','A7'};
end

%Fs=2.5e4;
new_rate = 1e3;

%vdq=timeseries2timetable(out.pmu7{2}.Values);
vdq=timeseries2timetable(out.pmu6{2}.Values);
vdq = retime(vdq, 'regular','mean','TimeStep', seconds(1/new_rate));
fig=figure(33); clf;
fig.Position = [600 300 1000 700];
plot(vdq.Time,vdq.Data(:,1:2))
plot(vdq.Time(27000:end),vdq.Data(27000:end,1:2))
legend({'vd12','vq12'})
title('dVd and dVq: A12')
fig=figure(11);clf;fig.Position =[431,327.6666666666666,990.6666666666665,429.3333333333334];
fig=figure(12);clf;fig.Position =[431,327.6666666666666,990.6666666666665,429.3333333333334];
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
    if DEF_case==1
        plot(uDEF_TT.Time(20000:end), uDEF_TT.Var1(20000:end,:));
    else
        plot(uDEF_TT.Time, uDEF_TT.Var1);
    end
    %legend(lengendstr);
    figure(11);
    subplot(2,4,i);
    %bar(uDEF_TT.Var1(end,:));

    data = uDEF_TT.Var1(end,:);
    data = data/(sum(abs(data(:))));
    X = categorical(lengendstr);
    X = reordercats(X,lengendstr);
    b=bar(X, data);
    
    bar_num=length(data);
    cmap=lines(bar_num);
    b.FaceColor='flat';
    for kk=1:bar_num
        b.CData(kk,:)=cmap(kk,:);
    end
    set(gca,'YLim',[-1,1]);
    %set(gca,'XTickLabel','(a)');
end
