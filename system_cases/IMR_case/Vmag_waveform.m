data_sel = 1;
switch data_sel
    case 1; load out_68_good.mat; data = out_68_good;
    case 2; load out_68_bad.mat; data = out_68_bad;
    case 3; load out_68_verybad.mat; data = out_68_verybad;
end

fs=6e3;
ts=1/fs;
tstart=5.5;
tend = 7;

v28= data.Vmag28_39{1}.Values.Data;
v39 = data.Vmag28_39{2}.Values.Data;
v59 = data.Vmag28_39{3}.Values.Data;

drange = (tstart*fs : tend*fs-1);

t = (tstart:ts:tend-ts);
figure(3)
clf;
plot(t, v28(drange),"LineWidth",1.2);
hold on;
plot(t, v39(drange),"LineWidth",1.2);
hold on;
plot(t, v59(drange),"LineWidth",1.2);
ylim([0.8 1.1])
xlabel('time (s)')

legend('Vmag28','Vmag39','Vmag59')
%title('Bus voltage')
grid on;

%%
FigSize = [0.1 0.1 0.5 0.6]*0.5;
set(gcf,'units','normalized','outerposition',FigSize);