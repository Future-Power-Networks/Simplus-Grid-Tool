Zpole_rad=pole(Zsys_restore);
Zpole_Hz=Zpole_rad/2/pi;



%% Residue Calculation
% step 1: find the mode index
mode_Hz=41.9; %freq of the selected mode
pole_imag_err=abs(imag(Zpole_Hz)-mode_Hz);
mk=find(pole_imag_err<0.1); % index of the mode
if length(mk)>1
    error('more than one mode is spoted');
end
if isempty(mk)
    error('mode not found');
end

% step 2: calculate the residue matrix
Rmat_rea = zeros(2*N_Bus,2*N_Bus);
[num,den] = tfdata(Zsys_restore);
for i=1:length(Rmat_rea)
    for j=1:length(Rmat_rea)       
        [r,p,~] = residue(num{i,j},den{i,j});
        if abs(p(mk)/2/pi-Zpole_Hz(mk))>1e-5
            error('mode mismatched!')
        else
            Rmat_rea(i,j)=r(mk);
        end
    end
end

% step 3: calculate the apparatus admittance YAr
YAr_est=evalfr(YAr,Zpole_rad(mk));

% step 4: participation analysis layer 1 and layer 2
clear restore_pf_results;
for k = 1:N_Bus
    
    if ApparatusType{k} ~= 100 %not a floating bus)

        Res_ = Rmat_rea(2*k-1:2*k,2*k-1:2*k);
        Ya = YAr_est(2*k-1:2*k,2*k-1:2*k);
        restore_pf_results.layer1(k)=sqrt(trace(Res_*Res_')) * sqrt(trace(Ya*Ya'));
        %restore_pf_results.layer1(k)= norm(-Rmat_rea(2*k-1:2*k,2*k-1:2*k)') * norm(YAr_est(2*k-1:2*k,2*k-1:2*k));
        restore_pf_results.layer2(k)= trace((-Rmat_rea(2*k-1:2*k,2*k-1:2*k))* YAr_est(2*k-1:2*k,2*k-1:2*k));
    else % floating bus, infinite bus...
        restore_pf_results.layer1(k)=0;
        restore_pf_results.layer2(k)=0;
    end
end

%restore_pf_results.layer1
%restore_pf_results.layer2

xx=figure(11); % pie chart for layer1
explode=ones(1,length(restore_pf_results.layer1));
p=pie(restore_pf_results.layer1);
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
for i=1:length(restore_pf_results.layer1)
    if str2double(erase(percentValues{i},'%')) >0.1 
        pText(i).String = strcat('A',num2str(i), '(',percentValues(i), ')');
    else   %for value smaller than 0.1%, don't show a string.
        pText(i).String = '';
    end
end
cmap=lines(length(restore_pf_results.layer1));
xx.Colormap=cmap;

% xx=figure(12);
% X = categorical(legend);
% X = reordercats(X,legend);
% b=bar(ax,X, data);
% 
% bar_num=length(data);
% cmap=lines(bar_num);
% b.FaceColor='flat';
% for i=1:bar_num
%     b.CData(i,:)=cmap(i,:);
% end