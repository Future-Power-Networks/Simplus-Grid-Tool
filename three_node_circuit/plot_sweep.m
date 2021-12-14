%% plot
%sw_eig_sel = [2,2,2];
%[R1  L1  C1  R2  L2  C2  R3  L3  C3  R12 L12 R13 L13 R23 L23];
sw_rho_sel = [6 10 11 12 9 1 2]; 
sw_eig_sel = [2 2  2  6  6 8 8];
pert_mag = [1.05 1.05 1.05 1.05 1.05 1.05 1.05];
point_num=15;


p_value_sw=p_value;
figure(1);
clf;
plot(Lambda,'x');
axis([-2 0 -1 1]);
hold on; 

colormap lines
newcolors = colormap;
color_n = length(newcolors);
clear dLambda;
clear P_Result;
for k=1:length(sw_rho_sel)
    param_k=sw_rho_sel(k);
    p_value_sw = p_value;
    lambda_k=sw_eig_sel(k);
    Lambda_sw1=Lambda;
    rho_sweep = linspace(1*p_value(param_k),pert_mag(k)*p_value(param_k),point_num);
    for rho = rho_sweep
        p_value_sw(param_k)=rho;
        Yn_det_sym_s_sw = subs(Yn_det_sym,param,p_value_sw); % substitute the values
        [num,den]=numden(Yn_det_sym_s_sw);
        Num=sym2poly(num);
        Den=sym2poly(den);
        Yn_det_tf_sw = tf(Num,Den);
        Lambda_sw = zero(Yn_det_tf_sw);      
        %re-order the eigen-values.
        d0=inf;
        Lambda_sw_= Lambda_sw;
        for j=1:length(Lambda)
            for n=1:length(Lambda)
                d1=abs(Lambda_sw(n)-Lambda_sw1(j));
                if d1<d0
                    Lambda_sw_(j)=Lambda_sw(n);
                    d0=d1;
                end
            end
        end
        Lambda_sw = Lambda_sw_;
        Lambda_sw1=Lambda_sw; 
        %plot
        p=plot([Lambda_sw(lambda_k),Lambda_sw(lambda_k+1)],'*');
        %p=plot(Lambda_sw,'*');
        p.Color = newcolors(k,:);
        hold on;       
    end
    
%     for i=1:length(Lambda)
%         dLambda(k).sweep_name = char(param(param_k));
%         dLambda(k).result(i).orig_lambda = Lambda(i);
%         dLambda(k).result(i).final_lambda = Lambda_sw(i);
%         dLambda(k).result(i).delta_lambda = Lambda_sw(i)-Lambda(i);
%         %dLambda(k).result(i).error = 
%     end
    P_Result(k).eig = Lambda(lambda_k);
    P_Result(k).param = char(param(param_k));
    P_Result(k).PlambdaPrho_pu = Layer3(lambda_k).Result(param_k).PlamdaPrho_pu;
    P_Result(k).perturbation = pert_mag(k)-1;
    P_Result(k).delta_lambda = Lambda_sw(lambda_k)-Lambda(lambda_k);
    a = P_Result(k).PlambdaPrho_pu * P_Result(k).perturbation; % expected value
    ar = P_Result(k).delta_lambda; % real value
    P_Result(k).err_real = (real(ar)-real(a))/real(a);
    P_Result(k).err_imag = (imag(ar)-imag(a))/imag(a);
    P_Result(k).error = sqrt((ar-a)*conj(ar-a)) / sqrt(a*conj(a));
    P_Result(k).error2 = norm((ar-a),2) / norm(a,2);
end
hold off;
grid on;

%axis([-0.16 -0.04 -0.06 0.06 ])

%C_result(1)
% Lambda_sw(8)-Lambda(8);
%set(gca,'ylim',[-1.5,1.5]);
%set(gca,'xlim',[-1.5,1.5]);