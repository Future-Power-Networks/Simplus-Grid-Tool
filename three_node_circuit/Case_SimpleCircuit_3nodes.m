% this is the 3-node simple circuit model in the paper. It ueses symbolic
% calculations as it's a simple circuit.
% Modes (eigenvalues) can be found from the variable 'Lambda'.
% S_lambda_c: critical admittance-eigenvalue sensitivity matrix;
% The missing coefficient xi can be found in 'xi', matched with modes.
% S_lambda_m1: eigenvalue sensitivity matrix calculated from equation (16);
% S_lambda_m2: eigenvalue sensitivity matrix calculated from equation (20);
% S_lambda_m1 and S_lambda_m2 are equal.
% Results of Table-I can be found from 'Sen_c' and 'Sen', matched with
% Lambda. 
% For results of Table-II, run 'plot_sweep.m' then go to
% 'P_Result' to check the results.

clear;
syms   s   R1   L1   C1  R2  L2 C2 R3 L3 C3 R12 L12 R13 L13 R23 L23
param =   [R1    L1  C1  R2  L2 C2 R3 L3 C3 R12 L12 R13 L13 R23 L23];
p_value = [1.2  5.8  7.9 2.2 2  4  5  5  6  0.5 0.3 0.6 0.2 1.5 0.8];
y1 = 1/(R1+s*L1)+s*C1;
y2 = 1/(R2+s*L2)+s*C2;
y3 = 1/(R3+s*L3)+s*C3;
y12 = 1/(R12+s*L12);
y13 = 1/(R13+s*L13);
y23 = 1/(R23+s*L23);

Yn = [y1+y12+y13    -y12        -y13;
      -y12          y2+y12+y23  -y23;
      -y13          -y23        y3+y13+y23;];

Yn_det_sym = det(Yn); %symbolic determinent
Yn_det_sym_s = subs(Yn_det_sym,param,p_value); % substitute the values

[num,den]=numden(Yn_det_sym_s);
Num=sym2poly(num);
Den=sym2poly(den);
Yn_det_tf = tf(Num,Den);
Lambda = zero(Yn_det_tf);
nl=length(Lambda);

%% Left-Right Critical Eigenvalue ----Layer 1
S_lambda_c=cell(nl,1);
for i = 1: nl
    lambda = Lambda(i);
    Ysys = eval(subs(Yn,...
    [s param],[lambda p_value])); % substitute all the values to the nodal matrix);
    [V,D]=eig(Ysys);
    W= inv(V);
    crit_num = find(abs(diag(D))<=1e-10); % locate the critical eigen-value
    if length(crit_num)>=2 || isempty(crit_num)
        error('critical eigenvalue not found');
    end
    S_lambda_c{i} = V(:,crit_num) * W(crit_num,:);    
end


%% xi --- the missing coefficient calculation
xi = zeros(nl,1);
S_lambda_m1 = cell(nl,1);
Ydet_p=diff(Yn_det_sym_s, s);  % calculate the partial differential Y'det(s)
for i = 1: nl
    lambda = Lambda(i);
    Ydet_p_lambda = eval(vpa(subs(Ydet_p,s,lambda)));%subsitute s=lambda to Y'det(lambda)
    Ysys = eval(vpa(subs(Yn,[s param],[lambda p_value]))); % substitute all the values to the nodal matrix);
    xi(i) = -trace(adjoint(Ysys))/Ydet_p_lambda;
    S_lambda_m1{i} = xi(i) * S_lambda_c{i};
end

%% residue method
Zsys_s = inv(subs(Yn,param,p_value));
S_lambda_m2 = cell(nl,1);
for i = 1:nl
    lambda = Lambda(i); %row
    for j = 1:length(Zsys_s) %column
        for k = 1:length(Zsys_s)
            [num,den]=numden(Zsys_s(j,k));
            Num=sym2poly(num);
            Den=sym2poly(den);
            Zsys_s_tf = tf(Num,Den);
            [R,P,K] = residue(Zsys_s_tf.num{1},Zsys_s_tf.den{1});
            r_num=find(P==lambda);
            S_lambda_m2{i}(j,k)=-1*R(r_num);
        end
    end
end

%% compare the results
for i = 1:nl
    if max(abs(S_lambda_m1{i}(:)-S_lambda_m2{i}(:))) >= 1e-10
        error('%d not equal', i)
    end
end

%% sensitivity with respect to shunt and series components
Sen=cell(nl,1);
Sen_c=cell(nl,1);
%Sen_layer1=cell(nl,1);
for i = 1:nl
    for j = 1:length(Zsys_s) %column
        for k = 1:length(Zsys_s)
            if(j==k)
                Sen{i}(j,k) = S_lambda_m1{i}(j,j);
                Sen_c{i}(j,k)= S_lambda_c{i}(j,j);
            else
                Sen{i}(j,k)= S_lambda_m1{i}(j,j) + S_lambda_m1{i}(k,k) -...
                    S_lambda_m1{i}(j,k) - S_lambda_m1{i}(k,j);
                Sen_c{i}(j,k)= S_lambda_c{i}(j,j) + S_lambda_c{i}(k,k) -...
                    S_lambda_c{i}(j,k) - S_lambda_c{i}(k,j);
            end
        end
    end
    Sen_layer1(i).EigVal = Lambda(i);
    Sen_layer1(i).xi = xi(i);
    for n=1:length(Yn(:))
        Sen_layer1(i).Result(n).Sen_c = Sen_c{i}(n);
        Sen_layer1(i).Result(n).Sen_c_abs_norm = abs(Sen_c{i}(n))*abs(eval(vpa(subs(Yn(n),[s param],[Lambda(i) p_value]))));
        Sen_layer1(i).Result(n).xi = xi(i);
        Sen_layer1(i).Result(n).Sen_sys = Sen{i}(n);
        Sen_layer1(i).Result(n).Sen_sys = Sen{i}(n);
        Sen_layer1(i).Result(n).Sen_sys_abs_norm = abs(Sen{i}(n))*abs(eval(vpa(subs(Yn(n),[s param],[Lambda(i) p_value]))));
    end
end
%% Layer3 calculation
%Layer3 = cell(length(Lambda),1);
eig_sel = 1:length(Lambda);
for sel=1:length(eig_sel)
    eigk = eig_sel(sel);
    lambda=Lambda(eigk);
    Layer3(sel).lambda = lambda;
    for j=1:length(param)
        rho = param(j);
        pYnprho = eval(vpa(subs(diff(Yn,rho),[s param],[lambda p_value])));
        Layer3(sel).Result(j).param_name = char(param(j));
        Layer3(sel).Result(j).PlamdaPrho=trace(S_lambda_m1{eigk} * pYnprho); %calculate p_lambda/p_rho.
        Layer3(sel).Result(j).PlamdaPrho_pu = Layer3(sel).Result(j).PlamdaPrho * p_value(j);
        a=Layer3(sel).Result(j).PlamdaPrho_pu;
        Layer3(sel).Result(j).PlamdaPrho_norm = sqrt(a*conj(a));
    end
end

run plot_sweep.m
