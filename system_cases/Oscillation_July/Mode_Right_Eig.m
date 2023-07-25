
clear Phi_vec_trim

A=GminSS.A;
[Phi,D]=eig(A);
Psi=inv(Phi); 

%ModeSel = 151; % case1
%ModeSel = 202; % case2

Phi_vec=Phi(:,ModeSel);

Psi_vec=Psi(ModeSel,:).';

StateString=GminStateStr;

i=1; sp=1;
for k =1: N_Apparatus
    if ApparatusType{k} <= 89  %apparatus
        %epsilon_x = Phi_vec(sp,1);%*D(ModeSel, ModeSel);
        %ex=real(epsilon_x);
        %ei=imag(epsilon_x);
        %id = Phi_vec(sp+1,1);
        %iq = Phi_vec(sp+2,1);
        %id_p = (cos(ex)*id-sin(ex)*iq);
        %iq_p = sin(ex)*id + cos(ex)*iq;
        %Phi_vec_trim(i:i+1) = exp(-ei)*[id_p,iq_p];
        Phi_vec_trim(i:i+1) = Phi_vec(sp+1:sp+2);
        Psi_vec_trim(i:i+1) = Psi_vec(sp+1:sp+2);
        sp = sp+length(ApparatusStateStr{k});
        i = i+2;      
    else %floating bus and passive load: not considered           
    end
end

Phi_vec_trim = Phi_vec_trim.'; % non-conjugate transpose.
Psi_vec_trim = Psi_vec_trim.';

Mode_Rad = D(ModeSel, ModeSel);
Mode_Hz = Mode_Rad/2/pi;
freq_sel = imag(Mode_Hz);
Phi_vec_trim(:,2) = abs(Phi_vec_trim(:,1));
Phi_vec_trim(:,3) = angle(Phi_vec_trim(:,1));