if DEF_case==3 % 1.5 Hz case
    fc_lpf=0.3; %Hz;
    k_lpf = 0.707; %damping ratio
    fpll = 3; % pll bandwidth
    Tdef_start=70;
    data_length=30; %save x~s data
    wc_lpf = 2*pi*fc_lpf;
    fc_lpf2=0.3;
    [lpf_A,lpf_B,lpf_C,lpf_D]=tf2ss(1,[1/(2*pi*0.3) 1]);
elseif DEF_case==2 % 18 Hz case
    fc_lpf=1; %Hz;
    k_lpf = 0.707; %damping ratio
    fpll = 3; % pll bandwidth
    Tdef_start=10;
    data_length=10; %save x~s data
    wc_lpf = 2*pi*fc_lpf;
    fc_lpf2=1;
    [lpf_A,lpf_B,lpf_C,lpf_D]=tf2ss(1,[1/(2*pi*fc_lpf2) 1]);
end