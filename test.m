clear; clc;
close all;
cmap = get(groot,'defaultAxesColorOrder');

%# common variables
w0 = 2*pi*50;
s = sym('s');
I = eye(2);
w = logspace(-2,3,1e5)*2*pi;
wd = [-flip(w),w];

%%
%# default parameters
para1.J = 3.5*2/w0^2;
para1.D = 1/w0^2;             
para1.L = 0.02/w0;
para1.R = 0.01;

para2.V_dc = 2.5;
para2.C_dc = 2*0.1*para2.V_dc^2;
para2.kp_v_dc = para2.V_dc*para2.C_dc*(10*2*pi);
para2.ki_v_dc = para2.kp_v_dc*(10*2*pi)/4;
para2.kp_pll = 2*2*pi;
para2.ki_pll = para2.kp_pll * (2*2*pi)/4; 
para2.tau_pll = 1/(2*pi*200);
para2.k_pf = 0;
para2.L = 0.03/(w0);
para2.R = 0.01;
para2.kp_i_dq = para2.L * (500*2*pi);
para2.ki_i_dq = para2.kp_i_dq *(500*2*pi)/4;

%%
%# Test: low inertia interact with pll
    
layout = 4;
sweep = 2;

if layout == 1
    %2 buses : test the match between symbolic and tf models     
    %-----------------------------------------------------------
    %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
    Bus     = [ 1     1     1.0     0     1.0    0    0.0    0     -1     1;
                2     3     1.0     0     1.0  -0.2   0.0    0     -1     1];
    %-----------------------------------------------------------
    %         |  From |  To   |   R   |   L   |   C   |   G   |
    %         |  Bus  |  Bus  |       |       |       |       |
    Line    = [  1       2      0.00     0.3    0.00      inf;
                 1       1        0       0     1e-1      2.0;
                 2       2        0       0     1e-1      0.0];     
elseif layout == 2
    %3 buses : test generator interaction        
    %-----------------------------------------------------------
    %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
    Bus     = [ 1     1     1.0     0     0.5    0    0.0    0     -1     1;
                2     2     1.0     0     0.5    0    0.0    0     -1     1;
                3     2     1.0     0     1.0    0    0.0    0     -1     1];            
    %-----------------------------------------------------------
    %         |  From |  To   |   R   |   L   |   C   |   G   |
    %         |  Bus  |  Bus  |       |       |       |       |
    Line    = [  1       2      0.01     0.3      0      inf;
                 2       3      0.01     0.3      0      inf;
                 3       1      0.01     0.3      0      inf;
                 1       1        0       0     2e-2     1.0;
                 2       2        0       0     2e-2     1.0;
                 3       3        0       0     2e-2     0.0]; 
elseif layout == 3      
    %4 buses : test generator converter interaction  
    %-----------------------------------------------------------
    %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
    Bus     = [ 1     2     1.0     0     0.5    0    0.0    0     -1     1;
                2     2     1.0     0     0.5    0    0.0    0     -1     1;
                3     1     1.0     0     0.5    0    0.0    0     -1     1;
                4     3     1.0     0     0.5  -0.2   0.0    0     -1     1];
    %-----------------------------------------------------------
    %         |  From |  To   |   R   |  wL   |  wC   |   G   |
    %         |  Bus  |  Bus  |       |       |       |       |
    Line    = [  1       2      0.01     0.3      0      inf;
                 2       3      0.01     0.3      0      inf;
                 3       1      0.01     0.3      0      inf;
                 3       4      0.01     0.3      0      inf;
                 1       1        0       0     1e-5     0.6;
                 2       2        0       0     1e-5     0.6;
                 3       3        0       0     1e-5     0.8;
                 4       4        0       0     1e-5      0];  
elseif layout == 4      
    %5 buses : test generator converter interaction in more meshed grid  
    %-----------------------------------------------------------
    %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
    Bus     = [ 1     2     1.0     0     0.5    0    0.0    0     -1     1;
                2     2     1.0     0     0.5    0    0.0    0     -1     1;
                3     1     1.0     0     0.5    0    0.0    0     -1     1;
                4     3     1.0     0     0.4  -0.2   0.0    0     -1     1;
                5     3     1.0     0     0.1   0.0   0.0    0     -1     1];
    %-----------------------------------------------------------
    %         |  From |  To   |   R   |  wL   |  wC   |   G   |
    %         |  Bus  |  Bus  |       |       |       |       |
    Line    = [  1       2      0.01     0.3      0      inf;
                 2       3      0.01     0.3      0      inf;
                 3       1      0.01     0.3      0      inf;
                 3       4      0.01     0.3      0      inf;
                 3       5      0.01     0.3      0      inf;
                 1       4      0.10     1.0      0      inf;
                 2       5      0.10     1.0      0      inf;
                 1       1        0       0     1e-5     0.6;
                 2       2        0       0     1e-5     0.6;
                 3       3        0       0     1e-5     0.8;
                 4       4        0       0     1e-5     0.0;
                 5       5        0       0     1e-5     0.0];
end

[~,~,Ang0,P0,Q0,V0]=PowerFlow(Bus,Line);

nbus = max(Bus(:,1));

[Ytf1,Yb1] = YbusCalcTF(Line(1:(end-nbus),:),w0);
[Ytf2,Yb2] = YbusCalcTF(Line((end-nbus+1):end,:),w0);
Yb1 = minreal(Yb1);
Zb2 = inv(Yb2);
Zb2 = minreal(Zb2);
Zb = feedback(Zb2,Yb1);

if layout == 1
    bandwidth_sweep = linspace(1,10,10);
elseif layout == 2
    bandwidth_sweep = linspace(1,10,10);
elseif layout == 3
    %bandwidth_sweep = linspace(5,20,10);
    bandwidth_sweep = 20;
elseif layout == 4
    bandwidth_sweep = 20;    
end

%bandwidth_sweep = 20;%debug

cpoint = 1; %colormap pointer
for bandwidth = bandwidth_sweep

    % reduce generator inertia
    para1_ = para1;
    para1_.J = para1_.J/10;
    para1__= para1;
    para1__.J = para1__.J/2;

    % change converter control gain
    para2_ = para2;
    if sweep == 1        %sweep pll gain        
        para2_.kp_pll = bandwidth*2*pi;
        para2_.ki_pll = para2_.kp_pll * (5*2*pi)/4; 
        para2_.kp_v_dc = para2_.V_dc*para2_.C_dc * (10*2*pi);
        para2_.ki_v_dc = para2_.kp_v_dc * (10*2*pi)/4;
    elseif sweep == 2    %sweep pll and dc-link gain
        para2_.kp_pll = bandwidth*2*pi;
        para2_.ki_pll = para2_.kp_pll * (bandwidth*2*pi)/4;                      
        para2_.kp_v_dc = para2_.V_dc*para2_.C_dc * (bandwidth*2*pi);
        para2_.ki_v_dc = para2_.kp_v_dc * (bandwidth*2*pi)/4;
    end

    % parameters
    if layout == 1
        type = {0,10};
        para = {para1_,para2_};
    elseif layout == 2
        % three generators
        type = {0,0,0};
        para = {para1,para1,para1};
    elseif layout == 3
        % four generators
        %type = {0,0,0,0};
        %para = {para1,para1,para1_,para1};
        % three generators and one wind farm
        type = {0,0,0,10};
        para = {para1,para1,para1_,para2_};
    elseif layout == 4
        type = {0,0,0,10,10};
        para = {para1,para1__,para1_,para2_,para2};
    end

    % model and link
    Gm = cell(1,nbus);
    Gc = cell(1,nbus);
    for n = 1:nbus
        [~,Gm{n},Gc{n},~] = MdlCreate('type', type{n} ,'flow',[-P0(n) -Q0(n) V0(n) Ang0(n) w0],'para',para{n});
    end
    Gm = MdlLink(Gm);
    dimu = length(Gm.B(1,:));
    Gsys = feedback(Gm,Zb,(dimu-2*nbus+1):dimu,(dimu-2*nbus+1):dimu);

    % manual link for two-node system
    if (layout == 1) && 1
        Ys1 = Gc{1}((end-1):end,(end-1):end);
        Ys2 = Gc{2}((end-1):end,(end-1):end);
        Ys1 = Ys1 + eye(2)*Line(2,6);
        Ys1 = Ys1 + [(1j + 1/w0*s) 0;0 (-1j + 1/w0*s)]*Line(2,5);
        Ys2 = Ys2 + [(1j + 1/w0*s) 0;0 (-1j + 1/w0*s)]*Line(3,5);                       
        Zs1 = Ys1^(-1);
        Zs1 = Zs1 + [(1j + 1/w0*s) 0;0 (-1j + 1/w0*s)]*Line(1,4);
        Ys1 = Zs1^(-1);
        Y = Ys1 + Ys2;
        Z = Y^(-1);
        pc = vpasolve(1/Z(1,1))/2/pi;
    end

    if 1    % pole plots
        psys = pole(Gsys)/2/pi;
        figure(layout);
        scatter(real(psys),imag(psys),'x','LineWidth',1.5);
        hold on; grid on;
        if (layout == 1) && 1
            scatter(real(pc),imag(pc),'o','LineWidth',1.5);            
        end
        xlabel('Real Part (Hz)');
        ylabel('Imaginary Part (Hz)');
        axis([-25,5,-80,80]);
    end
            
    if 1    % torque coefficient bode plots            
        for n = 1:length(type)
            %Gtw{n} = -tf2sym(tf(Gsys(2*n-1,2*n-1)));      %#ok<SAGROW>
            Gtw{n} = -ss2sym(Gsys(2*n-1,2*n-1));          %#ok<SAGROW>
            if type{n} < 10
                Htw{n} = 1/(para{n}.J*s);                 %#ok<SAGROW>
            elseif type{n} < 20
                Htw{n} = 1/(1+para{n}.tau_pll*s)*(para{n}.kp_pll + para{n}.ki_pll/s); %#ok<SAGROW>
            end           
            Kwt{n} = vpa(Gtw{n}^(-1) - Htw{n}^(-1));      %#ok<SAGROW>
            Gtt{n} = Kwt{n}* vpa(Htw{n});                 %#ok<SAGROW>
        end
   
        disp(['### test' num2str(cpoint) ': bandwidth=' num2str(bandwidth) ' ###']);
        for nplot = 1:length(type)
            figure(layout+10*nplot);
            if layout <= 3
                %isstable is very slow for layout >= 4, unsolved !!!
                if ~isstable(Kwt{nplot})
                    disp(['K' num2str(nplot) ' is not stable']);
                end
            end
            bodec(Gtt{nplot}*1j,1j*w,2*pi,'Color',cmap(mod(cpoint-1,length(cmap))+1,:));                    
            hold on;
            grid on;                           
        end
    end   
    
    if 1    %impedance bode plots
        Tj = [1 1j;1 -1j];  % real to complex
        for n = 1:length(type)
            Ytr{n}(1,1) = ss2sym(Gsys(length(type)*2 +2*n-1 ,length(type)*2 +2*n-1));  %#ok<SAGROW>
            Ytr{n}(1,2) = ss2sym(Gsys(length(type)*2 +2*n-1 ,length(type)*2 +2*n  ));  %#ok<SAGROW>
            Ytr{n}(2,1) = ss2sym(Gsys(length(type)*2 +2*n   ,length(type)*2 +2*n-1));  %#ok<SAGROW>
            Ytr{n}(2,2) = ss2sym(Gsys(length(type)*2 +2*n   ,length(type)*2 +2*n  ));  %#ok<SAGROW>
            Ytc{n} = Tj*Ytr{n}*Tj^(-1);  %#ok<SAGROW>
        end
   
        for nplot = 1:length(type)
            figure(layout+100);
            bodec(Ytc{nplot}(1,1),1j*wd,2*pi,'PhaseOn',0);                    
            hold on;
            grid on;                           
        end
    end
    
    cpoint = cpoint+1;
end   

if layout == 4
    fig = figure(layout+100);
    axis([1e0,1e3,1e-2,1e3]);
    fig.Position = [185 148 720 310];
    lgd = legend('Generator1','Generator2','Generator3','WindFarm1','WindFarm2');
    lgd.FontSize = 10;
    print(gcf,'Y.png','-dpng','-r600');
elseif layout == 1
    fig = figure(layout+10);
    fig.Position = [195 -28 560 590];
    fig.Children(1).FontSize = 8;
    fig.Children(2).FontSize = 8;
    print(gcf,'K.png','-dpng','-r600');
end
