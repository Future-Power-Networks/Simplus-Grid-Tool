

%% Back up for plotting
% % Plot admittance
% if Enable_PlotAdmittance
%     fprintf('Plotting admittance spectrum...\n')
%     Tj = [1 1j;     % real to complex transform
%           1 -1j];  
%     for k = 1:N_Bus
%         if DeviceType{k} <= 50
%             Gr_ss{k} = GminSS(Port_i([2*k-1,2*k]),Port_v([2*k-1,2*k]));
%             Gr_sym{k} = SimplusGT.ss2sym(Gr_ss{k});
%             Gr_c{k} = Tj*Gr_sym{k}*Tj^(-1);
%         end
%     end
%  	figure_n = figure_n+1;
%  	figure(figure_n);
%     CountLegend = 0;
%     VecLegend = {};
%     for k = 1:N_Bus
%         if DeviceType{k} <= 50
%             SimplusGT.bode_c(Gr_c{k}(1,1),1j*omega_pn,2*pi,'PhaseOn',0); 
%             CountLegend = CountLegend + 1;
%             VecLegend{CountLegend} = ['Bus',num2str(k)];
%         end
%     end
%     legend(VecLegend);
%     SimplusGT.mtit('Bode diagram: admittance');
% else
%     fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
% end

% % Plot w related
% if Enable_PlotSwing
%     fprintf('Plotting frequency-port dynamics...\n')
%     for k = 1:N_Bus
%         if floor(DeviceType{k}/10) == 0
%             Gt_ss{k} = GminSS(Port_w(k),Port_T_m(k));
%             Gt_sym{k} = -SimplusGT.ss2sym(Gt_ss{k});  % The negative sign is because of the load convention.
%         elseif floor(DeviceType{k}/10) == 1
%          	Gt_ss{k} = GminSS(Port_w(k),Port_ang_r(k));
%             Gt_sym{k} = -SimplusGT.ss2sym(Gt_ss{k});
%         end
%     end
%  	figure_n = figure_n+1;
%  	figure(figure_n);
%     CountLegend = 0;
%     VecLegend = {};
%     for k = 1:N_Bus
%         if (floor(DeviceType{k}/10) == 0) || (floor(DeviceType{k}/10) == 1)
%             SimplusGT.bode_c(Gt_sym{k},1j*omega_pn,2*pi,'PhaseOn',0);      
%          	CountLegend = CountLegend + 1;
%             VecLegend{CountLegend} = ['Bus',num2str(k)]; 
%         end
%     end
%     legend(VecLegend);
%     SimplusGT.mtit('Bode diagram: frequency-port transfer function');
%     
% else
%     fprintf('Warning: The default plot of omega-related spectrum is disabled.\n');
% end