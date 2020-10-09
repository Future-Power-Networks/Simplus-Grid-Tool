%% 
% =======================================================================
% Following codes are discarded currently, and should not be activated.
% =======================================================================
% if 0
% 
% %% Add device load into model
% Size_Load = [Size_Device(1),Size_Bus(2)];
% Shift_Load = Shift_Device;
% 
% Size_LoadGND = Size_DeviceGND;
% Shift_LoadGND = [-50,20];
% 
% % Add load device
% for i = 1:N_Device
%     if floor(DeviceType{i}/10) == 9 % This might be changed in the future
%         % Calculate R and L
%         switch DeviceType{i}
%             case 90
%                 P = PowerFlow{i}(1);
%                 Q = PowerFlow{i}(2);
%                 V = PowerFlow{i}(3);
%                 if Para{i}.Connection == 1
%                     S = P + 1i*Q;
%                     I = conj(S/V);
%                     Z = V/I;
%                     R = real(Z);
%                     L = imag(Z)/W0;
%                 elseif Para{i}.Connection == 2
%                     R = V^2/P;
%                     L = V^2/Q/W0;
%                 else
%                     error(['Error']);
%                 end
%             case 91
%                 R = Para{i}.R;
%                 L = Para{i}.L;
%             otherwise
%                 error(['Error']);
%         end
%         
%         % Add block
%         switch Para{i}.Connection
%             case 1
%                 Name_Load{i} = ['Series-Load' num2str(i)];
%                 FullName_Load{i} = [Name_Model '/' Name_Load{i}];
%                 add_block(['powerlib/Elements/Three-Phase Series RLC Branch'],FullName_Load{i});
%             case 2
%                 Name_Load{i} = ['Parallel-Load' num2str(i)];
%                 FullName_Load{i} = [Name_Model '/' Name_Load{i}];
%                 add_block(['powerlib/Elements/Three-Phase Parallel RLC Branch'],FullName_Load{i});
%             otherwise
%                 error(['Error']);
%                 
%         end
%         
%         % Set detailed load type
%         if R == 0 
%             set_param(gcb,'BranchType','L');
%         elseif L == 0
%             set_param(gcb,'BranchType','R');
%         else
%             set_param(gcb,'BranchType','RL');
%         end
%         
%         % Set customer data
%         set_param(gcb,'Resistance',num2str(R));
%         set_param(gcb,'Inductance',num2str(L));
%         
%       	% The position of device is set by referring to the position of correpsonding bus
%         Position_Load{i} = Position_Bus{i} + Shift_Load;
%         set_param(gcb,'position',[Position_Load{i},Position_Load{i}+Size_Load]);
%         
%       	% Connect load to bus
%      	add_line(Name_Model,...
%             {[Name_Load{i} '/Rconn1'],[Name_Load{i} '/Rconn2'],[Name_Load{i} '/Rconn3']},...
%             {[Name_Bus{i} '/Lconn1'],[Name_Bus{i} '/Lconn2'],[Name_Bus{i} '/Lconn3']},...
%             'autorouting','smart');
%         
%         % Connect the floating side of load to ground
%      	Name_LoadGND{i} = ['L-GND' num2str(i)];
%         FullName_LoadGND{i} = [Name_Model '/' Name_LoadGND{i}];
%         add_block('powerlib/Elements/Ground',FullName_LoadGND{i});
%         PortPosition_Load{i} = get_param(FullName_Load{i},'PortConnectivity');
%         Position_LoadGND{i} = PortPosition_Load{i}(2).Position;
%         Position_LoadGND{i} = Position_LoadGND{i} + Shift_LoadGND;
%         set_param(FullName_LoadGND{i},'position',[Position_LoadGND{i},Position_LoadGND{i}+Size_LoadGND]);
%         
%         % Connect load to ground
%         add_line(Name_Model,[Name_Load{i} '/LConn2'],[Name_LoadGND{i} '/LConn1'],...
%                 'autorouting','smart');
%         
%         % Form the Y confiration of the floating terminal
%         add_line(Name_Model,...
%             {[Name_Load{i} '/Lconn2'],[Name_Load{i} '/Lconn2']},...
%             {[Name_Load{i} '/Lconn1'],[Name_Load{i} '/Lconn3']});
%     end
% end
% 
% end