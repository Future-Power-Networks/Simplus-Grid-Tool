% This function links the device models

% Author(s): Yitong Li, Yunjie Gu

function Gobj = DeviceModelLink(GmObj)

    % Create a new obj
    Gobj = SimplexPS.Class.ModelBase;
    Gobj.SetDSS(Gobj,dss([],[],[],[],[]));
    Gobj.SetString(Gobj,{},{},{});
    
    % Append
    for n = 1:length(GmObj)
        Gobj = SimplexPS.ObjAppend(Gobj,GmObj{n});
    end

    %% Re-arrange input and output
%     % This part can be removed if using the strings to locate the required
%     % ports
%     
%     % Get strings
%     [StateStr,InputStr,OutputStr] = Gobj.ReadString(Gobj);
%     
%     % Find vdq and get index vector
%     CountIn1 = 0;
%     CountIn2 = 0;
%     for n = 1:length(InputStr)
%         if ( strcmp(InputStr{n},'v_d') || strcmp(InputStr{n},'v_q') ...
%           || strcmp(InputStr{n},'vd') || strcmp(InputStr{n},'vq') )
%             CountIn1 = CountIn1+1;
%             IndexIn1(CountIn1) = n;
%         else
%             CountIn2 = CountIn2+1;
%             IndexIn2(CountIn2) = n;
%         end
%     end
%     IndexIn = [IndexIn2,IndexIn1];
%     
%     % Find idq and get index vector
%     CountOut1 = 0;
%     CountOut2 = 0;
%  	for n = 1:length(OutputStr)
%         if ( strcmp(OutputStr{n},'i_d') || strcmp(OutputStr{n},'i_q') ...
%           || strcmp(OutputStr{n},'id') || strcmp(OutputStr{n},'iq') )
%             CountOut1 = CountOut1+1;
%             IndexOut1(CountOut1) = n;
%         else
%             CountOut2 = CountOut2+1;
%             IndexOut2(CountOut2) = n;
%         end
%     end
%     IndexOut = [IndexOut2,IndexOut1];
%     
%     % Re-arrange input and output for matrices
%     Gobj.ModelDSS.B = Gobj.ModelDSS.B(:,IndexIn);
%     Gobj.ModelDSS.C = Gobj.ModelDSS.C(IndexOut,:);
%     Gobj.ModelDSS.D = Gobj.ModelDSS.D(IndexOut,IndexIn);
%     
%     % Re-arrange input and output for string
%     InputStr = InputStr(:,IndexIn);
%     OutputStr = OutputStr(:,IndexOut);
%     Gobj.WriteString(Gobj,StateStr,InputStr,OutputStr);
end