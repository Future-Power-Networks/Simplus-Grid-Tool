% Author(s): Yitong Li

function ModalAnalysisStateSpace(ObjGsysSs,EigIndex,FigN)
    
    % Get the model and string    
    [~,GsysSs] = ObjGsysSs.GetSS(ObjGsysSs); 
    [GsysSsStateStr,~,~] = ObjGsysSs.GetString(ObjGsysSs);

    % Calculate the eigenvalue and eigenvector
    [PhiMat,EigMat] = eig(GsysSs.A);
    PsiMat = PhiMat^(-1);
    EigVec = diag(EigMat);
    EigVecHz = EigVec/2/pi;
    fprintf('The indices of eigenvalues whose real parts are larger than -10Hz: \n')
    IndexRiskEig = find( real(EigVecHz)>-10);
    SimplusGT.PrintArray(IndexRiskEig,10);
    
    % Initialize the participation factor matrix
    PfMat = zeros(length(GsysSs.A),0);
    for i = 1:length(EigVecHz)
        PfMatLambdai = PhiMat(:,i)*PsiMat(i,:);
        PfMatLambdai = diag(PfMatLambdai);
        PfMat = [PfMat,PfMatLambdai];
    end
    fprintf('Participation factor matrix: PfMat(k,i) \n')
    fprintf('where i -> ith eigenvalue, k -> kth state or a_kk in state matrix A. \n')
    
    for i =1:length(EigIndex)
        if EigIndex(i)<= 0
            error(['Error: The mode index ' num2str(EigIndex(i)) ' is wrong, which must be a positive integer.'])
        end
        if EigIndex(i) > length(EigVecHz)
            fprintf(['Warning: The mode index ' num2str(EigIndex(i)) ' is larger than the total number of eigenvalues ' num2str(length(EigVecHz)) ', please double check ModeIndex in UserMain.m \n'])
            EigIndex(i) = 1;
        end
        TargetMode(i) = EigVecHz(EigIndex(i));
        PfVec = PfMat(:,EigIndex(i));
        PfVec = abs(PfVec);
        PfVecSum = sum(PfVec);
        for PfMaxNum = 1:length(PfVec)
            [PfMaxValue,PfMaxIndex] = maxk(PfVec,PfMaxNum);
            PfRestIndex = setdiff([1:length(PfVec)],PfMaxIndex);
            PfRestValue = sum(PfVec(PfRestIndex));
            PfAllValue = [PfMaxValue;PfRestValue];
            
            % The rest of the states that participates less than certain
            % value.
            if (PfRestValue/PfVecSum<=0.1) || (PfMaxNum>=10) 
                break;
            end
        end
        fprintf('\n')
        fprintf(['The ' num2str(EigIndex(i)) 'th eigenvalue (' num2str(TargetMode(i)) ' Hz) is investigated...\n'])
        fprintf(['Top ' num2str(PfMaxNum) ' largest participation factor(s): \n'])
        
        for k = 1:length(PfMaxIndex)
            fprintf([num2str(PfMat(PfMaxIndex(k),EigIndex(i))) ' -> ' num2str(PfMaxIndex(k)) 'th state '  GsysSsStateStr{PfMaxIndex(k)} '\n']);
            PfStateStr{k} = strrep(GsysSsStateStr{PfMaxIndex(k)},'_','\_');     % For print in pie chart
        end
        
        % Plot pie chart
        if 1
            figure(FigN)
            PfStateStr{k+1} = 'Others';
            pie(PfAllValue,PfStateStr);
            title(['Participation of the ' num2str(EigIndex(i)) 'th Eigenvalue (' num2str(TargetMode(i)) ' Hz)'])
            FigN = FigN+1;
        end
        clear('PfStateStr')
    end
end