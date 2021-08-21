% bode plot for complex transfer function
% plot max singular value for transfer function matrix
% ### modified 2019.10.23: use matlabFunction to speed up symbolic eval
% ### modified 2019.10.25: anti wind up for phase plot
% ### modified 2019.10.27: add transfer function matrix plot (not just singular value)

function Xw = bodec(X,sbd,wbase,varargin)

    Option = LoadVar(0,'Option',varargin);
    PhaseOn = LoadVar(1,'PhaseOn',varargin);
    
    if (Option == 1)
        PlotOn = 1;
    elseif (length(X) == 1)
        PlotOn = 1;
    else
        PlotOn = 0;
    end
    PlotOn = LoadVar(PlotOn,'PlotOn',varargin);

    if (Option == 1)
        PhaseOn = 0;    %no phase for singular value
    end

    [M,N] = size(X);

    %X = simplify(X);
    %X = simplifyFraction(X);
    funcX = matlabFunction(X);

    if Option == 0           % bode for transfer function matrix
        Xw = zeros(M,N,length(sbd));
        for n = 1:length(sbd)
            try 
                Xw(:,:,n) = funcX(sbd(n));
            catch
                Xw(:,:,n) = funcX();
            end
            if(isnan(Xw(:,:,n)))
                Xw(:,:,n) = zeros(M,N);
            elseif(isinf(Xw(n)))
                Xw(:,:,n) = zeros(M,N);
            end
        end
    elseif Option == 1       % singular value for matrix transfer function
        Xw = zeros(1,1,length(sbd));
        for n = 1:length(sbd)
            try 
                Xw(1,1,n) = max(svd(funcX(sbd(n)))); 
            catch
                Xw(1,1,n) = max(svd(funcX()));
            end
            if(isnan(Xw(n)))
                Xw(1,1,n) = 0;
            elseif(isinf(Xw(n)))
                Xw(1,1,n) = 0;
            end
        end
    else
        return;
    end

    if PlotOn == 1
        plotc(Xw,imag(sbd)/wbase,'PhaseOn',PhaseOn,varargin);
    end
    
end
    
    


