function plotc(Xw,wbd,varargin)

    LineWidth = LoadVar(1.5,'LineWidth',varargin);
    LineStyle = LoadVar('-','LineStyle',varargin);
    Color = LoadVar([],'Color',varargin);
    PhaseOn = LoadVar(1,'PhaseOn',varargin);

    [M,N,W] = size(Xw);

    if wbd(1)*wbd(W) < 0
        % seperate positive and negative frequency
        wbdn = wbd(1:W/2);
        wbdp = wbd(W/2+1:W);
        Xwn  = Xw(:,:,1:W/2);
        Xwp  = Xw(:,:,W/2+1:W);

        % anti wind up
        Arg_wn = angle(Xwn);
        Arg_wp = angle(Xwp);
        for m = 1:M
            for n = 1:N
                Arg_wn(m,n,:) = flip(Arg_wn(m,n,:));
                for k = 1:(W/2-1)
                    while (Arg_wn(m,n,k+1) - Arg_wn(m,n,k) > 1.5*pi)
                        Arg_wn(m,n,k+1) = Arg_wn(m,n,k+1) - 2*pi;
                    end
                    while (Arg_wn(m,n,k+1) - Arg_wn(m,n,k) < -1.5*pi)
                        Arg_wn(m,n,k+1) = Arg_wn(m,n,k+1) + 2*pi;
                    end
                end
                Arg_wn(m,n,:) = flip(Arg_wn(m,n,:));

                for k = 1:(W/2-1)
                    while (Arg_wp(m,n,k+1) - Arg_wp(m,n,k) > 1.5*pi)
                        Arg_wp(m,n,k+1) = Arg_wp(m,n,k+1) - 2*pi;
                    end
                    while (Arg_wp(m,n,k+1) - Arg_wp(m,n,k) < -1.5*pi)
                        Arg_wp(m,n,k+1) = Arg_wp(m,n,k+1) + 2*pi;
                    end
                end
            end
        end

        if PhaseOn == 0
            for m = 1:M
                for n = 1:N
                    if M*N > 1
                        figure();
                    end
                    subplot(1,2,1);
                    p(1)= loglog(wbdn,abs(squeeze(Xwn(m,n,:))));
                    grid on;  hold on;

                    subplot(1,2,2);
                    p(2)= loglog(wbdp,abs(squeeze(Xwp(m,n,:))));
                    grid on;  hold on;
                end
            end
        else
            for m = 1:M
                for n = 1:N
                    if M*N > 1
                        figure();
                    end
                    subplot(2,2,1);
                    p(1)= loglog(wbdn,abs(squeeze(Xwn(m,n,:))));
                    grid on;  hold on;

                    subplot(2,2,3);
                    p(2)= semilogx(wbdn,squeeze(Arg_wn(m,n,:))*180/pi);
                    grid on;  hold on;

                    subplot(2,2,2);
                    p(3)= loglog(wbdp,abs(squeeze(Xwp(m,n,:))));
                    grid on;  hold on;

                    subplot(2,2,4);
                    p(4)= semilogx(wbdp,squeeze(Arg_wp(m,n,:))*180/pi);
                    grid on;  hold on;
                end
            end    
        end
    else

        % anti wind up
        Arg_w = angle(Xw);
        for m = 1:M
            for n = 1:N
                for k = 1:(W-1)
                    while (Arg_w(m,n,k+1) - Arg_w(m,n,k) > 1.5*pi)
                        Arg_w(m,n,k+1) = Arg_w(m,n,k+1) - 2*pi;
                    end
                    while (Arg_w(m,n,k+1) - Arg_w(m,n,k) < -1.5*pi)
                        Arg_w(m,n,k+1) = Arg_w(m,n,k+1) + 2*pi;
                    end
                end
            end
        end

        if PhaseOn == 0
            for m = 1:M
                for n = 1:N
                    if M*N > 1
                        figure();
                    end
                    p = loglog(wbd,abs(squeeze(Xw(m,n,:))));
                    grid on;  hold on; 
                end
            end 
        else
            for m = 1:M
                for n = 1:N
                    if M*N > 1
                        figure();
                    end
                    subplot(2,1,1);
                    p(1)= loglog(wbd,abs(squeeze(Xw(m,n,:))));
                    grid on;  hold on;

                    subplot(2,1,2);
                    p(2)= semilogx(wbd,squeeze(Arg_w(m,n,:))*180/pi);
                    grid on;  hold on;  
                end
            end    
        end
    end
    
    try 
        p; %#ok<VUNUS>
        for h = 1:length(p)
            p(h).LineWidth = LineWidth;
            p(h).LineStyle = LineStyle;
            if ~isempty(Color)
                p(h).Color = Color;
            end
        end
    catch
    end

end