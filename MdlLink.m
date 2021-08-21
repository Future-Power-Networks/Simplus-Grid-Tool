function Gm = MdlLink(varargin)

    if iscell(varargin{1})
        arg = varargin{1};
    else
        arg = varargin;
    end
    
    Gm = arg{1};
    for n = 2:length(arg)
        Gm = append(Gm,arg{n});
    end

    m = 0;    
    seq1 = [];
    seq2 = [];
    for n = 1:length(arg)
        k = length(arg{n}.B(1,:));
        seq1 = [seq1, (m + (1:(k-2)))]; %#ok<*AGROW>
        seq2 = [seq2, (m + ((k-1):k))];
        m = m + k;
    end

    pvect = [seq1,seq2];

    Gm.B = Gm.B(:,pvect);
    Gm.C = Gm.C(pvect,:);
    Gm.D = Gm.D(pvect,pvect);

end