% transform a state space model to complex frame

% Author(s): Yunjie Gu

function Gc = ss2cplx(Gr,ne)
    nu = length(Gr.B(1,:));
    ny = length(Gr.C(:,1));
    
    mu = nu - ne*2;
    my = ny - ne*2;
    
    Tj = [1 1j;1 -1j];  % real to complex
    Ti = Tj^(-1);       % complex to real

    TjCell = repmat({Tj}, 1, ne);
    Tjblk = blkdiag(TjCell{:});
    
    TiCell = repmat({Ti}, 1, ne);
    Tiblk = blkdiag(TiCell{:});

    Py = blkdiag(eye(my),Tjblk);
    Pu = blkdiag(eye(mu),Tiblk);
    
    Gc.A = Gr.A;
    Gc.B = Gr.B*Pu;
    Gc.C = Py*Gr.C;
    Gc.D = Py*Gr.D*Pu;
    
end