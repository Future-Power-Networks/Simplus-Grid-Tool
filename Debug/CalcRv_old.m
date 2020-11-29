     Ak = ModelSS.A;
     Ck = ModelSS.C;
     Bk = ModelSS.B;
     Wk = inv(eye(length(Ak)) - Ts/2*Ak);
     MatrixY = Ts/2*Ck*Wk*Bk;
     MatrixY = MatrixY(1:2,1:2);
     MatrixR = inv(MatrixY);
     DeviceDiscreDampingResistor = MatrixR(1,1);