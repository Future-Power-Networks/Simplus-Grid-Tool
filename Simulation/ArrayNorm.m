function Output = ArrayNorm(Input)
    [r,c] = size(Input);
    for i = 1:r
        Output(i) = norm(Input(i,:));
    end
end