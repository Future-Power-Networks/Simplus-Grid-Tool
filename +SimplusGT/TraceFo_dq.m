% Trace of two dq-based matrix.
% Out = Trace(a,b), where the elements is a d-q struct.
% down into each element, the product is calculated as <a(i,j),b(i,j)>fro
% Author: Yue Zhu

function TraceFinal = TraceFo_dq(a,b)

if size(a)~=size(b)
    error('size must be the same')
end
[row_num,col_num]=size(a);
if row_num~=col_num
    error('both must be square matrix')
end
TraceFinal =0;

for k=1:length(a)
    for j=1:length(a)
        TraceFinal= TraceFinal + SimplusGT.inner_product_dq(a(k,j),b(k,j));
    end
end

end