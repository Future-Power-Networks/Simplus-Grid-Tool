% Author: Yue Zhu

function yout = dqStrutMinus(a,b)

if size(a)~=size(b)
    error('size must be the same')
end
[row_num,col_num]=size(a);
if row_num~=col_num
    error('both must be square matrix')
end

n = length(a);
a_exp=zeros(2*n);
b_exp=zeros(2*n);
for k=1:n
    for j=1:n
        a_exp(2*k-1,2*j-1)=a(k,j).dd;
        a_exp(2*k-1, 2*j )=a(k,j).dq;
        a_exp(2*k, 2*j-1 )=a(k,j).qd;
        a_exp(2*k ,2*j   )=a(k,j).qq;
        
        b_exp(2*k-1,2*j-1)=b(k,j).dd;
        b_exp(2*k-1, 2*j )=b(k,j).dq;
        b_exp(2*k, 2*j-1 )=b(k,j).qd;
        b_exp(2*k ,2*j   )=b(k,j).qq;
    end
end
r_exp = a_exp-b_exp;
for k=1:n
    for j=1:n
        r(k,j).dd = r_exp(2*k-1,2*j-1);
        r(k,j).dq = r_exp(2*k-1,2*j);
        r(k,j).qd = r_exp(2*k,2*j-1);
        r(k,j).qq = r_exp(2*k,2*j);
    end
end
yout=r;
end