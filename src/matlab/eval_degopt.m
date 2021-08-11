function F=eval_degopt(degopt,A)
Ha=degopt.Ha;
Hb=degopt.Hb;
y=degopt.y;
m=size(Ha,1)

use_vpa= ~isnumeric(Ha(1,1));

n=size(A,1);
Q={};
F=zeros(n,n);
II=eye(n,n);
ZZ=zeros(n,n);
if use_vpa
    F=vpa(F);
    II=vpa(II);
    ZZ=vpa(ZZ);
end
for k=1:m
    Q{k}=ZZ;
end

for k=1:m
    Z={zeros(n,n), zeros(n,n)};
    if (use_vpa)
        Z{1}=vpa(Z{1});
        Z{2}=vpa(Z{2});
    end

    for j=1:2
        if (j==1)
            coefflist=Ha(k,1:k+1);
            matlist={II, A, Q{1:(k-1)}};
        else
            coefflist=Hb(k,1:k+1);
            matlist={II, A, Q{1:(k-1)}};
        end
        for i=1:length(coefflist)
            Z{j} = Z{j} + coefflist(i)*matlist{i};
        end
    end
    Q{k}=Z{1}*Z{2};
end
matlist={II, A, Q{1:m}};

for k=1:length(y);
    F=F+y(k)*matlist{k};
end
