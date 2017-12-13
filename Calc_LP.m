function [P,dPdth] = Calc_LP( N,costh )
%{
This function calculates P_j and its theta-derivative from 0 to N at a 
given theta by the recursive formula. In this function P(1) is P_0, P(2) 
is P_1, etc.

The recursive forluma is as following:
P(i)=(2*i-3)/(i-1)*costh*P(i-1)-(i-2)/(i-1)*P(i-2), P(1)=1, P(2)=costh;
dPdth(i)=(2*i-3)*P(i-1)+dPdth(i-2), dPdth(1)=0, dPdth(2)=-sqrt(1-costh^2);

If costh is close enough to 1 or -1 (within 1E-18), the function treats 
costh as 1 or -1, and use property of Legendre polynomial to give the value
of P_j and dP_j/dtheta at costh.
%}
P=zeros(N+1,1);
dPdth=zeros(N+1,1);
P(1)=1;
dPdth(1)=0;
if N>=1
tol=1E-18;
if abs(costh-1)<tol
    for i=2:N+1
        P(i)=1;
        dPdth(i)=0;
    end
else
    if abs(costh+1)<tol
        for i=2:N+1
            P(i)=power(-1,i+1);
            dPdth(i)=0;
        end
    else

        P(2)=costh;
        dPdth(2)=-sqrt(1-costh^2);
        for i=3:N+1
            P(i)=(2*i-3)/(i-1)*costh*P(i-1)-(i-2)/(i-1)*P(i-2);
            dPdth(i)=-sqrt(1-costh^2)*(2*i-3)*P(i-1)+dPdth(i-2);
        end
    end
end
end
end