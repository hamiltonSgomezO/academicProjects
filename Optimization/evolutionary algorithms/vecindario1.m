function [N1]=vecindario1(x,nSol,n)
N1(n,nSol)=0;
%% Vecindario 1
for j=1:nSol
    new=x;
    R1= randi(n,2,1);
    if(x(R1(1))==0)
        while(R1(2)==0)
            R1(2)= randi(n);
        end
        new(R1(1))=1;
        new(R1(2))=0;
    else
        while(R1(2)==1)
            R1(2)= randi(n);
        end
        new(R1(1))=0;
        new(R1(2))=1;
    end
    N1(:,j)= new;
end
end