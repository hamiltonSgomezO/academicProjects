function [S1]= perturbacion(x,nMov)
len=size(x,1);
%% Vecindario 1
for i=1:nMov
    R1= randi([1 len],2,1);
    if(x(R1(1))==0)
        while(R1(2)==0)
            R1(2)= randi([1 len]);
        end
        x(R1(1))=1;
        x(R1(2))=0;
    else
        while(R1(2)==1)
            R1(2)= randi([1 len]);
        end
        x(R1(1))=0;
        x(R1(2))=1;
    end
end
S1=x;
end