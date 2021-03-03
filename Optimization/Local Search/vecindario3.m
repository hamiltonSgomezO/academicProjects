function [N3]= vecindario3(x,nSol)
one= find(x==1);
zero= find(x==0);
len=size(x,1);
N3(1:len,nSol)=0;
R3(3)=0;
for j=1:nSol
    new=x;
    if(size(zero,1)>1)
        R3(1)= randi([1 size(zero,1)]);
        R3(2)= randi([1 size(zero,1)]);
        while(R3(1)==R3(2))
            R3(2)= randi([1 size(zero,1)]);
        end
        R3(3)= randi([1 size(one,1)],1,1);
        new(zero(R3(1)))=1;
        new(zero((R3(2))))=1;
        new(one(R3(3)))=0;
    end
    N3(:,j)= new;
end
end