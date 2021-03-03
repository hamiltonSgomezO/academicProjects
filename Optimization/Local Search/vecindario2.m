function [N2]= vecindario2(x,nSol)
one=find(x==1);
zero=find(x==0);
len=size(x,1);
N2(1:len,nSol)=0;
R2(3)=0;
for j=1:nSol
    new=x;
    if(size(one,1)>1)
        R2(1)= randi([1 size(one,1)]);
        R2(2)= randi([1 size(one,1)]);
        while(R2(1)==R2(2))
            R2(2)= randi([1 size(one,1)]);
        end
        R2(3)= randi([1 size(zero,1)]);
        new(one(R2(1)))=0;
        new(one(R2(2)))=0;
        new(zero(R2(3)))=1;
    end
    N2(:,j)= new;
end 
end