function L=Nueva_frontera()
clear
clc
for iter=1:20
    %% Lectura de datos
    format= 'I%d';  % Formato nombre hoja
    sheet = sprintf(format,iter);
    datos1=[];
    datos2=[];
    datos1 = xlsread('ResultadoAlgoritmo1',sheet); % Leer hoja Algoritmo 1
    datos2 = xlsread('ResultadoAlgoritmo2',sheet); % Leer hoja Algoritmo 2
    A=[];
    B=[];
    F=[];
    for i=2:size(datos1,1)
        l=length(datos1(i,:));
        while(isnan(datos1(i,l)))
            l=l-1;
        end
        A(i-1,:)= datos1(i,l-2:l);
    end
    for j=2:size(datos2,1)
        l=length(datos2(j,:));
        while(isnan(datos2(j,l)))
            l=l-1;
        end
        B(j-1,:)= datos2(j,l-2:l);
    end
    
    Un= [A;B];
    [F,nd] = dominancia(Un);
    int_A=0;

    for k=1:size(A,1)
        for m=1:size(F,1)
            A(k,:);
            F(m,:);
            if(A(k,:)==F(m,:))
                int_A= int_A + 1;
            end
        end
    end
    int_B=0;
    for n=1:size(B,1)
        for o=1:size(F,1)
            if(B(n,:)==F(o,:))
                int_B= int_B + 1;
            end
        end
    end
%     A
%     B
%     F
    L(iter,1)= int_A/size(F,1);
    L(iter,2)= int_B/size(F,1);
    L(iter,3)= int_A/size(A,1);
    L(iter,4)= int_B/size(B,1);
end
end

function [ND,nd] = dominancia(FO)
FO= unique(FO,'rows','stable');
nd=[];
Do=[];
for i=1:size(FO,1)
    for j=1:size(FO,1)
        if(FO(i,:)< FO(j,:))
            if(isempty(find(Do==i,1)))
                ind= length(Do);
                Do(1,ind+1)= i;
            end
        end
        if(FO(i,1)==FO(j,1))
            if(FO(i,2)==FO(j,2))
                if(FO(i,3)<FO(j,3))
                    if(isempty(find(Do==i,1)))
                        ind= length(Do);
                        Do(1,ind+1)= i;
                    end
                end
            elseif(FO(i,3)==FO(j,3))
                if(FO(i,2)<FO(j,2))
                    if(isempty(find(Do==i,1)))
                        ind= length(Do);
                        Do(1,ind+1)= i;
                    end
                end
            end
        elseif(FO(i,2)==FO(j,2))
            if(FO(i,3)==FO(j,3))
                if(FO(i,1)<FO(j,1))
                    if(isempty(find(Do==i,1)))
                        ind= length(Do);
                        Do(1,ind+1)= i;
                    end
                end
            end
        end
        if(FO(i,1)==FO(j,1))
            if(FO(i,2:3)<FO(j,2:3))
                if(isempty(find(Do==i,1)))
                    ind= length(Do);
                    Do(1,ind+1)= i;
                end
            end
        elseif(FO(i,2)==FO(j,2))
            if(FO(i,1)<FO(j,1) && FO(i,3)<FO(j,3))
                if(isempty(find(Do==i,1)))
                    ind= length(Do);
                    Do(1,ind+1)= i;
                end
            end
        elseif(FO(i,3)==FO(j,3))
            if(FO(i,1:2)<FO(j,1:2))
                if(isempty(find(Do==i,1)))
                    ind= length(Do);
                    Do(1,ind+1)= i;
                end
            end
        end
    end
    if(isempty(find(Do==i,1)))
        ind2= length(nd);
        nd(:,ind2+1)= i;
    end
end

for p=1:size(nd,2)
    ND(p,:) = FO(nd(p),:);
end
end