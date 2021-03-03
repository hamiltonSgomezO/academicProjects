%% Algoritmo Aleatorizado GRASP
% Basado en el métodos GRASP, el cual utiliza una
% Lista de Candidatos Restringida de posibles objetos
% a incluir en la solución.
% PARÁMETRO: alpha
% Hamilton Smith Gómez Osorio - 201810016101
% Versión 1 (25-09-2020)
function factibles= Aleatorio_GRASP()
clear
clc
factibles=[];
for iter=1:20
    %% Lectura de datos
    format= 'I%d';  % Formato nombre hoja
    sheet = sprintf(format,iter);
    datos = xlsread('Datos',sheet); % Leer hoja
    
    n= datos(1,1); % Número de variables
    m= datos(1,2); % Número de restricciones
    p= datos(1,3); % Número de Funciones objetivo(F.O)
    
    a= datos(2:m+1,1:n); % Coeficientes restricciones
    b= datos(2:m+1,n+1); % Coeficientes parte derecha
    z= datos(m+2:end,1:n); %Coeficientes F.O
    
    %% Algoritmo GRASP
    tic;
    alpha=0.2; %Parámetro lista RCL
    num= 3;
    sum1= sum(a);
    RCL=[];
    criterio= alpha*(max(sum1)-min(sum1)) + min(sum1);
    for i=1:n
        if(sum1(i) <= criterio)
            s= length(RCL);
            RCL(s+1)= i;
        end
    end
    x=[];
    for j=1:num*n
        x(:,j)= Construccion(RCL,n,m,a,b);
    end  
    [sol, VB]= factibilidad(a,x,b,z);
    
    for te=1:size(sol,1)
        FO(te,:)= (z*VB(:,te))';
    end
    if(a*VB(:,1)<=b)
        factibles(iter)= iter;
        if(size(sol,1)~=1)
            ND= dominancia(sol);
        else
            ND=sol;
        end
    else
        q= a*VB>b;
        q= double(q);
        sum1= sum(q);
        min1= min(sum1);
        sol2=[];
        for re=1:size(sol,1)
            if(sum1(re)==min1)
                ind= size(sol2,1);
                sol2(ind+1,:)= sol(re,:);
            end
        end
        if(size(sol2,1)~=1)
            ND=dominancia(sol2);
        else
            ND= sol2;
        end
    end
    X=[];
    for l=1:size(ND)
        aux= size(X,2);
        [row,col]= find(FO==ND(l,:),1);
        X(:,aux+1)= x(:,row);
    end
    time(iter)= toc;
    t= size(ND,1);
    xlswrite('ResultadoAlgoritmo2.xlsx',t,sheet,'A1');
    for y=1:size(ND,1)
        format= 'A%d';  % Formato nombre casilla
        pos = sprintf(format,y+1);
        Resultado= [sum(X(:,y)) find(X(:,y)==1)' (a*X(:,y))' ND(y,:)];
        xlswrite('ResultadoAlgoritmo2.xlsx',Resultado,sheet,pos);
    end
end
end

function sol= Construccion(RCL,n,m,a,b)
x= zeros(n,1); %Variable binarias
B_aux=zeros(m,1); %Pesos parciales
for i=1:m
    if (b(i)<0)
        for j=1:length(RCL)
            col= RCL(j);
            if(a(i,col)<0)
                B_aux(i)= B_aux(i)+ a(i,col);
            end
        end
        if(B_aux(i)>b(i))
            B_aux(i)=-1000000;
        end
    end
end
cont=0;
while(B_aux<=b)
    if(cont~= length(RCL))
        e= randi([1,length(RCL)]);
        if(x(e,1)~=1)
            for o=1:m
                if(b(o)<0 && a(o,e)<0)
                    b_2(o,1)= B_aux(o);
                else
                    b_2(o,1)= B_aux(o)+a(o,e);
                    if(b_2(o) > b(o))
                        break
                    end
                end
            end
            
            if(b_2<=b)
                B_aux= b_2;
                x(e,1)=1;
            else
                break
            end
            cont=cont+1;
        end
    else
        break;
    end
end
sol=x;
end

function [sol, VB]= factibilidad(a,x,b,z)
Fac=[];
Inf=[];
for i= 1: size(x,2)
    if(a*x(:,i)<=b)
        s= size(Fac,2);
        Fac(s+1,:)= (z*x(:,i))';
        VB1(:,s+1)= x(:,i);
    else
        s= size(Inf,2);
        Inf(s+1,:)= (z*x(:,i))';
        VB2(:,s+1)= x(:,i);
    end
    
    if(~isempty(Fac))
        sol= Fac;
        VB= VB1;
    else
        sol= Inf;
        VB= VB2;
    end
end
end

function ND = dominancia(FO)
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
    end
    if(isempty(find(Do==i,1)))
        ind2= length(nd);
        nd(:,ind2+1)= i;
    end
end
size(nd,2);
for p=1:size(nd,2)
    ND(p,:) = FO(nd(p),:);
end
end