%% Algoritmo Aleatorizado con Ruido
% Enfocado en el cumplimiento de las restricciones
% integrando los items con suma de pesos menores, en
% donde a estos se le han agregado ruido de una 
% distribución normal.
% PARÁMETRO: desv -> desviación estandar
% Hamilton Smith Gómez Osorio - 201810016101
% Versión 1 (25-09-2020)

function factibles= Aleatorio_Ruido()
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
    tic;
    desv= 200;  %Parámetro desviación estandar
    num=1;
    x=[];
    for i=1:num*n
        ruido= normrnd(0,desv,[1,n]);
        x(:,i)= solucion_iterativa(n,m,a,b,ruido);
    end
    [sol, VB]= factibilidad(a,x,b,z);
    for te=1:size(sol,1)
        FO(te,:)= (z*VB(:,te))';
    end
    if(a*VB(:,1)<=b)
         factibles(iter)=iter;
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
    xlswrite('ResultadoAlgoritmo1.xlsx',t,sheet,'A1');
    for y=1:size(ND,1)
        format= 'A%d';  % Formato nombre hoja
        pos = sprintf(format,y+1);
        Resultado= [sum(X(:,y)) find(X(:,y)==1)' (a*X(:,y))' ND(y,:)];
        xlswrite('ResultadoAlgoritmo1.xlsx',Resultado,sheet,pos);
    end
end

end

function x= solucion_iterativa(n,m,a,b,ruido)
%% Método centrado en restricciones

x= zeros(n,1); %Variable binarias
B_aux=zeros(m,1); %Pesos parciales
for i=1:m
    if (b(i)<0)
        for j=1:n
            if(a(i,j)<0)
                B_aux(i)= B_aux(i)+ a(i,j);
            end
        end
    end
end
M_select= sum(a) + ruido;
b_2=zeros(m,1);
for k=1:n
    [M,I]= min(M_select); %Citrerio de selección
    for o=1:m
        if(b(o)<0 && a(o,I)<0)
            b_2(o)= B_aux(o);
        else
            b_2(o)= B_aux(o)+a(o,I);
            if(b_2(o) > b(o))
                break
            end
        end
    end
    if(b_2<= b)
        B_aux= b_2;
        x(I,1)=1;
    else
        for q=1:m
            if(b(q)<0 && a(q,I)<0)
                B_aux(q) = B_aux(q) - a(q,I);
            end
        end
    end
    M_select(1,I)= 50000; %Número M >>0
end
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

for p=1:size(nd,2)
    ND(p,:) = FO(nd(p),:);
end
end