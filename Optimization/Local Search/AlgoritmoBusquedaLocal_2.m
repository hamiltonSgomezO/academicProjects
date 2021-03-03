function AlgoritmoBusquedaLocal_2()
clear
clc
for iter=1:20
    %% Lectura de datos
    format= 'I%d';  % Formato nombre hoja
    sheet = sprintf(format,iter);
    datos = xlsread('Datos',sheet); % Leer hoja
    
    n= datos(1,1); % Número de variables
    m= datos(1,2); % Número de restricciones
    
    a= datos(2:m+1,1:n); % Coeficientes restricciones
    b= datos(2:m+1,n+1); % Coeficientes parte derecha
    z= datos(m+2:end,1:n); %Coeficientes F.O
    
    nSol= 5*n; %Número de soluciones por vecindario
    nMov= 3; %Número de movimientos a hacer al perturbar
    t=300; % Tiempo máximo de corrida
    x= solucion_iterativa(n,m,a,b);
    xInicial=x;
    yes=[];
    tic;
    time=0;
    N=vecindario1(x,nSol);
    FP=[x N];
    [~,X1]= Busqueda_local(a,FP,b,z);
    FP= X1;
    iSol= size(X1,2);
    visit=[];
    visit(:,iSol)=0;
    for i=1:iSol
        check= 0;
        l= randi([1 size(X1,2)]);
        while(~isempty(find(l==visit, 1)))
            l= randi([1 size(X1,2)]);
        end
        visit(i)= l;
        inS= X1(:,l);
        yes= [];
        while(time<=t)
            if(check~=1)
                S1= perturbacion(inS,nMov);
                N=vecindario1(S1,nSol);
                FP=[FP S1 N];
                [ND,X]= Busqueda_local(a,FP,b,z);
                FP= X;
                yes= [yes inS];
                inS=[];
                for s=1:size(FP,2)
                    for y=1:size(yes,2)
                        if(yes(:,y)==FP(:,s))
                            break
                        elseif(y==size(yes,2))
                            inS= FP(:,s);
                        end
                    end
                end
                if(isempty(inS))
                    check=1;
                end
                time=toc;
            else
                break
            end
        end
        if(time>=t)
            break
        end
    end
    tim= size(X,2);
    xlswrite('ResultadoAlgoritmo2.xlsx',tim,sheet,'A1');
    for or=1:size(FP,2)
        format= 'A%d';  % Formato nombre hoja
        pos = sprintf(format,or+1);
        Resultado= [sum(X(:,or)) find(X(:,or)==1)' (a*X(:,or))' ND(or,:)];
        xlswrite('ResultadoAlgoritmo2.xlsx',Resultado,sheet,pos);
    end
end
end

function [ND,X]= Busqueda_local(a,FP,b,z)
[sol, VB]= factibilidad(a,FP,b,z);
if(a*VB(:,1)<=b)
    if(size(sol,1)~=1)
        [ND,nd]= dominancia(sol);
        pos=nd;
    else
        ND=sol;
        pos=1;
    end
else
    [sol2, indi]= mejor_Infac(a,VB,b,sol);
    if(size(sol2,1)~=1)
        [ND,nd]=dominancia(sol2);
        pos= indi(nd);
    else
        ND= sol2;
        pos= indi;
    end
end
X= VB(:,pos);
end

function [sol2, indi]= mejor_Infac(a,VB,b,sol)
q= double(a*VB>b);
sum1= sum(q);
min1= min(sum1);
sol2=[];
for re=1:size(sol,1)
    if(sum1(re)==min1)
        ind= size(sol2,1);
        sol2(ind+1,:)= sol(re,:);
        indi(ind+1,:)= re;
    end
end
end

function x= solucion_iterativa(n,m,a,b)
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
M_select= sum(a);
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
        s= size(Fac,1);
        Fac(s+1,:)= (z*x(:,i))';
        VB1(:,s+1)= x(:,i);
    else
        s= size(Inf,1);
        Inf(s+1,:)= ((z*x(:,i))');
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