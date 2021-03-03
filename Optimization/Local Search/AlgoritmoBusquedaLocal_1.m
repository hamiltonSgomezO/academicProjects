function AlgoritmoBusquedaLocal_1()
clear
clc
just= [1 2 3 4 7 8 9 10 15 16 20];
for qqq=1:1
    time=0;
    for ttt=1:50
        %for iter=1:20
            iter= just(qqq);
            %% Lectura de datos
            format= 'I%d';  % Formato nombre hoja
            sheet = sprintf(format,iter);
            datos = xlsread('Datos',sheet); % Leer hoja
            
            n= datos(1,1); % Número de variables
            m= datos(1,2); % Número de restricciones
            
            a= datos(2:m+1,1:n); % Coeficientes restricciones
            b= datos(2:m+1,n+1); % Coeficientes parte derecha
            z= datos(m+2:end,1:n); %Coeficientes F.O
            
            nSol= 3*n; %Número de soluciones por vecindario
            
            x= solucion_iterativa(n,m,a,b);
            xInicial=x;
            FP= x;
            check= 0;
            yes=[];
            n=1;
            change=0;
            tic;
            while(time<=300)
                if(check~=1)
                    if(n==1)
                        N=vecindario1(x,nSol);
                    elseif(N==2)
                        N=vecindario2(x,nSol);
                    else
                        N=vecindario3(x,nSol);
                    end
                    FP= [FP N];
                    [ND,X]= Busqueda_local(a,FP,b,z);
                    if(a*X(:,1)<=b)
                        if(size(ND,1)~=1)
                            FP= X;
                            n=1;
                            change=1;
                        elseif(x~=X)
                            FP= X;
                            n=1;
                            change=1;
                        else
                            if(n==1)
                                n=2;
                                change=0;
                            elseif(n==2)
                                n=3;
                                change=0;
                            else
                                FP=X;
                                n=1;
                                change=1;
                            end
                        end
                    else
                        SC1=sum(double(a*x<=b));
                        SC2=sum(double(a*X(:,1)<=b));
                        if(SC2>SC1)
                            FP=X;
                        else
                            if(n==1)
                                n=2;
                                change=0;
                            elseif(n==2)
                                n=3;
                                change=0;
                            else
                                FP=X;
                                n=1;
                                change=1;
                            end
                        end
                    end
                    if(change==1)
                        yes= [yes x];
                        x=[];
                        for s=1:size(FP,2)
                            for y=1:size(yes,2)
                                if(yes(:,y)==FP(:,s))
                                    break
                                elseif(y==size(yes,2))
                                    x= FP(:,s);
                                end
                            end
                        end
                        if(isempty(x))
                            check=1;
                        end
                        
                    end
                    time=toc;
                else
                    break
                end
            end
            time= time +toc;
            %     t= size(X,2);
            %     xlswrite('ResultadoAlgoritmo1.xlsx',t,sheet,'A1');
            %     for or=1:size(FP,2)
            %         format= 'A%d';  % Formato nombre hoja
            %         pos = sprintf(format,or+1);
            %         Resultado= [sum(X(:,or)) find(X(:,or)==1)' (a*X(:,or))' ND(or,:)];
            %         xlswrite('ResultadoAlgoritmo1.xlsx',Resultado,sheet,pos);
            %     end
        %end
    end
    timesito(qqq,1)= time/50;
end
timesito
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