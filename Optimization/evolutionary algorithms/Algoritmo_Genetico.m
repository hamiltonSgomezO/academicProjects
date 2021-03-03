function Algoritmo_Genetico()
clear
clc
periodos= [3 4 8 12 13 14 19 20];
for it=6:8
    iter= periodos(1,it)
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
    t_max=300;
    time=0;
    num=10; %número de individuos en la población
    C= 100; % Número de hijos
    p=0.2; % Probabilidad de mutación
    X0=[];
    X0(n,num)=0;
    for i=1:num
        S_i= Soluciones_Iniciales(n,m,a,b);
        if(S_i == X0)
            S_i= Soluciones_Iniciales(n,m,a,b);
        else
            X0(:,i)= S_i;
        end
    end
    sizX=size(X0,2);
    if(sizX<num)
        sol= rand(num);
        N1=vecindario1(sol,num-sizX,n);
        X0= [X0, N1];
    end
    %X0
    while(time<t_max)
        for new_C=1:C
            padres= Seleccion(X0,n,num,a,b,z);
            %padres
            hijos= Cruce(padres,n);
            %hijos
            %% Mutación por vecindarios
            for j=1:2
                r_mut= rand();
                if(r_mut<=p)
                    hijos(:,j)= vecindario1(hijos(:,j),1,n);
                end
            end
            G1= sum(double(z*hijos(:,2)<z*hijos(:,1)));
            G2= sum(double(z*hijos(:,1)<z*hijos(:,2)));
            if(G1>=G2)
                hijo= hijos(:,1);
            else
                hijo= hijos(:,2);
            end
            %hijo
            time=toc;
            if(time<t_max)
                FP= AlgoritmoBusquedaLocal_2(hijo,n,a,b,z,time);
                if(new_C==1)
                    X1= FP;
                else
                    X1=[X1,FP];
                end
                %X1
            else
                break
            end
        end
        
        %Actualización de la población
        P1=[X0,X1];
        %SizeP1= size(P1,2)
        P=[];
        for col=1:size(P1,2)-1
            for col2=col+1:size(P1,2)
                if(P1(:,col)==P1(:,col2))
                    break
                elseif(col2==size(P1,2))
                    sizP= size(P,2);
                    P(:,sizP+1)=P1(:,col);
                end
            end
        end
        %SizeP= size(P,2)
        [NP,X]= Actualizacion(a,P,b,z,num);
        %NP;
        X0=X;
        time=toc;
    end
    [Fac,VB_F,~,~]= factibilidad(a,P,b,z);
    [N,nd,~,~]= dominancia(Fac);
    ND= N;
    X= VB_F(:,nd);
    time(it,1)=toc;
    a*X(:,1)<=b
    t= size(X,2);
    %a*X<=b;
    xlswrite('ResultadoAlgoritmoGenetico7.xlsx',t,sheet,'A1');
    for or=1:size(X,2)
        format= 'A%d';  % Formato nombre de la fila
        pos = sprintf(format,or+1);
        Resultado= [sum(X(:,or)) find(X(:,or)==1)' (a*X(:,or))' ND(or,:)];
        xlswrite('ResultadoAlgoritmoGenetico7.xlsx',Resultado,sheet,pos);
    end
    
end
time
end

function padres= Seleccion(N,n,num,a,b,z)
padres(n,2)=0;
for i=1:2
    r= randi(num,2); %Individuos aleatorios
    while(r(1)==r(2))
        r(2)= randi(num);
    end
    FO1= a*N(:,r(1))<=b;
    FO2= a*N(:,r(2))<=b;
    if(~FO1)
        if(FO2)
            padres(:,i)= N(:,r(2));
        else
            if(sum(double(FO1)>sum(double(FO2))))
                padres(:,i)= N(:,r(1));
            elseif(sum(double(FO2)>sum(double(FO1))))
                padres(:,i)= N(:,r(2));
            else
                padres(:,i)= N(:,r(2));
            end
            %padres
        end
    else
        if(FO2)
            sum1=sum(double(z*N(:,r(2))<z*N(:,r(1))));
            sum2=sum(double(z*N(:,r(1))<z*N(:,r(2))));
            if(sum1==sum2)
                padres(:,i)= N(:,r(2));
            else
                if(sum1>sum2)
                    padres(:,i)= N(:,r(1));
                else
                    padres(:,i)= N(:,r(2));
                end
                %padres
            end
        else
            padres(:,i)= N(:,r(1));
        end
    end
end

end

function hijos= Cruce(padres,n)
Pc= randi(n);
hijos(:,1)= [padres(1:Pc,1);padres(Pc+1:end,2)];
hijos(:,2)= [padres(1:Pc,2);padres(Pc+1:end,1)];
end


function [ND,X]= Actualizacion(a,P,b,z,num)
[Fac,VB_F,Infac,VB_I]= factibilidad(a,P,b,z);
B_Fac= size(Fac,1);
ND=[];
X=[];
if(B_Fac>=num)
    while(size(ND,1) <num)
        [N,nd, DO, do]= dominancia(Fac);
        siz= size(N,1);
        if(siz>=num)
            ND= N(1:num,:);
            X= VB_F(:,nd(1:num));
        else
            if(size(ND,1)==0)
                ND= N;
                X= VB_F(:,nd);
            else
                ND= [ND; N];
                X= [X, VB_F(:,nd)];
            end
            Fac= DO;
            VB_F= VB_F(:,do);
        end
    end
else
    aun=num-B_Fac;
    ND= Fac;
    X= VB_F;
    size(P,2);
    size(ND,1);
    while(size(ND,1)<num)
        [solB_Inf, indB,solW_Inf,indW]= mejor_Infac(a,VB_I,b,Infac);
        if(size(solB_Inf,1)>=aun)
            if(size(ND,1)==0)
                ND=solB_Inf(1:aun,:);
                X= indB(:,1:aun);
            else
                ND= [ND; solB_Inf(1:aun,:)];
                X= [X, indB(:,1:aun)];
            end
        else
            if(size(ND,1)==0)
                ND= solB_Inf;
                X= indB;
            else
                ND= [ND; solB_Inf];
                X= [X, indB];
            end
            Infac= solW_Inf;
            VB_I= indW;
            aun= num-size(ND,1);
        end
    end
end
end

function [solB_Inf, indB,solW_Inf,indW]= mejor_Infac(a,VB_I,b,Infac)
q= double(a*VB_I>b);
sum1= sum(q);
min1= min(sum1);
solB_Inf=[];
solW_Inf=[];
indB=[];
indW=[];
for re=1:size(Infac,1)
    if(sum1(re)==min1)
        ind= size(solB_Inf,1);
        solB_Inf(ind+1,:)= Infac(re,:);
        indB(:,ind+1)= VB_I(:,re);
    else
        ind= size(solW_Inf,1);
        solW_Inf(ind+1,:)= Infac(re,:);
        indW(:,ind+1)= VB_I(:,re);
    end
end
end

function [Fac,VB_F,Infac,VB_I]= factibilidad(a,x,b,z)
Fac=[];
Infac=[];
VB_I=[];
VB_F=[];

for i= 1: size(x,2)
    if(a*x(:,i)<=b)
        s= size(Fac,1);
        Fac(s+1,:)= (z*x(:,i))';
        VB_F(:,s+1)= x(:,i);
    else
        s= size(Infac,1);
        Infac(s+1,:)= ((z*x(:,i))');
        VB_I(:,s+1)= x(:,i);
    end
end
end

function [ND,nd, DO, ddo] = dominancia(FO)
%FO= unique(FO,'rows','stable');
nd=[];
do=[];
for i=1:size(FO,1)
    for j=1:size(FO,1)
        if(FO(i,:)< FO(j,:))
            if(isempty(find(do==i,1)))
                ind= size(do,1);
                do(1,ind+1)= i;
            end
        end
        if(FO(i,1)==FO(j,1))
            if(FO(i,2)==FO(j,2))
                if(FO(i,3)<FO(j,3))
                    if(isempty(find(do==i,1)))
                        ind= size(do,1);
                        do(1,ind+1)= i;
                    end
                end
            elseif(FO(i,3)==FO(j,3))
                if(FO(i,2)<FO(j,2))
                    if(isempty(find(do==i,1)))
                        ind= size(do,1);
                        do(1,ind+1)= i;
                    end
                end
            end
        elseif(FO(i,2)==FO(j,2))
            if(FO(i,3)==FO(j,3))
                if(FO(i,1)<FO(j,1))
                    if(isempty(find(do==i,1)))
                        ind= size(do,1);
                        do(1,ind+1)= i;
                    end
                end
            end
        end
        if(FO(i,1)==FO(j,1))
            if(FO(i,2:3)<FO(j,2:3))
                if(isempty(find(do==i,1)))
                    ind= size(do,1);
                    do(1,ind+1)= i;
                end
            end
        elseif(FO(i,2)==FO(j,2))
            if(FO(i,1)<FO(j,1) && FO(i,3)<FO(j,3))
                if(isempty(find(do==i,1)))
                    ind= size(do,1);
                    do(1,ind+1)= i;
                end
            end
        elseif(FO(i,3)==FO(j,3))
            if(FO(i,1:2)<FO(j,1:2))
                if(isempty(find(do==i,1)))
                    ind= length(do);
                    do(1,ind+1)= i;
                end
            end
        end
    end
    if(isempty(find(do==i,1)))
        ind2= size(nd,2);
        nd(1,ind2+1)= i;
    end
end
ND=[];
for p=1:size(nd,2)
    ND(p,:) = FO(nd(p),:);
end
ddo=[];
DO=[];
if(isempty(nd))
    ddo= 1:size(FO,1);
    DO= FO;
else
    for j=1:size(FO,1)
        tl=size(ddo,2);
        if(isempty(find(nd==j, 1)))
            ddo(1,tl+1)= j;
            DO(tl+1,:)=FO(j,:);
        end
    end
end
end
