%% Algoritmo constructivo 2
% Enfocado en la búsqueda de una función objetivo mayor
% integrando a la solucion los items que producen una
% mayor proporción de ganancias en las funciones
% objetivo, sumando estos valores para cada una.
% Hamilton Smith Gómez Osorio -201810016101
%Versión 1 (04-09-2020)
function FO=Algoritmo_Constructivo2()
clear
clc
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
    x= primera_solucion(n,m,a,b,z);
    x=balance(a,b,x);
    g=0;
    if(a*x<=b)
        g=0;
    else
        g= sum(((-1)*abs(a*x-b))/3);
    end
    if(a*x<=b)
        ter(iter)=1;
    else
        ter(iter)=0;
    end
    FO_p= (z*x)+g; %Función objetivo penalizada
    FO(:,iter)= FO_p;
    time(iter)= toc;
    Resultado= [sum(x)' find(x==1)' (a*x)' FO_p'];
    xlswrite('ResultadoAlgoritmo2.xlsx',1,sheet,'A1');
    xlswrite('ResultadoAlgoritmo2.xlsx',Resultado,sheet,'A2');
end
time'
ter'
end

function x= primera_solucion(n,m,a,b,z)
%% Método centrado en restricciones

x= zeros(n,1); %Variable binarias
B_aux=zeros(m,1); %Pesos parciales
FO= 0;      %Función objetivo parcial
for i=1:m
    if (b(i)<0)
        for j=1:n
            if(a(i,j)<0)
                B_aux(i)= B_aux(i)+ a(i,j);
            end
        end
    end
end
again=[];
M_select= sum(z);
b_2=zeros(m,1);
for k=1:n
    [M,I]= max(M_select); %Citrerio de selección
    for o=1:m
        if(b(o)<0 && a(o,I)<0)
            b_2(o)= B_aux(o);
        else
            b_2(o)= B_aux(o)+a(o,I);
        end
    end
    b_2;
    if(b_2<= b)
        B_aux= b_2;
        x(I,1)=1;
        if(find(again,1)==I)
            again(find(again,1))=[];
        end
        if(a*x <= b)
            break;
        end
    else
        b_3=[];
        for q=1:m
            if(b(q)<0 && a(q,I)<0)
                b_3(q,1)= q;
                b_3(q,2)= B_aux(q) - a(q,I);
                if(b_3(q,2) > b(q))
                    break;
                end
            end
        end
        if(size(b_3(:,2))== size(b))
            if(b_3(:,2) <= b)
                for kt=1:size(b_3)
                    if(b_3(kt,1)~=0)
                        in= b_3(kt,1);
                        B_aux(in) = b_3(in,2);
                    end
                end
            end
        else
            break;
        end
    end
    B_aux;
    if(isempty(find(again==I,1)))
        M_select(1,I)= -50000; %Número M <<0
    end
    left= sum(M_select~=-50000);
end
end


function x= balance(a,b,x)
t= size(find(x==1));
q= size(find(x==0));
sum1(1,t)= 0;
sum0(1,q)= 0;

for op=1:t
    if(a*x<=b)
        break
    end
    one= find(x==1);
    zero= find(x==0);
    for i=1:size(one)
        for j=1:size(a,1)
            if(a(j,one(i))<0)
                sum1(1,i)= sum1(1,i)+a(j,one(i));
            end
        end
    end
    for k=1:size(zero)
        for pq=1:size(a,1)
            if(a(pq,zero(k))<0)
                sum0(1,k)= sum0(1,k)+a(pq,zero(k));
            end
        end
    end
    [M1,I1]= min(sum0);
    [M2,I2]= max(sum1);
    if(M1<M2)
        x(zero(I1))=1;
        x(one(I2))=0;
    else
        sum1(I2)=0;
    end
end
end
