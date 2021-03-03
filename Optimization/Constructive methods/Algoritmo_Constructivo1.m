%% Algoritmo constructivo 1
% Enfocado en el cumplimiento de las restricciones
% integrando los items con suma de pesos menor para
% cada una satisfacer cada una de ellas
% Hamilton Smith Gómez Osorio - 201810016101
%Versión 1 (04-09-2020) 

function FO=Algoritmo_Constructivo1()
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
    x= solucion_iterativa(n,m,a,b);
    g=0;
if(a*x<=b)
    g=0;
else
    g= sum(((-1)*abs(a*x-b))/3);
end
FO_p= (z*x)+g; %Función objetivo penalizada
(z*x)'
FO(:,iter)=FO_p;
time(iter)= toc;
Resultado= [sum(x)' find(x==1)' (a*x)' FO_p'];
xlswrite('ResultadoAlgoritmo1.xlsx',1,sheet,'A1');
xlswrite('ResultadoAlgoritmo1.xlsx',Resultado,sheet,'A2');

end
time'
x'
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
