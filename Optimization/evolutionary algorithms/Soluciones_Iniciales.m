function N= Soluciones_Iniciales(n,m,a,b)
%% Parámetros
%n: número de variables del problema
%m: número de restricciones
%a: matriz de costos
%b: vector de valores máximos

desv= 100; %Parámetro Aleatorio ruido
ruido= normrnd(0,desv,[1,n]);
N= solucion_iterativa(n,m,a,b,ruido);

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
