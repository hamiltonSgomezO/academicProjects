% Algoritmo de optimización estocástica
% Resolución por escenarios

%Tamaño de muestras para entrenamiento
N_g= [1e1 1e2 1e3 1e4];
%Tamaño de muestra para validación
M_g= 1e5;
for j=1:length(N_g)
    N= N_g(j);
    for n=1:N
        %Generación Parámetros estocásticos aleatorios
        n_1(n)= unifrnd(-0.8,0.8);
        n_2(n)= unifrnd(0,1.84);
        e_1(n)= unifrnd(-30.91,30.91);
        e_2(n)= unifrnd(-23.18,23.18);
    end
    % Cálculo parámetro ____
    n1= sum((1/N)*n_1);
    n2= sum((1/N)*n_2);
    e1= sum((1/N)*e_1);
    e2= sum((1/N)*e_2);
    
    %Variables a optimizar
    x1= optimvar('x1');
    x2= optimvar('x2');
    
    % Función objetivo
    linprob = optimproblem('Objective', x1*((2+n1)*2+...
        12) + x2*(6 + (3.4+n2)*2));
    
    %Resticciones
    linprob.Constraints.cons1= x1 >= 180 + e1;
    linprob.Constraints.cons2= x2 >= 162 + e2;
    
    %Solución problema lineal
    linsol(j)= solve(linprob);
    %Evaluacón costo mínimo
    opt(j)= evaluate(...
        linprob.Objective,linsol(j));
    
end
%Proceso de validación
[solVal optVal]= validacion(M_g);

%Reorganización de datos
for i=1: length(optVal)
    sol(:,i)=[opt(i);optVal(i)];
end
%Tabulación de reultados
data= {'Training', 'Validation'};
Solucion= table(sol(:,1),sol(:,2),...
    sol(:,3),sol(:,4),'RowNames',data)
Solucion({'Training','Validation'},:);
solVal

