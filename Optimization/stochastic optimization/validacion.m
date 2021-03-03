function [linsol opt]= validacion(M_g)
for j=1:4
    M= M_g;
    for m=1:M
        n_1(m)= unifrnd(-0.8,0.8);
        n_2(m)= unifrnd(0,1.84);
        e_1(m)= unifrnd(-30.91,30.91);
        e_2(m)= unifrnd(-23.18,23.18);
    end
    n1= sum((1/M)*n_1);
    n2= sum((1/M)*n_2);
    e1= sum((1/M)*e_1);
    e2= sum((1/M)*e_2);
    
    
    x1= optimvar('x1');
    x2= optimvar('x2');
    
    linprob = optimproblem('Objective', x1*((2+n1)*2 + 12) + x2*(6 + (3.4+n2)*2));
    
    linprob.Constraints.cons1= x1 >= 180 + e1;
    linprob.Constraints.cons2= x2 >= 162 + e2;
    
    linsol(j)= solve(linprob);
    opt(j)= evaluate(linprob.Objective,linsol(j));   
end
end