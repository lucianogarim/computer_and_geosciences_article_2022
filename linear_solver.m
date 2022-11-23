function lambdas = linear_solver(correlations)
    
    A = correlations(1:end-1,1:end-1);
    B = correlations(1:end-1,end);
    lambdas = linsolve(A,B);
    
end