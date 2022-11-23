function [fitness_population] = fitness(population,depth,shadow,dx,window_score,facies)

aux1 = size(facies.Inputs(1, 1).MembershipFunctions,2);
aux2 = size(facies.Inputs(1, 2).MembershipFunctions,2);
fitness_population = zeros(size(population,1),1);

cont = 0;
for k = 1:size(population,1)
    
    waveE = wave(depth, population(k,1), shadow, population(k,2), dx);
    waveE = waveE - min(waveE(:));
    waveE = waveE ./ max(waveE(:));
    
    for i = 1:aux1-1
        facies.Inputs(1, 1).MembershipFunctions(1, i).Parameters(1,end-1:end) = population(k,i+2);
        facies.Inputs(1, 1).MembershipFunctions(1, i+1).Parameters(1,1:2) = population(k,i+2);
    end
    
    for i = 1:aux2-1
        facies.Inputs(1, 2).MembershipFunctions(1, i).Parameters(1,end-1:end) = population(k,i+aux1+1);
        facies.Inputs(1, 2).MembershipFunctions(1, i+1).Parameters(1,1:2) = population(k,i+aux1+1);
    end
    
    facies_aux = evalfis(facies,[depth(:),waveE(:)]);
    facies_aux = reshape(facies_aux,size(depth));

warning('off','all')
    for i = 1:size(facies_aux,1)
        for j = 1:size(facies_aux,2)
            score = window_score(i,j,floor(facies_aux(i,j)));
            cont = score+cont;
        end
    end
    
    fitness_population(k) =  cont/numel(facies_aux);
    cont = 0;
end