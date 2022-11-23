function indicators = create_indicators(facies)
    
    num_facies = max(unique(facies));
    indicators = zeros(length(facies), num_facies);
    
    for i = 1 : length(facies)
        pos = facies(i);
        indicators(i,pos) = 1;
    end

end





    

