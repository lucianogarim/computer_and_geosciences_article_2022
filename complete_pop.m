function generation = complete_pop(generation,nPopulationSize,num_member_func_depth)

int_dir = [0,45,90,135,180,225,270,315];
dx = 0.01:0.01:0.1;

while( size(generation,1) < nPopulationSize )
   
    sample = generation(randi(size(generation,1), 1),:);
    alt = sample(2);
    x = sample(3:num_member_func_depth+1);
    y = sample(num_member_func_depth+2:end);
    
    if rand>0.5
        alt = alt+dx(randi(length(dx), 1))*alt;
        new_x(:) = x(:)+randi([-1,1])*dx(randi(length(dx), 1))*x(:);
        new_y(:) = y(:)+randi([-1,1])*dx(randi(length(dx), 1))*y(:);
    else
        alt = alt-dx(randi(length(dx), 1))*alt;
        new_x(:) = x(:)-randi([-1,1])*dx(randi(length(dx), 1))*x(:);
        new_y(:) = y(:)-randi([-1,1])*dx(randi(length(dx), 1))*y(:);
    end
    
    generation = vertcat(generation,[int_dir(randi(length(int_dir), 1)) alt new_x new_y]);
end