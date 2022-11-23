function [window_score] = movewindow(seismic,depth)

weight_facies = [0 0 0.5 1; 0 0.6 1 0.3; 0.3 1 0.6 0; 1 0.6 0.3 0]; 

num_cell_seismic = numel(seismic);
num_cell_facies = numel(depth);

size_window = sqrt(num_cell_seismic/num_cell_facies);

step1 = (1:size_window:size(seismic,1));
step2 = (size_window:size_window:size(seismic,1));

weight_S1 = 0; weight_S2 = 0; weight_S3 = 0; weight_S4 = 0;

for i = 1:size(depth,1)
    for j = 1:size(depth,2)
        for i2 = step1(i):step2(i)
            for j2 = step1(j):step2(j)
                switch seismic(i2,j2)
                    case 1 % poderia ser 1 2 3 4
                        weight_S1 = weight_S1 +1;
                    case 2
                        weight_S2 = weight_S2 +1;
                    case 3
                        weight_S3 = weight_S3 +1;
                    case 4
                        weight_S4 = weight_S4 +1;
                end
            end
        end
        weight_S1 = weight_S1/size_window^2;
        weight_S2 = weight_S2/size_window^2;
        weight_S3 = weight_S3/size_window^2;
        weight_S4 = weight_S4/size_window^2;
        weight_S = [weight_S1; weight_S2 ; weight_S3; weight_S4];
 
        window_score (i,j,1:4)= max(weight_S.*weight_facies);
        
        weight_S1 = 0; weight_S2 = 0; weight_S3 = 0; weight_S4 = 0;
    end
end
end