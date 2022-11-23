%% Script manuscript
clear all

% load colormaps
facies_CM = load('Facies_CM.mat').Facies_CM;
seismic_facies_CM = load('Seismic_CM.mat').Facies_CM;

% load seismic data
seismic_facies = import_seismic_facies('113_0_Unit_A_Upper_RIFE_SF+AF.grd', 6, 645);
y_seismic = (7134000:50:7166000-50);  x_seismic = (672000:50:691000-50);
[X_seismic,Y_seismic] = meshgrid(x_seismic,y_seismic);

%load paleo data
paleo = importapaleo("113.0_Base_Sal_Paleo_x500.grd", [6, Inf]);
y_paleo = (7134000:500:7166000-500); x_paleo = (672000:500:691000-500);
[X_paleo,Y_paleo] = meshgrid(x_paleo,y_paleo);

depth =  paleo; shadow = 1; dx = 500; opt = 1;

[window_score] = movewindow(seismic_facies,depth);

%% Initial parameters
dir = 90; h0 = 1;
waveE = wave(depth, dir, shadow, h0, dx);
facies_aux = readfis('Facies.fis');
facies = evalfis(facies_aux,[depth(:),waveE(:)]);
facies_initial = reshape(round(facies),size(depth));
num_member_func_depth = size(facies_aux.Inputs(1, 1).MembershipFunctions,2);
num_member_func_energy = size(facies_aux.Inputs(1, 2).MembershipFunctions,2);

for i = 1 : (num_member_func_depth-1)
    x(i) = facies_aux.Inputs(1, 1).MembershipFunctions(1, i).Parameters(1,3);
end

list=[];
for i = 1 : num_member_func_energy
    for j = 1 : 4
        list = horzcat(list,facies_aux.Inputs(1, 2).MembershipFunctions(1,i).Parameters(1,j));
    end
end

list = unique(list);
list(end) = max(max(waveE));
y = list(2:end-1)/max(max(waveE));

facies_aux.Inputs(1, 2).Range(1,end) = 1;
facies_aux.Inputs(1, 2).MembershipFunctions(1, 1).Parameters(1,1:2) = 0;
facies_aux.Inputs(1, 2).MembershipFunctions(1, num_member_func_energy).Parameters(1,end-1:end) = 1;

%% Genetic algorithm hiperparameters
nPopulationSize = 30; n_gen = 10; n_sel = 16;
probabilityOfCrossOver = 0.5; p_cross = 2;
sol = [];

%% Initial population
[population] = first_pop(dir, h0, x, y ,nPopulationSize);
generation = population;

%% Main Loop
for i = 1:n_gen

    [fitness_population] = fitness(generation, depth, shadow, dx, window_score, facies_aux);
    generation = horzcat(generation, fitness_population);
    generation = sortrows(generation, size(generation,2),'descend');
    
    while ( nPopulationSize < size(generation,1))
        generation(end,:)=[];
    end
   
    sol = vertcat(sol, generation);
    best = unique(sol, 'rows');
    best = sortrows(best, size(population,2)+1);
    best = best(end-2:end, 1:size(population,2)+1);
   
    [selections] = selection(generation, n_sel);
    
    [generation] = reproduction(selections, probabilityOfCrossOver, p_cross);
    generation = unique(generation, 'rows');
    
    generation = complete_pop(generation, nPopulationSize,num_member_func_depth);
end

best = best(end,:); dir = best(1); h0 = best(2);
waveE = wave(depth,dir, shadow, h0, dx);

best(num_member_func_depth+2:end-1) = max(max(waveE)).*best(num_member_func_depth+2:end-1);
facies_aux.Inputs(1, 2).Range(1,end) = max(max(waveE));
facies_aux.Inputs(1, 2).MembershipFunctions(1, num_member_func_energy).Parameters(1,end-1:end) = max(max(waveE));

for i = 1:num_member_func_depth-1
    facies_aux.Inputs(1, 1).MembershipFunctions(1, i).Parameters(1,end-1:end) = best(i+2);
    facies_aux.Inputs(1, 1).MembershipFunctions(1, i+1).Parameters(1,1:2) = best(i+2);
end

for i = 1:num_member_func_energy-1
    facies_aux.Inputs(1, 2).MembershipFunctions(1, i).Parameters(1,end-1:end) = best(i+num_member_func_depth+1);
    facies_aux.Inputs(1, 2).MembershipFunctions(1, i+1).Parameters(1,1:2) = best(i+num_member_func_depth+1);
end

facies = evalfis(facies_aux,[depth(:),waveE(:)]);
facies = reshape(floor(facies),size(depth));
facies_genetic = facies;
figure(3)
surf(X_paleo,Y_paleo,-depth,facies,'EdgeColor','none')
axis image
title('Simulated Facies')
c = colorbar;
c.TickLabels = {'Reworked','Stromatolic','Transitional','Non-reservoir'};
c.Ticks =  [1,2,3,4];
colormap(facies_CM)
caxis([0 5])
set(gcf,'color','w');

figure(4)
surf(X_paleo,Y_paleo,-depth,facies_initial,'EdgeColor','none')
axis image
title('Initial Facies')
c = colorbar;
c.TickLabels = {'Reworked','Stromatolic','Transitional','Non-reservoir'};
c.Ticks =  [1,2,3,4];
colormap(facies_CM)
caxis([0 5])
set(gcf,'color','w');

figure(5)
surf(X_seismic,Y_seismic,seismic_facies,'EdgeColor','none')
axis image
title('Seismic Facies')
d = colorbar;
d.Ticks =  [1,2,3,4];
d.TickLabels = {'S1','S2','S3','S4'};
colormap(seismic_facies_CM)
set(gcf,'color','w');

%% Well adjustment
field1 = 'pocos';
value1 = {"2D";"3";"4D";"5";"6";"7D";"8";"10";"11";"12DA";"13";"15D";"17";"19D";"20D";"23";"55";"69";"77A";"82A";"95";"97"};
field2 = 'X';
value2 = {677859;675904;679061;677848;684656;682515;680980;683623;677112;677850;685128;680886;684517;679908;683110;681190;677915;683110;678651;678471;685582;679100};
field3 = 'Y';
value3 = {7140905;7138953;7149235;7146823;7162911;7156620;7152820;7161727;7136745;7139675;7155027;7152002;7161201;7148113;7157356;7156260;7144541;7157656;7148241;7138793;7157336;7145771};
field4 = 'AF'; % 1-Reword 2-Stromatolic  3-Transitionals 4-non-reservoir
af113 = {1;3;2;1;2;2;3;2;2;2;2;2;3;1;1;3;1;2;2;2;2;2};
af114 = {3;2;2;2;2;2;2;2;2;3;3;2;1;3;2;3;2;2;3;3;3;2};
af115 = {2;3;2;3;2;3;3;2;4;4;3;3;2;3;3;3;3;3;3;3;3;3};

wells = struct(field1,value1,field2,value2,field3,value3,field4,af113);

X = [wells(:).X];
Y = [wells(:).Y];
well_facies = [wells(:).AF];

aux_X = X; aux_Y = Y;
aux = well_facies;
I = [X_paleo(:)  Y_paleo(:)];

indicators = create_indicators(well_facies);

%Loop
for i = 1: length(I)

    I_u = I(i,:); 

    X = [X I_u(1)];
    Y = [Y I_u(2)];

    distance = distance_calculation(X,Y);
    correlations = correlation(1000, 2000, 75, 75, distance, 'gau');
    
    lambdas = linear_solver(correlations); 
    [~,new_facies] = prob_calculation(lambdas, indicators);
    
    aux = [aux new_facies];
    aux_X = [aux_X I_u(1)];
    aux_Y = [aux_Y I_u(2)];
    
    X = [wells(:).X];
    Y = [wells(:).Y];   
end
i = length(I);
while i ~= 1
    if mod(aux_X(:,i),500) ~= 0
        aux_X(:,i) = [];
        aux_Y(:,i) = [];
        facies(:,i) = [];
        aux(:,i) = [];
    end
    i=i-1;
end

for i = 1:size(X_paleo(:),1)
    [krig_one, ~] = Simple_kriging([X_paleo(i) Y_paleo(i)], [X ;Y]', ones(size(1,22)), 0, 1, 'gau', 1000, 2000, 75, 75);
    One(i) = krig_one;
end

Aux = [];
k = 1;
for i = 1: size(X_paleo,2)
	for j = 1:size(Y_paleo,1)
		Aux(j,i) = One(k);
		k = k+1;
	end
end
One = Aux;
aux(end)=[];
facies = reshape(aux,[64,38]);
Final = round(One.*facies + (1-One).*facies_genetic);

figure(6)
surf(X_paleo,Y_paleo,-depth,Final,'EdgeColor','none')
axis image
title('Final Model')
c = colorbar;
c.TickLabels = {'Reworked','Stromatolic','Transitional','Non-reservoir'};
c.Ticks =  [1,2,3,4];
colormap(facies_CM)
caxis([0 5])
set(gcf,'color','w');
hold on

function [Krig, Variance] = Simple_kriging(coord_celll, coord_samples, sample_values, mean_variable, variance, type, ...
    min_corr, max_corr, azim, theta)

    warning('off','all')
    dist = squareform(pdist([coord_celll; coord_samples]));
    dist_vector = dist(2:end,1);
    dist_matrix = dist(2:end,2:end);
    
    vetor_krig = variance*correlation(min_corr, max_corr, azim, theta, dist_vector, type);
    matrix_krig = variance*correlation(min_corr, max_corr, azim, theta, dist_matrix, type);
    
    Weight = matrix_krig\vetor_krig;
    Krig = mean_variable+sum(Weight.*(sample_values-mean_variable));
    Variance = variance-sum(Weight.*vetor_krig);
end

function C = correlation(ComprCorrelMin, ComprCorrelMax, Azim, Theta, Dist, Tipo)
    
    switch Tipo
        case 'exp'
            C = ExpCov(Dist, radial_correlation(ComprCorrelMin, ComprCorrelMax, Azim, Theta));    
        case 'gau'
            C = GauCor(Dist, radial_correlation(ComprCorrelMin, ComprCorrelMax, Azim, Theta));
        case 'sph'
            C = SphCor(Dist, radial_correlation(ComprCorrelMin, ComprCorrelMax, Azim, Theta));
    end  
end

function L = radial_correlation(ComprCorrelMin, ComprCorrelMax, Azim, Theta)
  
  L = sqrt((ComprCorrelMin^2*ComprCorrelMax^2)./(ComprCorrelMax^2*(sin(Azim-Theta)).^2+ComprCorrelMin^2*(cos(Azim-Theta)).^2));

end

function C = GauCor(Dist,L)
    
    C = exp(-3*Dist.^2/L.^2);
    
end

function C = SphCor(Dist,L)

    C = zeros(size(Dist));   
    C(Dist<=L) = 1-3/2*Dist(Dist<=L)/L+1/2*Dist(Dist<=L).^3/L^3;

end

function C = ExpCor(Dist,L)
    
    C = exp(-3*Dist/L);
    
end