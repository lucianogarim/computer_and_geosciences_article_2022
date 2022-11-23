function source = propagation(dir,depth)

source = -2*ones(size(depth)); 
if dir == 0   
    source(end,:) = 0;
elseif dir == 90 
    source(:,1) = 0;
elseif dir == 180 
    source(1,:) = 0; 
elseif dir == 270 
    source(:,end) = 0;
elseif dir > 0 && dir < 90 
    source(end,1) = 0;
elseif dir > 90 && dir < 180 
    source(1,1) = 0;
elseif dir > 180 && dir < 270 
    source(1,end) = 0;
elseif dir > 270 
    source(end,end) = 0;  
end
source(depth<0) = -2;
