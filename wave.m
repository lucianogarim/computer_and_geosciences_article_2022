function  Eb = wave(depth,dir,shadow,h0,dx)

inland = ones(size(depth));
[seawater] = double((depth>0));
inland = inland.*~seawater;
shoaling = 0.99; grav = 9.8; rhow = 1000;
source = propagation(dir,depth); 

[~,~,travel,waveH] = airymodel(dx,shoaling,h0,depth,source,inland,shadow,dir);

for i = 1:size(waveH,1)
    for j = 1:size(waveH,2)
        if (travel(i,j)<0 && depth(i,j)>0) 
            waveH(i,j) = waveH(i,j)*0.05;
        elseif waveH(i,j) > h0*1.25
                waveH(i,j) = h0*1.25;
        end
    end
end
     
Hb = 0.78.*depth;
Hb(depth<0) = 0;

waveH = imgaussfilt(waveH,1);     
for i = 1:size(waveH,1)
    for j = 1:size(waveH,2) 
        if (depth(i,j) <= 0)
           waveH(i,j) = 0;
        end
        if(waveH(i,j) > Hb(i,j))
            waveH(i,j) = Hb(i,j);
        end
    end
end
waveH(depth<0) = 0;

E = 0.5*grav*rhow.*waveH.^2;
Eb = E.*tanh(exp(-0.24.*depth));
Eb = imgaussfilt(Eb,1);


