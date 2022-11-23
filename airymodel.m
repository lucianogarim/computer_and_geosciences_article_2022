function [waveC,waveL,travel,waveH] = airymodel(dx,shoalC,h0,depth,source,inland,shadow,dir)

waveL=zeros(size(depth));
kk=1; 
n_row = size(inland,1);
n_col = size(inland,2);
constant_s = 0.9;

grav  = 9.8;
pi2   = pi*2;
onpi2 = 1/pi2;

iradius=[-2,-1,1,2,-2,-1,0,1,2,-1,0,1,-2,-1,0,1,2,-1,0,-1];
jradius=[0,0,0,0,1,1,1,1,1,2,2,2,-1,-1,-1,-1,-1,-2,-2,-2];

dist=[2,1,1,2,sqrt(5),sqrt(2),1,sqrt(2),sqrt(5),sqrt(5),2,sqrt(5),sqrt(5),sqrt(2),...
    1,sqrt(2),sqrt(5),sqrt(5),2,sqrt(5)];

waveC = zeros(size(depth));
ks=zeros(size(depth));
ks(1,1)=1; 
waveH = zeros(size(depth));

% calculate wave length (deep water)
tperiod0 = max(0.47*h0+6.76, pi*pi*sqrt(h0/grav));
L0 = grav*tperiod0^2*onpi2;
c0 = grav*tperiod0*onpi2;

waveL(1,1) = L0;
for j = 1: n_col
    for i = 1: n_row
        TM = L0;
        if(inland(i,j)==0)
            while   kk<10
                MN = 0.5*(waveL(i,j)+TM);
                TM = waveL(i,j);
                waveL(i,j) = L0*tanh(pi2*depth(i,j)/MN);
                if(abs(waveL(i,j)-TM)<1.0e-3 ) 
                    break
                end
            end 
            waveC(i,j) = c0*waveL(i,j)/L0;
            kh = depth(i,j)*pi2/waveL(i,j);
            tmp = 1+2*kh/sinh(2*kh);
            waveH(i,j) = h0/sqrt(tanh(kh)*tmp);
            ks(i,j) = sqrt(c0/(tmp*waveC(i,j))); 
        end
    end
end
travel = source;

keeploop = 1;
while  keeploop~=0
    keeploop = 0;
    for j = 1: n_col
        for i = 1: n_row
            if(travel(i,j)>=0)
                for k = 1:20
                    ix = i+iradius(k);
                    jx = j+jradius(k);
                    if(ix>0 && ix<=n_row && jx>0 && jx<=n_col)
                        if(inland(ix,jx)==1)
                            travel(ix,jx) = -1;
                        elseif(travel(ix,jx)<0)
                            travel(ix,jx) = travel(i,j)+dist(k)*dx/waveC(i,j);
                            keeploop = 1;
                            if(depth(i,j)/waveL(i,j)<0.5) 
                                frac = 2*(1-shoalC)*depth(i,j)/waveL(i,j)+shoalC;
                                if(waveH(ix,jx)>frac*waveH(i,j))
                                    waveH(ix,jx) = frac*waveH(i,j);
                                end
                            else
                                if(shadow==1 && waveH(ix,jx)>waveH(i,j))
                                    waveH(ix,jx) = waveH(i,j);
                                end
                            end
                        else
                            dt = travel(i,j)+(dist(k)*dx)/waveC(i,j);
                            if(travel(ix,jx)>dt && dt>0)
                                travel(ix,jx) = dt;
                                if(depth(i,j)/waveL(i,j)<0.5) 
                                    frac = 2*(1-shoalC)*depth(i,j)/waveL(i,j)+shoalC;
                                    if(waveH(ix,jx)>frac*waveH(i,j))
                                        waveH(ix,jx) = frac*waveH(i,j);
                                    end
                                else
                                    if(shadow==1 && waveH(ix,jx)>waveH(i,j))
                                        waveH(ix,jx) = waveH(i,j);
                                    end
                                end
                                keeploop = 1;
                            end 
                        end  
                    end  
                end 
            end 
        end 
    end 
    if(keeploop == 0) break
    end
end 

waveH = waveH.* ks;


if (dir==0)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i-1,j)>depth(i,j) && travel(i-1,j)>travel(i,j))
                waveH(i-1,j)=waveH(i-1,j)*constant_s;
            end
        end
    end
end

if (dir==90)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i,j+1)>depth(i,j) && travel(i,j+1)>travel(i,j))
                waveH(i,j+1)=waveH(i,j+1)*constant_s;
            end
        end
    end
end

if (dir==180)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i+1,j)>depth(i,j) && travel(i+1,j)>travel(i,j))
                waveH(i+1,j)=waveH(i+1,j)*constant_s;
            end
        end
    end
end

if (dir==270)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i,j-1)>depth(i,j) && travel(i,j-1)>travel(i,j))
                waveH(i,j-1)=waveH(i,j-1)*constant_s;
            end
        end
    end
end

if (dir==45)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i-1,j+1)>depth(i,j) && travel(i-1,j+1)>travel(i,j))
                waveH(i-1,j+1)=waveH(i-1,j+1)*constant_s;
            end
        end
    end
end

if (dir==135)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i+1,j+1)>depth(i,j) && travel(i+1,j+1)>travel(i,j))
                waveH(i+1,j+1)=waveH(i+1,j+1)*constant_s;
            end
        end
    end
end

if (dir==225)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i+1,j-1)>depth(i,j) && travel(i+1,j-1)>travel(i,j))
                waveH(i+1,j-1)=waveH(i+1,j-1)*constant_s;
            end
        end
    end
end



if (dir==315)
    for j = 2: n_col-1
        for i = 2: n_row-1
            if (depth(i-1,j-1)>depth(i,j) && travel(i-1,j-1)>travel(i,j))
                waveH(i-1,j-1)=waveH(i-1,j-1)*constant_s;
            end
        end
    end
end
