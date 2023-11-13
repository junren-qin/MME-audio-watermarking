function r = SDCVP(y, H, radius)
%CVP depth first
%for solving cvp problem
%the symbol set is real
%return the x coefficient, r

%   written by Shanxiang Lyu (s.lyu14@imperial.ac.uk)
%   Last updated on Mar. 14 /2017

if nargin == 2
    radius = realmax;
end

if size(H, 1) < size(H, 2)
	H = [H; zeros(size(H, 2) - size(H, 1), size(H, 2))];%make the matrix square
end

n = size(H,2);

[Q, R] = qr(H, 0);%economy version

z = Q'*y;



% add examine this variable before make it global
global SPHDEC_RADIUS;
global RETVAL;
global x;
global NUMX;
SPHDEC_RADIUS = radius;

RETVAL        = zeros(n, 1);
x        = zeros(n, 1);
NUMX    = 0;

sphdec_core(z, R, n, 0);

if NUMX > 0
    r = RETVAL;%retrival is set to r only this time
else
    r = zeros(n, 1);%no vector inside, return zeros
end


clear SPHDEC_RADIUS RETVAL SYMBSETSIZE NUMX;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sphdec_core(z, R, layer, dist)
%dist and d are both current radius
global SPHDEC_RADIUS;
global RETVAL;
global x;
global NUMX;
% disp(SPHDEC_RADIUS)
n=size(R,2);
           if layer==n
                 zi=z;%the last y is the first y
           else
                 zi=z - R(:,layer+1:end)*x(layer+1:end);%minius coef in prev rounds
           end
                c=round(zi(layer)/R(layer,layer));
                x(layer) = c;
              d = abs(z(layer) - R(layer,layer:end)*x(layer:end))^2 + dist;
                
                if (d <= SPHDEC_RADIUS) %while the current radius haven't exceed
                    if  layer==1
                 RETVAL        =  x;
                SPHDEC_RADIUS =  d;
                NUMX    =  NUMX + 1;
                    else
                    sphdec_core(z, R,  layer-1, d);
                    end
                end
                
                        delta=0;
 
        while  d<=SPHDEC_RADIUS % cumulated distance smaller than current outer radius
                delta=delta+1;
                for k=1:2
                    ci=c+delta*(-1)^k;
                        x(layer) = ci;
                        d = abs(z(layer) - R(layer,layer:end)*x(layer:end))^2 + dist;
                    if (d <= SPHDEC_RADIUS) %while the current radius haven't exceed
                        if  layer==1
                         RETVAL        =  x;
                        SPHDEC_RADIUS =  d;
                        NUMX    =  NUMX + 1;
%                         disp(RETVAL);
                        else
                        sphdec_core(z, R,  layer-1, d);
                        end
                    end
                    
                end
                
  
       
        end



