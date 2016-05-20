function [src,ns,ns2] = defsrc(fmax,dt)

fc = fmax / 2.5;
ns2 = floor(1+2*(4*7.5/fmax/2.5/dt+0.5));
ns = floor((ns2-1)/2);
src = (1:ns2)*0;

for is = -ns:ns
    a1 = is*fc*dt*pi;
    a2 = a1*a1;
    src(is+ns+1)  = (1-2*a2)*exp(-a2);
end

end