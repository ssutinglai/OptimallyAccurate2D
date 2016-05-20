function sig = defparamcpml(ax,vm,npml)

nx = length(ax);
dx = ax(2) - ax(1);

L = (npml-1)*dx;
R = 5.0;
vpml = 3 * vm / 2 / L * R;

nxe = nx + 2*npml;
sig = zeros(nxe,1);

for ind = 1:npml
    val = vpml*(ind/npml)^2;
    sig(npml - ind + 1)  = val;
    sig(nx + npml + ind) = val;
end