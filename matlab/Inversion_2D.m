% FWI inversion 

%%
% clear all
close all
%---------------------modelling set up-------------------------------------
%dt=0.0013;
fmax=40;
vmin = 1500;
vmax = 2500;
 NX=210;
 NZ=420;
%NX=400;
%NZ=240;
%NX=800;
%NZ=480;   %720 1200
%NX=20*3;
%NZ=120*3;

% Deduce dx and dt (in 2d)
% Dispersion condition (5 points per wavelength)
disp_z = vmin / fmax / 5;
disp_x = vmin / fmax / 5;

% Stabiblity condition
stability = disp_x / vmax / sqrt(2) * 0.80;

% choosing dt, dx,dz
dt= 0.0013/2;
dx=6/2;     %
dz=6/2;

if(dx>disp_x || dz>disp_z)
    disp('violation of the dispersion condition')
    return
end

if(dt>stability)
    disp('violation of the stability condition')
    return
end

rho=zeros(NZ+1,NX+1)+1000;
parameter.rho=rho;
parameter.dt=dt;
parameter.dz=dx;
parameter.dx=dz;

parameter.NX=NX;
parameter.NZ=NZ;
parameter.NSTEP=365*6;
NSTEP=365*6;

model=1;
if (model==1)
    % exact velocity model  
parameter.cs=zeros(NZ+1,NX+1);     
parameter.cs(1:floor(NZ/2),:)=1500;   
parameter.cs(floor(NZ/2)+1:NZ+1,:)=2500;

    % initial model 
parameter.cs1=zeros(NZ+1,NX+1);     
parameter.cs1(1:floor(NZ/2),:)=1500;   
parameter.cs1(floor(NZ/2)+1:NZ+1,:)=1500;
    
end

ax =(0:parameter.NX)*dx;
parameter.ax=ax;

az=(0:parameter.NZ)*dz;
parameter.az=az;
at =(0:parameter.NSTEP-1)*dt;


% density exact model
parameter.rho=zeros(NZ+1,NX+1)+1000;

% density starting model
parameter.rho1=zeros(NZ+1,NX+1)+1000;

%  exact stiffness model
parameter.mu=(parameter.rho.*parameter.cs.*parameter.cs);

%starting  stiffness  model
parameter.mu1=(parameter.rho1.*parameter.cs1.*parameter.cs1);

% Ricker source
[parameter.src,ns,ns2] = defsrc(fmax,dt);
parameter.model_exact=true;

zsrc = 3;
parameter.zsrc=2*3;  %12
parameter.xsrc=105;   %360
parameter.pmlx=80*2;    %(50 3)
parameter.pmlz=130*2;  % 140 (6) 240(3)
ixrec = (1:NX); % line of receivers
izrec = 40;   %240
parameter.ixrec=ixrec;
parameter.izrec=izrec;
parameter.matricial=true;
%--------------------------------------------------------------------------

% selecting the modelling scheme
parameter.choice = menu('Please choose a method','CONV 2','CONV 4', 'OPT 2', 'OPT 4','STAGG2','STAGG4') ;

%call to th modelling engine
tic
[ u_exact] = twoD_operators(parameter);
toc
%[ u_exact] = FWI_2D_mod(parameter );
%[A_con2,K_con2x, A_opt2,K_con4x, u_calc, value_mu_du_dxx] = FWI_2D_mod(parameter );
% parameter.model_exact=false;
% [ u_cal] = twoD_operators(parameter);
%[ u_cal] = FWI_2D_mod(parameter );

%observed and numerical data
dobs = squeeze(u_exact(izrec,ixrec,:));
%dobs = squeeze(u_exact_bis(izrec,ixrec,:));
% dcal= squeeze(u_cal(izrec,ixrec,:));
%return
%% Plotting

% Define the colormap
set(gcf,'Color','w');

cmap1 = colormap(bone);
cmap2 = colormap(hot);

cmap = zeros(64,3);
cmap(1:32,:)  = cmap1(1:2:end,:);
cmap(33:64,:) = cmap2(end:-2:1,:);
colormap(cmap); colorbar;

%   % adjoint source
%   adjoint=((dobs-dcal))';
%   adjoint=adjoint(NSTEP:-1:1,:);
%   parameter.adjoint=adjoint;
%   
%   % adjoint modelling
%   [u_adjoint] = ajoint_2D_mod(parameter);
%   u_adjoint=u_adjoint(:,:,NSTEP:-1:1);


% Display
%dmax = max(abs(dobs(:)))/5;
umax = max(abs(u_exact(:)))/5;
% umax1 = max(abs(u_adjoint(:)))/5;
for it = 2:NSTEP
%     subplot(1,2,1);
    imagesc(ax,az,squeeze(u_exact(:,:,it)));
    %title(sprintf('t[%d over %d] = %2.3f s',it,nt,ft + (it-1)*dt));
    axis image;
    %colorbar;
    set(gca,'FontSize',16);
    caxis([-umax,umax]);
    
%     subplot(1,2,2);
%     imagesc(ax,az,squeeze(u_adjoint(:,:,it)));
%     %title(sprintf('t[%d over %d] = %2.3f s',it,nt,ft + (it-1)*dt));
%     axis image;
%     %colorbar;
%     set(gca,'FontSize',16);
%     caxis([-umax1,umax1]);
       
%     subplot(1,2,2);
%     imagesc(ax,at(1:it),dobs(:,1:it)');
%     axis([ax(1),ax(end),at(1),at(end)]);
%     colorbar;
%     set(gca,'FontSize',16);
%      caxis([-dmax,dmax]);
    
   % colormap(cmap);
%    f=getframe;
drawnow;
% 
% G=sum(p_initial.*dlambda ,3);

end

 
%   % second order derivative of calculated wavefield
%   der_dcal= dcal_derivative( u_cal,dt,NZ,NX,NSTEP);
% 
%   % Gradient  
%   G= gradient_2D( u_adjoint,der_dcal);
%   
%   figure
%   imagesc(G)







