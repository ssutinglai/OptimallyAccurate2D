function [ u] = twoD_operators(parameter)

%%
% TAKING VARIABLES FROM THE DATA STRUCTURE

if (parameter.model_exact==true)
rho_ini=parameter.rho;
mu_ini=parameter.mu;
end

if (parameter.model_exact==false)
rho_ini=parameter.rho1;
mu_ini=parameter.mu1;
end


dx=parameter.dx;
dz=parameter.dz;
dz2=dz*dz;
dx2=dx*dx;
dt=parameter.dt;
dt2=dt*dt;
NX=parameter.NX+1;
NZ=parameter.NZ+1;
NSTEP=parameter.NSTEP;
choice=parameter.choice;
src=parameter.src;
nxext=parameter.pmlx;
nzext=parameter.pmlz;
xsrc=parameter.xsrc+nxext;
zsrc=parameter.zsrc+nzext;
matricial=parameter.matricial;
%
%%
% Extend the model mu------------------------------
nxe = NX + 2*nxext;
nze = NZ + 2*nzext;
% corners
mu = zeros(nze,nxe);
mu(1:nzext,1:nxext) = mu_ini(1,1); % upper left pml block
mu(1:nzext,nxe-nxext+1:nxe) = mu_ini(1,NX); % right upper block
mu(nze-nzext+1:nze,1:nxext) = mu_ini(NZ,1); % lower left block
mu(nze-nzext+1:nze,nxe-nxext+1:nxe) = mu_ini(NZ,NX); % lower right block

% edges
for iz = 1:nzext
    mu(iz,1+nxext:nxe-nxext) = mu_ini(1,:);
    mu(nze-nzext+iz,1+nxext:nxe-nxext) = mu_ini(NZ,:);
end

for ix = 1:nxext
    mu(1+nzext:nze-nzext,ix) = mu_ini(:,1);
    mu(1+nzext:nze-nzext,nxe-nxext+ix) = mu_ini(:,NX);
end
% inside
mu(nzext+1:nze-nzext,nxext+1:nxe-nxext) = mu_ini;

%%
% Extend the model rho-----------------------------
nxe = NX + 2*nxext;
nze = NZ + 2*nzext;
% corners
rho= zeros(nze,nxe);
rho(1:nzext,1:nxext) = rho_ini(1,1); % upper left pml block
rho(1:nzext,nxe-nxext+1:nxe) = rho_ini(1,NX); % right upper block
rho(nze-nzext+1:nze,1:nxext) = rho_ini(NZ,1); % lower left block
rho(nze-nzext+1:nze,nxe-nxext+1:nxe) = rho_ini(NZ,NX); % lower right block

% edges
for iz = 1:nzext
    rho(iz,1+nxext:nxe-nxext) = rho_ini(1,:);
    rho(nze-nzext+iz,1+nxext:nxe-nxext) = rho_ini(NZ,:);
end

for ix = 1:nxext
    rho(1+nzext:nze-nzext,ix) = rho_ini(:,1);
    rho(1+nzext:nze-nzext,nxe-nxext+ix) = rho_ini(:,NX);
end
% inside
rho(nzext+1:nze-nzext,nxext+1:nxe-nxext) = rho_ini;

%%
%pml
ax=parameter.ax;
az=parameter.az;
vm=(1500+2500)/2;
% For PML in x and z ----------------------------
[sigpmlxa] = defparamcpml1d(ax,vm,nxext);
[sigpmlza] = defparamcpml1d(az,vm,nzext);

%%
NX=nxe;
NZ=nze;
%% selecting the modelling operator
if choice ==1 ;
    CONV2 = true ;
    OPT2=false;
    disp('CONV2')
end

if choice==3
    OPT2=true;
    CONV2=false;
    disp('OPT2 choice')
end

%%
if(matricial==true)
% operator CONV2 and OPT2
A_con2=zeros(9,3,NZ,NX);
K_con2z=zeros(9,3,NZ,NX);
A_opt2=zeros(9,3,NZ,NX);
K_opt2z=zeros(9,3,NZ,NX);
K_opt2x=zeros(9,3,NZ,NX);
K_con2x=zeros(9,3,NZ,NX);


for i= 2:NX-1
 for j=2:NZ-1
    
    A_con2(1:9,1:3,j,i) = (1/dt2)*[0,0,0; 0, rho(j,i), 0; 0,0,0; ...
       0,0,0; 0,-2*rho(j,i),0;0,0,0;...
       0,0,0; 0,rho(j,i),0; 0,0,0];
    
    K_con2z(1:9,1:3,j,i) = (1/(dz2))*[0, 0,0; 0,0,0; 0,0,0;...
        0,(mu(j,i)+mu(j+1,i+1))/2,0; 0, -(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2,0; 0, (mu(j-1,i)+mu(j,i))/2,0 ;...
        0,0,0; 0,0,0,;0,0,0];
    
    A_opt2(1:9,1:3,j,i) = (1/(144*dt2))*[rho(j,i), 10*rho(j,i), rho(j,i); 10*rho(j,i), 100*rho(j,i),10*rho(j,i); rho(j,i), 10*rho(j,i), rho(j,i);...
        -2*rho(j,i), -20*rho(j,i),-2*rho(j,i); -20*rho(j,i), -200*rho(j,i), -20*rho(j,i); -2*rho(j,i), -20*rho(j,i), -2*rho(j,i);...
        rho(j,i), 10*rho(j,i), rho(j,i); 10*rho(j,i), 100*rho(j,i), 10*rho(j,i); rho(j,i), 10*rho(j,i), rho(j,i)];
    
    K_opt2z(1:9,1:3,j,i)= (1/(144*dz2))*[ (mu(j,i)+mu(j+1,i+1))/2, 10*(mu(j,i)+mu(j+1,i+1))/2, (mu(j,i)+mu(j+1,i+1))/2;...
         -(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2, -10*(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2,    -(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2;...
         (mu(j-1,i)+mu(j,i))/2, 10*(mu(j-1,i)+mu(j,i))/2, (mu(j-1,i)+mu(j,i))/2;...
        10*(mu(j,i)+mu(j+1,i+1))/2, 100*(mu(j,i)+mu(j+1,i+1))/2, 10*(mu(j,i)+mu(j+1,i+1))/2;...
         -10*(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2, -100*(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2, -10*(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2;...
         10*(mu(j-1,i)+mu(j,i))/2, 100*(mu(j-1,i)+mu(j,i))/2, 10*(mu(j-1,i)+mu(j,i))/2;...
         (mu(j,i)+mu(j+1,i+1))/2, 10*(mu(j,i)+mu(j+1,i+1))/2, (mu(j,i)+mu(j+1,i+1))/2;...
         -(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2, -10*(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2,-(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2;...
         (mu(j-1,i)+mu(j,i))/2, 10*(mu(j-1,i)+mu(j,i))/2, (mu(j-1,i)+mu(j,i))/2 ]; 
  
 end  
end


for j= 2:NZ-1
  for i= 2:NX-1
    
    K_con2x(1:9,1:3,j,i) = (1/(dx2))*[0, 0 ,0; 0,0,0; 0,0,0;...
        0,0,0;(mu(j,i-1)+mu(j,i))/2, -(mu(j,i-1)+2*mu(j,i)+mu(j,i+1))/2, (mu(j,i)+mu(j,i+1))/2 ; 0 0 0;...
        0,0,0;0,0,0;0,0,0 ];
    
    K_opt2x(1:9,1:3,j,i)= (1/(144*dx2))*[(mu(j,i)+mu(j+1,i+1))/2, 10*(mu(j,i)+mu(j+1,i+1))/2, (mu(j,i)+mu(j+1,i+1))/2;...
        -((mu(j-1,i)+2*mu(j,i))+mu(j+1,i+1)/2), -10*((mu(j-1,i)+2*mu(j,i))/2+mu(j+1,i+1))/2,  -((mu(j-1,i)+2*mu(j,i))/2+mu(j+1,i+1))/2;...
          (mu(j-1,i)+mu(j,i))/2, 10*(mu(j-1,i)+mu(j,i))/2, (mu(j-1,i)+mu(j,i))/2;...
          10*(mu(j,i)+mu(j+1,i+1))/2, 100*(mu(j,i)+mu(j+1,i+1))/2, 10*(mu(j,i)+mu(j+1,i+1))/2;...
        -10*(mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2, -100*((mu(j-1,i)+2*mu(j,i))/2+mu(j+1,i+1))/2,  -10*((mu(j-1,i)+2*mu(j,i)+mu(j+1,i+1))/2);...
          10*(mu(j-1,i)+mu(j,i))/2, 100*(mu(j-1,i)+mu(j,i))/2, 10*(mu(j-1,i)+mu(j,i))/2;...
          (mu(j,i)+mu(j+1,i+1))/2, 10*(mu(j,i)+mu(j+1,i+1))/2, (mu(j,i)+mu(j+1,i+1))/2;...
        -((mu(j-1,i)+2*mu(j,i))+mu(j+1,i+1)/2), -10*((mu(j-1,i)+2*mu(j,i))/2+mu(j+1,i+1))/2,  -((mu(j-1,i)+2*mu(j,i))/2+mu(j+1,i+1))/2;...
          (mu(j-1,i)+mu(j,i))/2, 10*(mu(j-1,i)+mu(j,i))/2, (mu(j-1,i)+mu(j,i))/2];
          

   end  
end

% left hand side corrector
K(:,:,:,:)=K_opt2x(:,:,:,:)+K_opt2z(:,:,:,:);
dAK(:,:,:,:)=A_opt2(:,:,:,:)-K(:,:,:,:);
end

%additional arrays
u=zeros(NZ,NX, NSTEP);
g=zeros(NZ,NX);
UI2=zeros(3,3, NZ,NX);
%apr=zeros(3,3);
%spr1=zeros(4,3);

if(matricial==false)
    disp('Non matricial computation')
    mu_dxx_left=zeros(NZ,NX);
    mu_dxx_right=zeros(NZ,NX);
    mu_dzz_left=zeros(NZ,NX);
    mu_dzz_right=zeros(NZ,NX);
    
    
    for i=2:NX-1
        mu_dxx_left(:,i)=(((mu(:,i-1)+mu(:,i))/2)./(rho(:,i)))*(dt2/dx2);
        mu_dxx_right(:,i)=(((mu(:,i)+mu(:,i+1))/2)./(rho(:,i)))*(dt2/dx2);
    end
    
    
    for j=2:NZ-1
        mu_dzz_left(j,:)=(((mu(j-1,:)+mu(j,:))/2)./(rho(j,:)))*(dt2/dz2);
    end
    
    for j=2:NZ-1
        for i=2:NX-1
        mu_dzz_right(j,i)=(((mu(j,i)+ mu(j+1,i+1))/2)./(rho(j,i)))*(dt2/dz2);
        end
    end
    
%        for j=2:NZ-1
%        % for i=2:NX-1
%         mu_dzz_right(j,:)=(((mu(j,:)+ mu(j+1,:))/2)./(rho(j,:)))*(dt2/dz2);
%     end
end

  value_mu_du_dxx=zeros(NZ,NX);
  value_mu_du_dzz=zeros(NZ,NX);
%%
for it = 2:NSTEP-1
     %%

%src=adjoint_source;
% CONV2
%--------------------------------------------------------------------------
 if(matricial==true)
if CONV2
      if(it==2)
      disp('CONV2 computation')
      end
      
       disp('iteration')
       fprintf('%d\n',it-1);
      
      value_mu_du_dxx(:,:)=0;
      value_mu_du_dzz(:,:)=0;

           
       %u(:) = 0.d0 ;
       % Predictor
       
         for i = 3:NX-2
             value_mu_du_dxx(:,i)=u(:,i-1,it).*squeeze(K_con2x(5,1,:,i))+u(:,i,it).*squeeze(K_con2x(5,2,:,i))+u(:,i+1,it).*squeeze(K_con2x(5,3,:,i));
         end
     
         for j = 30:NZ-100
             value_mu_du_dzz(j,:)=u(j-1,:,it).*squeeze(K_con2z(4,2,j,:))'+u(j,:,it).*squeeze(K_con2z(5,2,j,:))'+u(j+1,:,it).*squeeze(K_con2z(6,2,j,:))';
         end         
         
         
        % TIME INTEGRATION 
        for i=3:NX-2
            for j=30:NZ-100
                           
           u(j,i,it+1) = (value_mu_du_dxx(j,i)+ value_mu_du_dzz(j,i)-u(j,i,it)*A_con2(5,2,j,i)-u(j,i,it-1)*A_con2(8,2,j,i))/A_con2(2,2,j,i)...
               -sigpmlza(iz)* u(j,i,it)*dt-sigpmlxa(ix)* u(j,i,it)*dt;
              
            end
        end
        
         if(it <= size(src,2))
       u(zsrc,xsrc,it+1)=u(zsrc,xsrc,it+1)+src(1,it);
          end
     

 elseif OPT2
     
     if(it==2)
      disp('OPT2 computation')
     end
      disp('iteration');  fprintf('%d\n',it-1);
      
           %u(:) = 0.d0 ;
       % Predictor
       
         for i = 3:NX-2
             value_mu_du_dxx(:,i)=u(:,i-1,it).*squeeze(K_con2x(5,1,:,i))+u(:,i,it).*squeeze(K_con2x(5,2,:,i))+u(:,i+1,it).*squeeze(K_con2x(5,3,:,i));
         end
     
         for j = 10:NZ-250
             value_mu_du_dzz(j,:)=u(j-1,:,it).*squeeze(K_con2z(4,2,j,:))'+u(j,:,it).*squeeze(K_con2z(5,2,j,:))'+u(j+1,:,it).*squeeze(K_con2z(6,2,j,:))';
         end         
         
       
         
         
        % TIME INTEGRATION 
        for i=3:NX-2
            for j=10:NZ-250
           u(j,i,it+1) = (value_mu_du_dxx(j,i)+ value_mu_du_dzz(j,i)-u(j,i,it)*A_con2(5,2,j,i)-u(j,i,it-1)*A_con2(8,2,j,i))/A_con2(2,2,j,i)...
             -sigpmlza(iz)* u(j,i,it)*dt-sigpmlxa(ix)* u(j,i,it)*dt;
         
           
            end
        end
        
              
         UI2(:,:,:,:)=0;
           g(:,:)=0;
       for j = nzext+1: NZ-nzext
           for i=nxext+1:NX-nxext
           UI2(1:9,1:3,j,i) = [u(j+1,i-1,it+1) u(j+1,i,it+1) u(j+1,i+1,it+1);...%first 3 lines stand for time t+1
                                u(j,i-1,it+1) u(j,i,it+1) u(j,i+1,it+1);...          
                                u(j-1,i-1,it+1) u(j-1,i,it+1) u(j-1,i+1,it+1);
                                u(j+1,i-1,it) u(j+1,i,it) u(j+1,i+1,it);...    %next 3 lines stand for time t+1
                                u(j,i-1,it) u(j,i,it) u(j,i+1,it);...          
                                u(j-1,i-1,it) u(j-1,i,it) u(j-1,i+1,it);
                                u(j+1,i-1,it-1) u(j+1,i,it-1) u(j+1,i+1,it-1);...    %next 3 lines stand for time t+1
                                u(j,i-1,it-1) u(j,i,it-1) u(j,i+1,it-1);...          
                                u(j-1,i-1,it-1) u(j-1,i,it-1) u(j-1,i+1,it-1)];    
                            
                          if(i == xsrc) && (j == zsrc)
                                dirac = 1;
                          else
                               dirac = 0;
                          end
                          if(it <= size(src,2))
                           addsrc = src(1,it);
                          else
                           addsrc = 0;
                    
                          end
       
           g(j,i) = (-sum(sum((dAK(1:9,1:3,j,i)).*UI2(1:9,1:3,j,i)))/A_con2(2,2,j,i))+addsrc*dirac ;
           end
       end
        
       
         for i = nxext+1:NX-nxext
           for j=nzext+1:NZ-nzext
           u(j,i,it+1) = u(j,i,it+1)+(g(j,i)/A_con2(2,2,j,i)); %update: c = c0 + kdc
           end
         end      
        
        
         if(it <= size(src,2))
       u(zsrc,xsrc,it+1)=u(zsrc,xsrc,it+1)+src(1,it);
         end

end
end
   
%     
%       if(matricial == false)
%           
%           %PREDICTED WAVEFIELD
%           if CONV2
%            if(it==2)
%              disp('CONV2 computation')
%            end 
%            
%           for i = 2:NX-1
%             for j = 2:NZ-1
% 
%                   
%            u(j,i,it+1) = 2*u(j,i,it)-u(j,i,it-1)+ mu_dxx_left(j,i)*(u(j,i-1,it)-u(j,i,it))-mu_dxx_right(j,i)*(u(j,i,it)-u(j,i+1,it))+...
%                 mu_dzz_left(j,i)*(u(j-1,i,it)-u(j,i,it))- mu_dzz_right(j,i)*(u(j,i,it)-u(j+1,i,it));
%             
%               if(it <= size(src,2))
%                   if(i == xsrc) && (j == zsrc)
%                    u(zsrc,xsrc,it+1)=u(zsrc,xsrc,it+1)+src(1,it);
%                   end
%              end
%             end
%           end
%          
%           
%           
%   elseif OPT2
%      if(it==2)
%       disp('OPT2 computation')
%      end 
%       % PREDICTED WAVEFIELD
%           for i = 2:NX-1
%             for j = 2:NZ-1
% 
%            u(j,i,it+1) = 2*u(j,i,it)-u(j,i,it-1)+ mu_dxx_left(j,i)*(u(j,i-1,it)-u(j,i,it))-mu_dxx_right(j,i)*(u(j,i,it)-u(j,i+1,it))+...
%                 mu_dzz_left(j,i)*(u(j-1,i,it)-u(j,i,it))- mu_dzz_right(j,i)*(u(j,i,it)-u(j+1,i,it));
%             
%               if(it <= size(src,2))
%                   if(i == xsrc) && (j == zsrc)
%                    u(zsrc,xsrc,it+1)=u(zsrc,xsrc,it+1)+src(1,it);
%                   end
%              end
%             end
%           end
%           
%           
%           
%           %CORRECTED WAVEFIELD
%            for i = 2:NX-2
%             for j = 2:NZ-1
% 
%      apr(1:10,1) = [u(j,i,it+1)-2*u(j,i,it)+u(j,i,it-1);...             %a*p'r'  1
%                    u(j-1,i-1,it+1)-2*u(j-1,i-1,it)+u(j-1,i-1,it-1);...  %a*(p'-1)(r'-1) 2
%                    u(j+1,i-1,it+1)-2*u(j+1,i-1,it)+u(j+1,i-1,it-1);...  %a*(p'-1)(r'+1) 3
%                    u(j-1,i+1,it+1)-2*u(j-1,i+1,it)+u(j-1,i+1,it-1); ... %a*(p'+1)(r'-1) 4
%                    u(j+1,i+1,it+1)-2*u(j+1,i+1,it)+u(j+1,i+1,it-1);...  %a*(p'+1)(r'+1)  5
%                    u(j,i-1,it+1)-2*u(j,i-1,it)+u(j,i-1,it-1);...        %a*(p'-1)r'          6
%                    u(j-1,i,it+1)-2*u(j-1,i,it)+u(j-1,i,it-1);...        %a*p'(r'-1)      7
%                    u(j+1,i,it+1)-2*u(j+1,i,it)+u(j+1,i,it-1);...        %a*p'(r'+1)'       8
%                    u(j,i+1,it+1)-2*u(j,i+1,it)+u(j,i+1,it-1);...        %a*(p'+1)r'           9
%                    u(j+1,i+2,it+1)-2*u(j+1,i+2,it)+u(j+1,i+2,it-1)];    %a*(p'+2)(r'+1)    10
%                
%                
%               Cp(1:10,1) = [apr(1,1)+12*u(j,i,it);...      %C*p'r'=a*p'r'+12Cp'r'   1
%                          apr(6,1)+12*u(j,i-1,it);...       %C*(p'-1)r'=a*(p'-1)r'+12C(p'-1)r'   2
%                          apr(2,1)+12*u(j-1,i-1,it);...     %C*(p'-1)(r'-1)=a*(p'-1)(r'-1)+12*C(p'-1)(r'-1) 3
%                          apr(7,1)+12*u(j-1,i,it);...       %C*p'(r'-1)=a*p'(r'-1)+12Cp'(r'-1)  4
%                          apr(8,1)+12*u(j+1,i,it);...       %C*p'(r+1)=a*p'(r'+1)+12*Cp'(r'+1)  5
%                          apr(5,1)+12*u(j+1,i+1,it);...     %C*(p'+1)(r'+1)=a*(p'+1)(r'+1)+12*C(p'+1)(r'+1)6
%                          apr(4,1)+12*u(j-1,i+1,it);...     %C*(p'+1)(r'-1)=a*(p'+1)(r'-1)+12Cp(r-1)7
%                          apr(9,1)+12*u(j,i+1,it);...       %C*(p'+1)r'=a*(p'+1)r'+12*C(p'+1)r' 8
%                          apr(10,1)+12*u(j+1,i+2);...       %C*(p'+2)(r'+1)=a*(p'+2)(r'+1)+12*C(p'+2)(r'+1) 9
%                          apr(3,1)+12*u(j+1,i-1,it)];       %C*(p'-1)(r'+1)=a*(p'-1)(r'+1)+12*C(p'-1)(r'+1) 10
%                      
%                      
%               S1_pr(:,1) = [Cp(2,1)-Cp(1,1);...           %S1*p'r'=C*(p'-1)r'-C*p'r'
%                               Cp(3,1)-Cp(4,1);...         %S1*p'(r'-1)=C*(p'-1)(r'-1)-C*p'(r'-1) 
%                               Cp(5,1)-Cp(6,1)];           %S1*(p'+1)(r'+1)=C*p'(r'+1)-C*(p'+1)(r'+1) 
%                        
%               S3_pr(:,1) = [Cp(4,1)-Cp(1,1);...           %S3*p'r'=C*p'(r'-1)-C*p'r'
%                               Cp(3,1)-Cp(2,1);...         %S3*(p'-1)r'=C*(p'-1)(r'-1)-C*(p'-1)r'
%                               Cp(7,1)-Cp(8,1)];           %S3*(p'+1)r'=C*(p'+1)(r'-1)-C*(p'+1)r'
%                           
%               S1p1r(:,1) = [Cp(4,1)-Cp(7,1);...           %s1*(p'+1)(r'-1)=C*p'(r'-1)-C*(p'+1)(r'-1)
%                               Cp(1,1)-Cp(8,1);...         %s1*(p'+1)r'=C*p'r'-C*(p'+1)r'
%                               Cp(4,1)-Cp(7,1)] ;          %s1*(p'+2)(r'+1)=C*(p'+1)(r'+1)-C*(p'+2)(r'+1)
%                                                            
%                                                            
%                                                            
%              S3p1r(:, 1) = [Cp(2,1)-Cp(10,1);...           %s3*(p'-1)(r'+1)=C*(p'-1)r'-C*(p'-1)(r'+1)
%                                Cp(1,1)-Cp(5,1);...        %s3*p'(r'+1)=C*p'r'-C*p'(r'+1)
%                                Cp(8,1)-Cp(6,1)];          %s3*(p'+1)(r'+1)=C*(p'+1)r'-C*(p'+1)(r'+1)
%                           
%                      
%                      
%                   S1_pr_dble_aster  =  S1_pr(2,1)+10*S1_pr(1,1)+S1_pr(3,1);
%                   S3_pr_dble_aster  =   S3_pr(2,1)+10*S3_pr(1,1)+S3_pr(3,1);
%                   S1_pr1_dble_aster =  S1p1r(1,1)+10*S1p1r(2,1)+S1p1r(3,1);
%                   S3_pr1_dble_aster =  S3p1r(1,1)+10*S3p1r(2,1)+S3p1r(3,1);
%                    
%                    
%               
%              u(j,i,it+1) = u(j,i,it+1)+(-1/144)*((apr(2,1)+apr(3,1)+apr(4,1)+apr(5,1))...
%                           +10*(apr(6,1)+apr(7,1)+apr(8,1)+apr(9,1))+100*apr(1,1))...
%                           +(mu_dxx_left(j,i)/144)*S1_pr_dble_aster-(mu_dxx_right(j,i)/144)*S1_pr1_dble_aster...
%                           +(mu_dzz_left(j,i)/144)*S3_pr_dble_aster-(mu_dzz_right(j,i))*S3_pr1_dble_aster;
%                    
%                      if(it <= size(src,2))
%                         if(i == xsrc) && (j == zsrc)
%                            u(zsrc,xsrc,it+1) = u(zsrc,xsrc,it+1)+src(1,it);
%                         end
%                      end
%                                                                                                                  
%          
%                
%             end
%            end
%               
%           end      
      
%       end
      
      
end

u=u(nzext:nze-nzext,nxext:nxe-nxext,:);
%save ('test','u','-v7.3')
end
























