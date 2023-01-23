% The code is the code of Fast Implementation of GEP-LMR-S-FDTD-II for solving 2D Maxwell Eq using local space mesh refinement
% which preserves global energy;
% 
% electric field:  E_x, E_y
% magnetic filed: H_z
% induced magnetic current: K_z
% induced electric current: J_x, J_y
%
% E , H , K, J are energy conserved in this version
clear; 

tic

format long
%% domain is [xa,xb]x[ya,yb]
xa=0; xb=1;     
yc=0; yd=1;

%% parameters
mu=1;     %% constant  magnetic permeability 
epsilon=1; %% constant  electric permittivity
omega_pe = 1; %% electric plasma frequency
omega_pm = omega_pe; %% magnetic plasma frequency

%% metamaterial parameters
re = (epsilon*omega_pe^2)^(-1);
gamma_m = 1; %% magnetic damping frequency
rm = (mu*omega_pm^2)^(-1);
gamma_e = 1; %% electric damping frequency

kx=1;                                   %% wave number
ky=1;
omg=abs(sqrt((kx^2+ky^2)/epsilon/mu));  %% angular frequency of the wave

T=1;            %% the end time
h1=1/250;        %% Coarse Space step size in X direction   
h2=h1/2;        %% Fine Space step size in X direction  
hyc=h1;         %% space stepsize in y direction on the coarse 
hyf=h2;         %% stepsize in y on the fine grid

tau=1/500;       %% Time Space Step size with step size tau = 1/640
N=T/tau;        %% Maximum number of time steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for meta-refine-FDTD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pa = (rm/tau + rm*gamma_m*.5)^-1;
pb = (rm/tau - rm*gamma_m*.5);
pc = ((mu/tau) + .25*pa)^-1;
pd = ((mu/tau) - .25*pa);
pe = (re/tau + re*gamma_e*.5)^-1;
pf = (re/tau - re*gamma_e*.5);
pex =(re/tau + re*gamma_e*.25)^-1;
pfx =(re/tau - re*gamma_e*.25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t=0:tau:T;
X_c=0.5;                                                     %% Interface
Mxc=(X_c-xa)/2/h2;                                           %% grid points on the coarse grid on x direction
Mxf=(xb-(X_c))/h2;                                           %% grid points on the fine grid on x direction
x=[xa:2*h2:X_c (X_c+h2):h2:xb];                              %% integer points on x direction
xx=[(xa+0.5*h1):h1:(X_c-h2) (X_c+0.5*h2):h2:(xb-0.5*h2)];    %% half points on the x direction E_x, H_z location in x direction;
Myc=(yd-yc)/hyc;                                             %% grid points on the coarse grid on y direction
Myf=(yd-yc)/hyf;                                             %% grid points on the fine grid on y direction
y1=yc:2*h2:yd;                                               %% grid points in y direction on the coarse grid
y2=yc:h2:yd;                                                 %% grid points in y direction on the fine grid
y3=(yc+0.5*h1):h1:(yd-h2);                                   %% half points on the coarse grid in y diretion       
y4=[(yc+0.5*h2):h2:(yd-0.5*h2)];                             %% half points on the fine grid in y direction;                                     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize partial energy W_c1 and W_f1
% in x direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_c1 = zeros(1,Mxc); 
W_f1 = zeros(1,Mxc+Mxf+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize partial energy W_c2 and W_f2
% in y direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_f2 = zeros(1,Mxc+Mxf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize partial energy W_c3 and W_f3
% at interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_c3 = 0; 
W_f3 = 0;



%%%%%
% with metameterial 
%% coefficient matrix in x direction
      a1(1:Mxc-2) = 1 + pc*tau/(8*h2^2*epsilon) + pe*tau/(4*epsilon);
      a2 = 1 + pc*tau/(8*h2^2*epsilon) + pe*tau/(4*epsilon);
      a3= 3*pc/(2*h2) + 6*epsilon*h2/tau + 1.5*pe*h2;
      a4(1:Mxf-1)=(1 + pc*tau/(2*h2^2*epsilon) + pe*tau/(4*epsilon))*2;

     a = [a1 a2 a3 a4];

      b1(1:Mxc-3) = -pc*tau/(16*h2^2*epsilon);
      b2 = -pc*tau/(16*h2^2*epsilon);
      b3 = -pc/(2*h2);
      b4(1:Mxf-1)= -pc*tau/(2*h2^2*epsilon);

     b = [b1 b2 b3 b4];
    
    

     c1(1:(Mxc-2)) = -pc*tau/(16*h2^2*epsilon);
     c2 = -pc*tau/(16*h2^2*epsilon);
     c3 = -pc/h2;
     c4(1:(Mxf-2))= -pc*tau/(2*h2^2*epsilon);

     c = [c1 c2 c3 c4];
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% coefficient matrix in x direction on the fine grid by solving "U_form structure"

a11=[ones(1,Mxf-1)*(1+pc*tau/(2*h2^2*epsilon)+pe*tau/(4*epsilon)) (-1)*(epsilon*h2/tau+pc/(2*h2)+pe*h2/4) 2 (-1)*(epsilon*h2/tau+pc/(2*h2)+pe*h2/4) ones(1,Mxf-1)*(1+pc*tau/(2*h2^2*epsilon)+pe*tau/(4*epsilon))];
a22=[ones(1,Mxf-2)*(-1)*(pc*tau/(4*h2^2*epsilon))  pc/(2*h2)  -1  -2*epsilon*h2/tau-pc/(4*h2)-pe*h2/2 ones(1,Mxf-1)*(-1)*(pc*tau/(4*h2^2*epsilon))];
a33=[ones(1,Mxf-1)*(-1)*(pc*tau/(4*h2^2*epsilon))  -2*epsilon*h2/tau-pc/(4*h2)-pe*h2/2  -1 pc/(2*h2) ones(1,Mxf-2)*(-1)*(pc*tau/(4*h2^2*epsilon))];






%% coefficient matrix in y direction

for j=1:Myc-1
    aa(j)=1+tau^2/32/mu/epsilon/h2^2 + tau*pex/16/epsilon;
end

for j=1:Myc-2
    bb(j)=-tau^2/64/mu/epsilon/h2^2;
end

for j=1:Myc-2
    cc(j)=-tau^2/64/mu/epsilon/h2^2;
end


%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:Myf-1
    aaa(j)=1+tau^2/8/mu/epsilon/h2^2 + tau*pex/16/epsilon;
end

for j=1:Myf-2
    bbb(j)=-tau^2/16/mu/epsilon/h2^2;
end

for j=1:Myf-2
    ccc(j)=-tau^2/16/mu/epsilon/h2^2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initialize variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize variables
E_xc(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
E_xf(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);
E_yc(1:Mxc+1,1:Myc) = zeros(Mxc+1,Myc);
E_yf(Mxc+1:Mxc+Mxf+1,1:Myf) = zeros(Mxf+1,Myf);

E_yf1(Mxc+1:Mxc+Mxf+1,1:Myc) = zeros(Mxf+1,Myc);
J_yf1(Mxc+1:Mxc+Mxf+1,1:Myc) = zeros(Mxf+1,Myc);

H_zc(1:Mxc,1:Myc) = zeros(Mxc,Myc);
H_zf(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);

J_xc(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
J_xf(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);
J_yc(1:Mxc+1,1:Myc) = zeros(Mxc+1,Myc);
J_yf(Mxc+1:Mxc+Mxf+1,1:Myf) = zeros(Mxf+1,Myf);
K_zc(1:Mxc,1:Myc) = zeros(Mxc,Myc);
K_zf(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);

H_z1f(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);
H_zf1(Mxc+1:Mxc+Mxf,1:Myc) = zeros(Mxf,Myc);
K_zf1(Mxc+1:Mxc+Mxf,1:Myc) = zeros(Mxf,Myc);

H_z1c(1:Mxc,1:Myc) = zeros(Mxc,Myc);
H_zc1(1:Mxc,1:Myc) = zeros(Mxc,Myc);

E_xc_next(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
E_xf_next(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);
E_yc_next(1:Mxc+1,1:Myc) = zeros(Mxc+1,Myc);
E_yf_next(Mxc+1:Mxc+Mxf+1,1:Myf) = zeros(Mxf+1,Myf);
H_zc_next(1:Mxc,1:Myc) = zeros(Mxc,Myc);
H_zf_next(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);

J_xc_next(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
J_xf_next(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);
J_yc_next(1:Mxc+1,1:Myc) = zeros(Mxc+1,Myc);
J_yf_next(Mxc+1:Mxc+Mxf+1,1:Myf) = zeros(Mxf+1,Myf);
K_zc_next(1:Mxc,1:Myc) = zeros(Mxc,Myc);
K_zf_next(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);

H_z1c_next(1:Mxc,1:Myc) = zeros(Mxc,Myc);
H_z1f_next(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);
H_zc1_next(1:Mxc,1:Myc) = zeros(Mxc,Myc);

%% 3 steps method intermediates
E_x1c(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
E_x1c_next(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
E_x1f(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);
E_x1f_next(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);
J_x1c(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
J_x1c_next(1:Mxc+1,1:Myc+1) = zeros(Mxc+1,Myc+1);
J_x1f(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);
J_x1f_next(Mxc+1:Mxc+Mxf+1,1:Myf+1) = zeros(Mxf+1,Myf+1);

H_zf2(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);
H_z2c(1:Mxc,1:Myc) = zeros(Mxc,Myc);
H_z2f(Mxc+1:Mxc+Mxf,1:Myf) = zeros(Mxf,Myf);
%%%% boundary condition

    for i=1:Mxc
        E_xc(i,1)=0; 
        E_xc(i,Myc+1)=0;  
        E_x1c(i,1)=0; 
        E_x1c(i,Myc+1)=0; 
    end

%%%%%%%%%%%%

    for i=Mxc+1:Mxc+Mxf
            E_xf(i,1)=0;
            E_xf(i,2*Myc+1)=0;
            E_x1f(i,1)=0;
            E_x1f(i,2*Myc+1)=0;
    end


    for j=1:Myc
           E_yc(1,j)=0;            
    end

    for j=1:2*Myc
           E_yf(Mxc+Mxf+1,j)=0;
    end

    for i=1:Mxc
        J_xc(i,1)=0; 
        J_xc(i,Myc+1)=0;
        J_x1c(i,1)=0; 
        J_x1c(i,Myc+1)=0;
    end

    for i=Mxc+1:Mxc+Mxf
            J_xf(i,1)=0;
            J_xf(i,2*Myc+1)=0;
            J_x1f(i,1)=0;
            J_x1f(i,2*Myc+1)=0;
    end

    for j=1:Myc
           J_yc(1,j)=0;            
    end
      
    for j=1:2*Myc
           J_yf(Mxc+Mxf+1,j)=0;
    end

% initial value for  E_x
for i=1:Mxc
    for j=1:Myc+1
        E_xc(i,j)=ky/epsilon/sqrt(mu)/omg*cos(kx*pi*xx(i))*sin(ky*pi*y1(j));       
    end
end

for i=Mxc+1:Mxc+Mxf
    for j=1:2*Myc+1
        E_xf(i,j)=ky/epsilon/sqrt(mu)/omg*cos(kx*pi*xx(i))*sin(ky*pi*y2(j));
    end
end
%%%%%%%%%%% initial value for E_y
   for i=1:Mxc+1
       for j=1:Myc
           E_yc(i,j)=-kx/epsilon/sqrt(mu)/omg*sin(kx*pi*x(i))*cos(ky*pi*y3(j));
       end
   end
        
      for i=Mxc+1:Mxc+Mxf+1
          for j=1:2*Myc
              E_yf(i,j)=-kx/epsilon/sqrt(mu)/omg*sin(kx*pi*x(i))*cos(ky*pi*y4(j));
          end
      end
      %%%%%%%% initial value for H_z
      
      for i=1:Mxc
          for j=1:Myc
              H_zc(i,j)=0;
          end
      end
      
      
      for i=Mxc+1:Mxc+Mxf
          for j=1:2*Myc
              H_zf(i,j)=0;
          end
      end
%%%%%%%%% initial value for J x direction
for i=1:Mxc
    for j=1:Myc+1
         J_xc(i,j)= 0;
    end
end

for i=Mxc+1:Mxc+Mxf
    for j=1:2*Myc+1
         J_xf(i,j)= 0;
    end
end   
%%%%%%%%% initial value for J y direction
   for i=1:Mxc+1
       for j=1:Myc
            J_yc(i,j) = 0;
       end
   end
        
   for i=Mxc+1:Mxc+Mxf+1
       for j=1:2*Myc
            J_yf(i,j) = 0;
       end
   end
%%%%%%%%% initial value for K_z

   for i=1:Mxc
       for j=1:Myc
            K_zc(i,j) = 0;
       end
   end
      
      
   for i=Mxc+1:Mxc+Mxf
        for j=1:2*Myc
             K_zf(i,j) = 0;
        end
   end
      
%% initial value of auxiliary  variable
% H E K J variables used in x direction matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=Mxc+1:Mxc+Mxf
        for j=1:Myc
         H_zf1(i,j)=(H_zf(i,2*j)+H_zf(i,2*j-1))/2;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K avg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=Mxc+1:Mxc+Mxf
        for j=1:Myc
         K_zf1(i,j)=(K_zf(i,2*j)+K_zf(i,2*j-1))/2;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E avg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=Mxc+1:Mxc+Mxf+1
       for j=1:Myc
         E_yf1(i,j)=(E_yf(i,2*j)+E_yf(i,2*j-1))/2;
       end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J avg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=Mxc+1:Mxc+Mxf+1
       for j=1:Myc
         J_yf1(i,j)=(J_yf(i,2*j)+J_yf(i,2*j-1))/2;
       end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
%% how to define  the sum variables on the fine grid .......
%%%%% numerical solutions on the coarse grid: E_yc, E_xc, H_zc, H1_zc( intermediate variable)

  for k=1:N
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  solve in y direction   stage 1        
       for i=1:Mxc+Mxf
           if i>=1&&i<=Mxc
          RHS0(1)= E_xc(i,2)+tau/4/epsilon/h2*(H_zc(i,2)- H_zc(i,1) )...
              +tau^2/64/epsilon/mu/h2^2*(E_xc(i,3)-2*E_xc(i,2)+E_xc(i,1))...
              -tau*(pex*pfx+1)/4/epsilon*J_xc(i,2)...
              -tau*pex/16/epsilon*E_xc(i,2);
           for j=3:Myc-1
               RHS0(j-1)= E_xc(i,j)+tau/4/epsilon/h2*(H_zc(i,j)- H_zc(i,j-1))...
                   +tau^2/64/epsilon/mu/h2^2*(E_xc(i,j+1)-2*E_xc(i,j)+E_xc(i,j-1))...
                   -tau*(pex*pfx+1)/4/epsilon*J_xc(i,j)...
                   -tau*pex/16/epsilon*E_xc(i,j);
           end
           RHS0(Myc-1)= E_xc(i,Myc)+tau/4/epsilon/h2*(H_zc(i,Myc)- H_zc(i,Myc-1))...
               +tau^2/64/epsilon/mu/h2^2*(E_xc(i,Myc+1)-2*E_xc(i,Myc)+E_xc(i,Myc-1))...
               -tau*(pex*pfx+1)/4/epsilon*J_xc(i,Myc)...
               -tau*pex/16/epsilon*E_xc(i,Myc);

           X = chase(aa,bb,cc,RHS0);
           for j=2:Myc
            E_x1c_next(i,j)=X(j-1);
           end
         
           %%%%%% solve explicitly H_z^*

           for j=1:Myc
               H_z1c(i,j)=H_zc(i,j)+...
                   tau/8/mu/h2*(E_x1c_next(i,j+1)...
                   -E_x1c_next(i,j) + E_xc(i,j+1)-E_xc(i,j));
           end

           %%%%%%% solve explicitly J_x^*
           for j=2:Myc
                J_x1c_next(i,j)=pex*pfx*J_xc(i,j)...
                 + pex/4*(E_x1c_next(i,j) + E_xc(i,j));
           end

           %%%%%%%%%%%%%%%%%%%%%%%
           
           
           %%%%%%%
           elseif i>=Mxc+1&&i<=(Mxc+Mxf)
               RHS1(1)=E_xf(i,2)+tau/2/epsilon/h2*(H_zf(i,2)-H_zf(i,1))...
               +tau^2/16/epsilon/mu/h2^2*(E_xf(i,3)-2*E_xf(i,2)+E_xf(i,1))...
               -tau*(pex*pfx+1)/4/epsilon*J_xf(i,2)...
               -tau*pex/16/epsilon*E_xf(i,2);
              
               for j=3:Myf-1
               RHS1(j-1)=E_xf(i,j)+tau/2/epsilon/h2*(H_zf(i,j)-H_zf(i,j-1))...
                      +tau^2/16/epsilon/mu/h2^2*(E_xf(i,j+1)-2*E_xf(i,j)+E_xf(i,j-1))...
                      -tau*(pex*pfx+1)/4/epsilon*J_xf(i,j)...
                      -tau*pex/16/epsilon*E_xf(i,j);
               end
               RHS1(Myf-1)=E_xf(i,Myf)+tau/2/epsilon/h2*(H_zf(i,Myf)-H_zf(i,Myf-1))...
                   +tau^2/16/epsilon/mu/h2^2*(E_xf(i,Myf+1)-2*E_xf(i,Myf)+E_xf(i,Myf-1))...
                   -tau*(pex*pfx+1)/4/epsilon*J_xf(i,Myf)...
                   -tau*pex/16/epsilon*E_xf(i,Myf);


              XX=chase(aaa,bbb,ccc,RHS1);
              
              for j=2:Myf
              E_x1f_next(i,j)=XX(j-1);
              end

               %%%%%% solve explicitly H_z1f
           for j=1:Myf
               H_z1f(i,j)=H_zf(i,j)...
                   +tau/4/mu/h2*(E_x1f_next(i,j+1)-E_x1f_next(i,j)...
                   + E_xf(i,j+1)-E_xf(i,j) );
           end
               %%%%%% solve explicitly J_x1f
           for j=2:Myf
                J_x1f_next(i,j)=pex*pfx*J_xf(i,j)...
                 + pex/4*(E_x1f_next(i,j) + E_xf(i,j));
           end

           end
       end
       
% %    %%%%%%%%%% the sum of value on the fine grid;   \bar{H_z^f}^{*}_{i,2j+1}  

        for i=Mxc+1:Mxc+Mxf
            for j=1:Myc
                 H_zf2(i,j)=(H_z1f(i,2*j)+H_z1f(i,2*j-1))/2;
            end
        end
          
      
      
      
      
      
      
      
      
      
      
      
      
      
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Metamaterial with Js and Ks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for j=1:Myc          
             for i=2:Mxc-1
                RHS2(i-1)=(1- (pe*tau)/(4*epsilon))*E_yc(i,j)...
                -tau*(pc*pd+1)/4/epsilon/h2*(H_z1c(i,j)-H_z1c(i-1,j))...
                +tau*pc/16/epsilon/h2^2*(E_yc(i+1,j)-2*E_yc(i,j)+E_yc(i-1,j))...
                -tau*(pe*pf+1)/2/epsilon*J_yc(i,j)...
                +tau*(pa*pb+1)*pc/8/epsilon/h2*(K_zc(i,j) - K_zc(i-1,j));
             end    


 RHS2(Mxc-1)=(1- (pe*tau)/(4*epsilon))*E_yc(Mxc,j)...
     -tau*(pc*pd+1)/4/epsilon/h2*(H_z1c(Mxc,j)-H_z1c(Mxc-1,j) )...
     +tau/16/epsilon/h2^2*pc*(E_yc(Mxc+1,j)-2*E_yc(Mxc,j)+E_yc(Mxc-1,j) )...
     +tau/16/epsilon/h2^2*pc*(E_yf1(Mxc+1,j)-E_yc(Mxc+1,j))...
     -tau*(pe*pf+1)/2/epsilon*J_yc(Mxc,j)...
     +tau*(pa*pb+1)*pc/8/epsilon/h2*(K_zc(Mxc,j) - K_zc(Mxc-1,j));
 
     %%%% x_0    interface line
     
      RHS2(Mxc)=-(2*epsilon*h2/tau + 3*pc/(2*h2)+pe*h2*1.5)*E_yf1(Mxc+1,j)...
           +pc/(2*h2)*E_yc(Mxc,j) + 8*h2/tau*E_yc(Mxc+1,j)...
           +pc/h2*E_yf1(Mxc+2,j)...
           -2*(pc*pd + 1)*(H_zf2(Mxc+1,j) - H_z1c(Mxc,j))...
           -1*h2*(pe*pf+1)*(J_yf1(Mxc+1,j)+2*J_yc(Mxc+1,j))...
           +pc*(pa*pb+1)*(K_zf1(Mxc+1,j) - K_zc(Mxc,j));

     for i=Mxc+2:Mxc+Mxf
            RHS2(i-1)=2*(1 - (pe*tau)/(4*epsilon))*E_yf1(i,j)...
            -2*tau*(pc*pd+1)/2/epsilon/h2*(H_zf2(i,j)-H_zf2(i-1,j))...
            +2*tau/4/epsilon/h2^2*pc*(E_yf1(i+1,j)-2*E_yf1(i,j)+E_yf1(i-1,j))...
            -2*tau*(pe*pf+1)/2/epsilon*J_yf1(i,j)...
            +2*tau*(pa*pb+1)*pc/4/epsilon/h2*(K_zf1(i,j) - K_zf1(i-1,j));
     end 
 
       XXX=chase(a,b,c,RHS2);

            E_yc_next(2:Mxc,j)=XXX(1:Mxc-1);


       for i=1:Mxc-1
             H_z2c(i,j) = pc*pd*H_z1c(i,j)...
             -pc/(4*h2)*(E_yc_next(i+1,j) - E_yc_next(i,j))...
             -pc/(4*h2)*(E_yc(i+1,j) - E_yc(i,j))...
             -pc*(pa*pb+1)/2*K_zc(i,j);
       end
       
 %%%%%% solve J_yc and K_zc explicitly
        for i=1:Mxc-1
           K_zc_next(i,j) = pa*pb*K_zc(i,j)...
             + (pa/2)*(H_z2c(i,j) + H_z1c(i,j));
       end
       

       for i=2:Mxc
           J_yc_next(i,j) = pe*pf*J_yc(i,j)... 
              + pe/2*(E_yc_next(i,j) + E_yc(i,j));
       end
 
        for i=Mxc+Mxf:-1:Mxc+2
             RHS4(Mxc+Mxf+1-i)= (1- (pe*tau)/(4*epsilon))*E_yf(i,2*j)...
             -tau*(pc*pd+1)/2/epsilon/h2*(H_z1f(i,2*j)-H_z1f(i-1,2*j))...
             +tau/4/epsilon/h2^2*pc*(E_yf(i+1,2*j)-2*E_yf(i,2*j)+E_yf(i-1,2*j))...
             -tau*(pe*pf+1)/2/epsilon*J_yf(i,2*j)...
             +tau*(pa*pb+1)*pc/4/epsilon/h2*(K_zf(i,2*j)- K_zf(i-1,2*j));
        end


        RHS4(Mxf)=...
        -pc/(4*h2)*E_yc_next(Mxc,j) + (-2*epsilon*h2/tau + pc/(4*h2) + pe*h2/2)*E_yc(Mxc+1,j)...
        -pc/(4*h2)*E_yc(Mxc,j)...
        -(epsilon*h2/tau - pc/(2*h2) - pe*h2/4)*E_yf(Mxc+1,2*j)...
        -pc/(2*h2)*E_yf(Mxc+2,2*j)...
        +(pc*pd+1)*(H_z1f(Mxc+1,2*j) - H_z1c(Mxc,j))... 
        -pc/2*(pa*pb+1)*(K_zf(Mxc+1,2*j) - K_zc(Mxc,j))...
        +h2/2*(pe*pf+1)*(J_yf(Mxc+1,2*j)+2*J_yc(Mxc+1,j));    
 
       
       RHS4(Mxf+1)= E_yf(Mxc+1,2*j)+E_yf(Mxc+1,2*j-1)-2*E_yc(Mxc+1,j);

       
       RHS4(Mxf+2)=...
        -pc/(4*h2)*E_yc_next(Mxc,j) + (-2*epsilon*h2/tau + pc/(4*h2) + pe/2*h2)*E_yc(Mxc+1,j)...
        -pc/(4*h2)*E_yc(Mxc,j)...
        -(epsilon*h2/tau - pc/(2*h2) - pe*h2/4)*E_yf(Mxc+1,2*j-1)...
        -pc/(2*h2)*E_yf(Mxc+2,2*j-1)...
        +(pc*pd+1)*(H_z1f(Mxc+1,2*j-1) - H_z1c(Mxc,j))...
        -pc/2*(pa*pb+1)*(K_zf(Mxc+1,2*j-1) - K_zc(Mxc,j))...
        +h2/2*(pe*pf+1)*(J_yf(Mxc+1,2*j-1) + 2*J_yc(Mxc+1,j)); 
     
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Lower level fine mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       for i=Mxc+2:Mxc+Mxf-1
            RHS4(Mxf+1+i-Mxc)=...
            (1 -pe*tau/(4*epsilon))*E_yf(i,2*j-1)...
           +tau*pc/4/epsilon/h2^2*(E_yf(i+1,2*j-1)-2*E_yf(i,2*j-1)+E_yf(i-1,2*j-1))...
           -tau*(pc*pd+1)/2/epsilon/h2*(H_z1f(i,2*j-1)-H_z1f(i-1,2*j-1))...
           +tau*(pa*pb + 1)*pc/4/epsilon/h2*(K_zf(i,2*j-1)- K_zf(i-1,2*j-1))...
           -tau*(pe*pf+1)/2/epsilon*J_yf(i,2*j-1);
       end
            RHS4(2*Mxf+1)=...
            (1-pe*tau/(4*epsilon))*E_yf(Mxc+Mxf,2*j-1)...
           +tau*pc/4/epsilon/h2^2*(E_yf(Mxc+Mxf+1,2*j-1)-2*E_yf(Mxc+Mxf,2*j-1)+E_yf(Mxc+Mxf-1,2*j-1))...
           -tau*(pc*pd+1)/2/epsilon/h2*(H_z1f(Mxc+Mxf,2*j-1)-H_z1f(Mxc+Mxf-1,2*j-1))...
           +tau*(pa*pb + 1)*pc/4/epsilon/h2*(K_zf(Mxc+Mxf,2*j-1)- K_zf(Mxc+Mxf-1,2*j-1))...
           -tau*(pe*pf+1)/2/epsilon*J_yf(Mxc+Mxf,2*j-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Solve E_yf in the U structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        XXXX=chase(a11,a22,a33,RHS4);
         for i=Mxf+Mxc:-1:Mxc+1
              E_yf_next(i,2*j)=XXXX(Mxc+Mxf+1-i);
         end
         
        E_yc_next(Mxc+1,j)=XXXX(Mxf+1);
        %%%%% Solve J_yc on Mxc + 1
        J_yc_next(Mxc+1,j) = (pe*pf)*J_yc(Mxc+1,j)...
                           + pe/2*(E_yc_next(Mxc+1,j) + E_yc(Mxc+1,j));
        
        for i=Mxc+1:Mxc+Mxf   
              E_yf_next(i,2*j-1)=XXXX(i-Mxc+Mxf+1);
        end
            
 
        
         

        H_z2c(Mxc,j) = pc*pd*H_z1c(Mxc,j)...
               -pc/(4*h2)*( E_yc_next(Mxc+1,j) -E_yc_next(Mxc,j))...
               -pc/(4*h2)*( E_yc(Mxc+1,j) - E_yc(Mxc,j))...
               -pc*(pa*pb + 1)/2*K_zc(Mxc,j);
           
        K_zc_next(Mxc,j) = pa*pb*K_zc(Mxc,j)... 
                        + pa/2*(H_z2c(Mxc,j) + H_z1c(Mxc,j));

           
         for i=Mxc+1:Mxc+Mxf
           H_z2f(i,2*j-1) = pc*pd*H_z1f(i,2*j-1)...
           -pc/(2*h2)*( E_yf_next(i+1,2*j-1)-E_yf_next(i,2*j-1))...
           -pc/(2*h2)*( E_yf(i+1,2*j-1) - E_yf(i,2*j-1) )...
           -pc*(pa*pb+1)/2*K_zf(i,2*j-1); 
         end
 
        %%%%%%%%%%%%%%%% Solve for K_zf         
        for i=Mxc+1:Mxc+Mxf
            K_zf_next(i,2*j-1) = pa*pb*K_zf(i,2*j-1)...
                 + pa/2*(H_z2f(i,2*j-1) + H_z1f(i,2*j-1));
        end
        %%%%%%%%%%%%%%%% Solve for J_yf
        for i=Mxc+1:Mxc+Mxf   
              J_yf_next(i,2*j-1) = (pe*pf)*J_yf(i,2*j-1)...
                    + pe/2*(E_yf_next(i,2*j-1) + E_yf(i,2*j-1));  
        end
        
        for i=Mxc+1:Mxc+Mxf   
              J_yf_next(i,2*j) = (pe*pf)*J_yf_next(i,2*j)...
                    + pe/2*(E_yf_next(i,2*j) + E_yf(i,2*j));  
        end
        
       for i=Mxc+1:Mxc+Mxf
            H_z2f(i,2*j) = pc*pd*H_z1f(i,2*j)...
                -pc/(2*h2)*(E_yf_next(i+1,2*j)-E_yf_next(i,2*j))...
                -pc/(2*h2)*(E_yf(i+1,2*j) - E_yf(i,2*j))...
                -pc/2*(pa*pb+1)*K_zf(i,2*j);
       end
       
       %%%%%%%%%%%%%%%% Solve for K_zf         
        for i=Mxc+1:Mxc+Mxf
            K_zf_next(i,2*j) = pa*pb*K_zf(i,2*j)...
                 + pa/2*(H_z2f(i,2*j) + H_z1f(i,2*j));
        end 
    
      
      end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  solve in y direction     stage3     
      for i=1:Mxc+Mxf
           if i>=1&&i<=Mxc
          RHS5(1)= E_x1c_next(i,2)...
              +tau/4/epsilon/h2*(H_z2c(i,2)- H_z2c(i,1))...
              +tau^2/64/epsilon/mu/h2^2*(E_x1c_next(i,3)-2*E_x1c_next(i,2)+E_x1c_next(i,1))...
              -tau*(pex*pfx+1)/4/epsilon*J_x1c_next(i,2)...
              -tau*pex/16/epsilon*E_x1c_next(i,2);
           for j=3:Myc-1
               RHS5(j-1)= E_x1c_next(i,j)...
                   +tau/4/epsilon/h2*(H_z2c(i,j)- H_z2c(i,j-1))...
                   +tau^2/64/epsilon/mu/h2^2*(E_x1c_next(i,j+1)-2*E_x1c_next(i,j)+E_x1c_next(i,j-1))...
                   -tau*(pex*pfx+1)/4/epsilon*J_x1c_next(i,j)...
                   -tau*pex/16/epsilon*E_x1c_next(i,j);
           end
           RHS5(Myc-1)= E_x1c_next(i,Myc)+tau/4/epsilon/h2*(H_z2c(i,Myc)- H_z2c(i,Myc-1))...
               +tau^2/64/epsilon/mu/h2^2*(E_x1c_next(i,Myc+1)-2*E_x1c_next(i,Myc)+E_x1c_next(i,Myc-1))...
               -tau*(pex*pfx+1)/4/epsilon*J_x1c_next(i,Myc)...
               -tau*pex/16/epsilon*E_x1c_next(i,Myc);
           X5=chase(aa,bb,cc,RHS5);
           for j=2:Myc
            E_xc_next(i,j)=X5(j-1);
           end
         
           %%%%%% solve explicitly H_zc at time level $k+1$

           for j=1:Myc
               H_zc_next(i,j)=H_z2c(i,j)...
                   +tau/8/mu/h2*( E_xc_next(i,j+1)-E_xc_next(i,j)...
                   + E_x1c_next(i,j+1)-E_x1c_next(i,j));
           end

           %%%%%% solve explicitly J_xc at time level $k+1$
           %%%% Careful
           for j=2:Myc
                J_xc_next(i,j)=pex*pfx*J_x1c_next(i,j)...
                 + pex/4*(E_xc_next(i,j) + E_x1c_next(i,j));
           end
           %%%%%%%%%%%%%%%%%%%%%%%
           
           
           %%%%%%%
           elseif i>=Mxc+1&&i<=(Mxc+Mxf)
               RHS6(1)=E_x1f_next(i,2)+tau/2/epsilon/h2*(H_z2f(i,2)-H_z2f(i,1))...
                   +tau^2/16/epsilon/mu/h2^2*(E_x1f_next(i,3)-2*E_x1f_next(i,2)+E_x1f_next(i,1))...
                   -tau*(pex*pfx+1)/4/epsilon*J_x1f_next(i,2)...
                   -tau*pex/16/epsilon*E_x1f_next(i,2);
              
               for j=3:Myf-1
                  RHS6(j-1)=E_x1f_next(i,j)+tau/2/epsilon/h2*(H_z2f(i,j)-H_z2f(i,j-1))...
                      +tau^2/16/epsilon/mu/h2^2*(E_x1f_next(i,j+1)-2*E_x1f_next(i,j)+E_x1f_next(i,j-1))...
                      -tau*(pex*pfx+1)/4/epsilon*J_x1f_next(i,j)...
                      -tau*pex/16/epsilon*E_x1f_next(i,j);
               end
               RHS6(Myf-1)=E_x1f_next(i,Myf)+tau/2/epsilon/h2*(H_z2f(i,Myf)-H_z2f(i,Myf-1))...
                   +tau^2/16/epsilon/mu/h2^2*(E_x1f_next(i,Myf+1)-2*E_x1f_next(i,Myf)+E_x1f_next(i,Myf-1))...
                   -tau*(pex*pfx+1)/4/epsilon*J_x1f_next(i,Myf)...
                   -tau*pex/16/epsilon*E_x1f_next(i,Myf);

               X6=chase(aaa,bbb,ccc,RHS6);
              for j=2:Myf
               E_xf_next(i,j)=X6(j-1);
              end
%                 
           %%%%%% solve explicitly H_z1f at time level $k+1$
           for j=1:Myf
               H_zf_next(i,j)=H_z2f(i,j)+tau/4/mu/h2*(E_xf_next(i,j+1)...
                   -E_xf_next(i,j) + E_x1f_next(i,j+1)-E_x1f_next(i,j) );
           end
           %%%%%% solve explicitly J_xf at time level $k+1$
           for j=2:Myf
               J_xf_next(i,j) = pex*pfx*J_x1f_next(i,j)...
               +pex/4*(E_xf_next(i,j) + E_x1f_next(i,j));
           end
      
           
           end
      end
        
       for i=Mxc+1:Mxc+Mxf
            for j=1:Myc
                H_zf1_next(i,j)=(H_zf_next(i,2*j)+H_zf_next(i,2*j-1))/2; 
                K_zf1_next(i,j)=(K_zf_next(i,2*j)+K_zf_next(i,2*j-1))/2;             
            end
       end
       for i=Mxc+1:Mxc+Mxf+1
            for j=1:Myc
                E_yf1_next(i,j)=(E_yf_next(i,2*j)+E_yf_next(i,2*j-1))/2;
                J_yf1_next(i,j)=(J_yf_next(i,2*j)+J_yf_next(i,2*j-1))/2;            
            end
       end

  % Discrete Energy
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute partial energy
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j1 = 1:Mxc
           W_c1(j1) = W_c1(j1) + (1/4*tau*gamma_e*epsilon^-1*omega_pe^-2*...
             norm(J_xc_next(j1,1:Myc+1) + J_x1c_next(j1,1:Myc+1))^2 ...
             +1/4*tau*gamma_e*epsilon^-1*omega_pe^-2*...
             norm(J_x1c_next(j1,1:Myc+1) + J_xc(j1,1:Myc+1))^2 ...
             +1/2*tau*gamma_e*epsilon^-1*omega_pe^-2*...
             norm(J_yc_next(j1,1:Myc) + J_yc(j1,1:Myc))^2 ...
             + 1/2*tau*gamma_m*mu^-1*omega_pm^-2*...
             norm(K_zc_next(j1,1:Myc) + K_zc(j1,1:Myc))^2)*4*h2*h2;
  end
  
  for j2 = Mxc+2:Mxc+Mxf+1
  W_f1(j2) = W_f1(j2) + (1/2*tau*gamma_e*epsilon^-1*omega_pe^-2*...
             norm(J_yf_next(j2,1:Myf) + J_yf(j2,1:Myf))^2)*h2*h2;    
  end
  
  for j3 = Mxc+1:Mxc+Mxf
  W_f2(j3) = W_f2(j3) + (1/4*tau*gamma_e*epsilon^-1*omega_pe^-2*...
             norm(J_xf_next(j3,1:Myf+1) + J_x1f_next(j3,1:Myf+1))^2 ...
             + 1/4*tau*gamma_e*epsilon^-1*omega_pe^-2*...
             norm(J_x1f_next(j3,1:Myf+1) + J_xf(j3,1:Myf+1))^2 ...
             + 1/2*tau*gamma_m*mu^-1*omega_pm^-2*...
             norm(K_zf_next(j3,1:Myf) + K_zf(j3,1:Myf))^2)*h2*h2;
  end
  
  W_c3 = W_c3 + (1/2*tau*gamma_e*epsilon^-1*omega_pe^-2*...
            norm(J_yc_next(Mxc+1,1:Myc) + J_yc(Mxc+1,1:Myc))^2)*2*h2*h2;
  
  W_f3 = W_f3 + (1/2*tau*gamma_m*epsilon^-1*omega_pe^-2*...
            norm(J_yf_next(Mxc+1,1:Myf) + J_yf(Mxc+1,1:Myf))^2)*.5*h2*h2;      

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Swap for E_x,E_y,H_z
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  E_xf(Mxc+1:Mxc+Mxf+1,1:Myf+1) = E_xf_next(Mxc+1:Mxc+Mxf+1,1:Myf+1);
  E_xc(1:Mxc+1,1:Myc+1) = E_xc_next(1:Mxc+1,1:Myc+1);
  E_yf(Mxc+1:Mxc+Mxf+1,1:Myf) = E_yf_next(Mxc+1:Mxc+Mxf+1,1:Myf);
  E_yc(1:Mxc+1,1:Myc) = E_yc_next(1:Mxc+1,1:Myc);
  H_zf(Mxc+1:Mxc+Mxf,1:Myf) = H_zf_next(Mxc+1:Mxc+Mxf,1:Myf);
  H_zc(1:Mxc,1:Myc) = H_zc_next(1:Mxc,1:Myc);
  H_zf1(Mxc+1:Mxc+Mxf,1:Myc) = H_zf1_next(Mxc+1:Mxc+Mxf,1:Myc);
  E_yf1(Mxc+1:Mxc+Mxf,1:Myc) = E_yf1_next(Mxc+1:Mxc+Mxf,1:Myc);
  
  J_xf(Mxc+1:Mxc+Mxf+1,1:Myf+1) = J_xf_next(Mxc+1:Mxc+Mxf+1,1:Myf+1);
  J_xc(1:Mxc+1,1:Myc+1) = J_xc_next(1:Mxc+1,1:Myc+1);
  J_yf(Mxc+1:Mxc+Mxf+1,1:Myf) = J_yf_next(Mxc+1:Mxc+Mxf+1,1:Myf);
  J_yc(1:Mxc+1,1:Myc) = J_yc_next(1:Mxc+1,1:Myc);
  K_zf(Mxc+1:Mxc+Mxf,1:Myf) = K_zf_next(Mxc+1:Mxc+Mxf,1:Myf);
  K_zc(1:Mxc,1:Myc) = K_zc_next(1:Mxc,1:Myc);
  K_zf1(Mxc+1:Mxc+Mxf,1:Myc) = K_zf1_next(Mxc+1:Mxc+Mxf,1:Myc);
  J_yf1(Mxc+1:Mxc+Mxf,1:Myc) = J_yf1_next(Mxc+1:Mxc+Mxf,1:Myc);
  
  J_x1c(1:Mxc+1,1:Myc+1) = J_x1c_next(1:Mxc+1,1:Myc+1);
  J_x1f(Mxc+1:Mxc+Mxf+1,1:Myf+1) = J_x1f_next(Mxc+1:Mxc+Mxf+1,1:Myf+1);
  end % end of k loop
  
        for j1=1:Mxc
      W_c1(j1) = W_c1(j1) + 4*h2*h2*epsilon*norm(E_yc(j1,1:Myc))^2 ... 
      +4*h2*h2*epsilon*norm(E_xc(j1,1:Myc+1))^2 ...
      +4*h2*h2*mu*norm(H_zc(j1,1:Myc))^2 ...
      +4*h2*h2*mu^-1*omega_pm^(-2)*norm(K_zc(j1,1:Myc))^2 ...
      +4*h2*h2*epsilon^(-1)*omega_pe^(-2)*norm(J_xc(j1,1:Myc+1))^2 ...
      +4*h2*h2*epsilon^(-1)*omega_pe^(-2)*norm(J_yc(j1,1:Myc))^2;          
      end
      
      for j2=Mxc+2:Mxc+Mxf+1
            W_f1(j2) = W_f1(j2) + h2*h2*epsilon*norm(E_yf(j2,1:Myf))^2 ...
                      +h2*h2*epsilon^-1*omega_pe^-2*...
                      norm(J_yf(j2,1:Myf))^2;
      end
      
      for j3=Mxc+1:Mxc+Mxf
          W_f2(j3) = W_f2(j3) + h2*h2*epsilon*norm(E_xf(j3,1:Myf+1))^2 ...
          +h2*h2*mu*norm(H_zf(j3,1:Myf))^2 ...
          +h2*h2*epsilon^-1*omega_pe^-2*norm(J_xf(j3,1:Myf+1))^2 ...
          +h2*h2*mu^-1*omega_pm^-2*norm(K_zf(j3,1:Myf))^2;
      end
      

        W_c = (sum(W_c1))+(2*h2*h2*epsilon*norm(E_yc(Mxc+1,1:Myc))^2)...
                          +(2*h2*h2*epsilon^-1*omega_pe^-2*norm(J_yc(Mxc+1,1:Myc))^2);             
                          
        W_f = (sum(W_f1)+sum(W_f2))+(0.5*h2*h2*epsilon*norm(E_yf(Mxc+1,1:Myf))^2)...
                                    +(0.5*h2*h2*epsilon^-1*omega_pe^-2*norm(J_yf(Mxc+1,1:Myf))^2);
        W_c = W_c + W_c3;
        W_f = W_f + W_f3; 
        W = W_c + W_f;

fprintf("The relative energy error is %d \n",abs(W - .25)/.25);
sprintf('\n');
endtime = toc;
fprintf("The CPU time is %d",endtime); 
sprintf('\n'); 

   
  
