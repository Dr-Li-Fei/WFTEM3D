%+ Iteration and output results.
%
function [t_iteration_H,DBZ_Rx] = Sub_Iteration(i0,j0,k0,XI,YJ,ZK,dx,dy,...
    dz,Rx_st,Rx_end,iter_n_max,Alpha,model_EC,EX,EY,EZ,HX,HY,HZ,t1_E,t1_H)
%
%--------------------------------------------------------------------------
% In: i0,j0,k0,XI,YJ,ZK,dx,dy,dz,Rx_st,Rx_end,iter_n_max,Alpha,model_EC,
%     EX,EY,EZ,HX,HY,HZ,t0_E,t0_H
% Out: t_iteration_H,DBZ_Rx
% Description:
% Iteration and output results.
%
% Method:
% Step Maxwell¡¯s equations in time using a staggered grid and
% a modified DuFort-Frankel method. See 3D Finite-difference Transient
% Electromagnetic Modeling with Whole-space Initial Field for detail.
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. <Fei Li and Jiulong Cheng>
%
% Declarations:
Permeability_vac = 4*pi*10^(-7); % Magnetic permeability of vacuum.
% i0, j0, k0      % Tx loop center is at the bottom face center
% of the cell (i0,j0,k0).
% XI, YJ, ZK      % Number of cells in the x-, y-, and z-directions.
% dx, dy, dz      % Grid size in the x-, y-, and z-directions.
% Rx_st, Rx_end   % Receivers are at the bottom face center of cells
% from (i0,Rx_st,k0) to (i0,Rx_end,k0).
% iter_n_max      % Maximum number of iterations.
% Alpha           % Coefficient in the time step size calculation equation.
% model_EC        % Model conductivity.
% EX(:,:,:,:)     % Ex.
% EY(:,:,:,:)     % Ey.
% EZ(:,:,:,:)     % Ez.
% HX(:,:,:,:)     % Hx.
% HY(:,:,:,:)     % Hy.
% HZ(:,:,:,:)     % Hz.
% t0_E            % Initial time of electric field.
% t0_H            % Initial time of magnetic field.
DBZ_Rx = zeros(iter_n_max,Rx_end-Rx_st+1);% dBz/dt at receivers.
t_iteration_E = zeros(iter_n_max,1);% Time of electric field.
t_iteration_H = zeros(iter_n_max,1);% Time of magnetic field.
EC_x = zeros(XI,YJ,ZK);             % Conductivity of nodes.
EC_y = zeros(XI,YJ,ZK);             % ditto.
EC_z = zeros(XI,YJ,ZK);             % ditto.
dt = zeros(2);                      % Time step size.
%integer::i,j,k,p                   % Temporary loop variables.
%integer::iter_n                    % Temporary iteration variable.
%character(len = 256)::n_Rx         % Number of receivers.
%real(8)::gamma                     % Artificial dielectric constant.
%- End of header ----------------------------------------------------------

%--------------------------------------------------------------------------
% [1.0] Calculate conductivity of grid nodes:
%--------------------------------------------------------------------------
% Equation 13 in paper:
for k = 2:ZK
    for j = 2:YJ
        for i = 1:XI
            EC_x(i,j,k) = (dy(j-1)*dz(k-1)*model_EC(i,j-1,k-1) ...
                +dy(j)*dz(k-1)*model_EC(i,j,k-1)               ...
                +dy(j-1)*dz(k)*model_EC(i,j-1,k)               ...
                +dy(j)*dz(k)*model_EC(i,j,k))                  ...
                /((dy(j-1)+dy(j))*(dz(k-1)+dz(k)));
        end
    end
end
% Equation 15 in paper:
for k = 2:ZK
    for j = 1:YJ
        for i = 2:XI
            EC_y(i,j,k) = (dx(i-1)*dz(k-1)*model_EC(i-1,j,k-1) ...
                +dx(i)*dz(k-1)*model_EC(i,j,k-1)               ...
                +dx(i-1)*dz(k)*model_EC(i-1,j,k)               ...
                +dx(i)*dz(k)*model_EC(i,j,k))                  ...
                /((dx(i-1)+dx(i))*(dz(k-1)+dz(k)));
        end
    end
end
% Equation 17 in paper:
for k = 1:ZK
    for j = 2:YJ
        for i = 2:XI
            EC_z(i,j,k) = (dx(i-1)*dy(j-1)*model_EC(i-1,j-1,k) ...
                +dx(i)*dy(j-1)*model_EC(i,j-1,k)               ...
                +dx(i-1)*dy(j)*model_EC(i-1,j,k)               ...
                +dx(i)*dy(j)*model_EC(i,j,k))                  ...
                /((dx(i-1)+dx(i))*(dy(j-1)+dy(j)));
        end
    end
end
%--------------------------------------------------------------------------
% [2.0] Get initial time:
%--------------------------------------------------------------------------
t_iteration_E(1) = t1_E;%Initial time of electric field.
t_iteration_H(1) = t1_H; %Initial time of magnetic field.
% Calculate the time step size using equation 28 in paper:
dt(1) = Alpha*dz(k0)*sqrt(Permeability_vac*model_EC(i0,j0,k0)*t1_E/6);
% Calculate dBz/dt at receivers at initial time using equation 21 in paper:
for p = Rx_st:1:Rx_end
    DBZ_Rx(1,p-Rx_st+1) = ((EX(i0,p+1,k0+1,1)-EX(i0,p,k0+1,1))/dy(p) ...
        -(EY(i0+1,p,k0+1,1)-EY(i0,p,k0+1,1))/dx(i0));
end
%--------------------------------------------------------------------------
% [3.0] Iteration:
%--------------------------------------------------------------------------
iter_n = 2
while iter_n <= iter_n_max
    % Calculate the electric field time of the iter_nth iteration:
    t_iteration_E(iter_n) = t_iteration_E(iter_n-1)+dt(1);
    % Calculate the time step size using equation 23 in paper:
    dt(2) = Alpha*dz(k0)*sqrt(Permeability_vac ...
        *model_EC(i0,j0,k0)*t_iteration_E(iter_n)/6);
    % Calculate the magnetic field time of the iter_nth iteration:
    t_iteration_H(iter_n) = t_iteration_H(iter_n-1)+(dt(1)+dt(2))/2;
    % Calculate artificial dielectric constant using equation 22 in paper:
    gamma = (4/Permeability_vac)*(dt(1)/dz(k0))^2;
    % Update the value of Ex using equation 12 in paper:
    for k = 2:ZK
        for j = 2:YJ
            for i = 1:XI
                EX(i,j,k,2) = (2*gamma-EC_x(i,j,k)*dt(1))         ...
                    /(2*gamma+EC_x(i,j,k)*dt(1))*EX(i,j,k,1)      ...
                    +4*dt(1)/(2*gamma+EC_x(i,j,k)*dt(1))          ...
                    *((HZ(i,j,k,1)-HZ(i,j-1,k,1))/(dy(j-1)+dy(j)) ...
                    -(HY(i,j,k,1)-HY(i,j,k-1,1))/(dz(k-1)+dz(k)));
            end
        end
    end
    % Update the value of Ey using equation 14 in paper:
    for k = 2:ZK
        for j = 1:YJ
            for i = 2:XI
                EY(i,j,k,2) = (2*gamma-EC_y(i,j,k)*dt(1))         ...
                    /(2*gamma+EC_y(i,j,k)*dt(1))*EY(i,j,k,1)      ...
                    +4*dt(1)/(2*gamma+EC_y(i,j,k)*dt(1))          ...
                    *((HX(i,j,k,1)-HX(i,j,k-1,1))/(dz(k-1)+dz(k)) ...
                    -(HZ(i,j,k,1)-HZ(i-1,j,k,1))/(dx(i-1)+dx(i)));
            end
        end
    end
    % Update the value of Ez using equation 16 in paper:
    for k = 1:ZK
        for j = 2:YJ
            for i = 2:XI
                EZ(i,j,k,2) = (2*gamma-EC_z(i,j,k)*dt(1))         ...
                    /(2*gamma+EC_z(i,j,k)*dt(1))*EZ(i,j,k,1)      ...
                    +4*dt(1)/(2*gamma+EC_z(i,j,k)*dt(1))          ...
                    *((HY(i,j,k,1)-HY(i-1,j,k,1))/(dx(i-1)+dx(i)) ...
                    -(HX(i,j,k,1)-HX(i,j-1,k,1))/(dy(j-1)+dy(j)));
            end
        end
    end
    % Dirichlet boundary condition:
    for i = 1:XI
        for j = 1:YJ+1
            EX(i,j,1,2) = 0;
            EX(i,j,ZK+1,2) = 0;
        end
    end
    for i = 1:XI
        for k = 1:ZK+1
            EX(i,1,k,2) = 0;
            EX(i,YJ+1,k,2) = 0;
        end
    end
    for i = 1:XI+1
        for j = 1:YJ
            EY(i,j,1,2) = 0 ;
            EY(i,j,ZK+1,2) = 0 ;
        end
    end
    for k = 1:ZK+1
        for j = 1:YJ
            EY(1,j,k,2) = 0;
            EY(XI+1,j,k,2) = 0;
        end
    end
    for i = 1:XI+1
        for k = 1:ZK
            EZ(i,1,k,2) = 0;
            EZ(i,YJ+1,k,2) = 0;
        end
    end
    for k = 1:ZK
        for j = 1:YJ+1
            EZ(1,j,k,2) = 0;
            EZ(XI+1,j,k,2) = 0;
        end
    end
    % Update the value of Hx using equation 18 in paper:
    for k = 1:ZK
        for j = 1:YJ
            for i = 1:XI+1
                HX(i,j,k,2) = HX(i,j,k,1)                 ...
                    +(dt(1)+dt(2))/(2*Permeability_vac)       ...
                    *((EY(i,j,k+1,2)-EY(i,j,k,2))         ...
                    /dz(k)-(EZ(i,j+1,k,2)-EZ(i,j,k,2))/dy(j));
            end
        end
    end
    % Update the value of Hy using equation 19 in paper:
    for k = 1:ZK
        for j = 1:YJ+1
            for i = 1:XI
                HY(i,j,k,2) = HY(i,j,k,1)                 ...
                    +(dt(1)+dt(2))/(2*Permeability_vac)       ...
                    *((EZ(i+1,j,k,2)-EZ(i,j,k,2))         ...
                    /dx(i)-(EX(i,j,k+1,2)-EX(i,j,k,2))/dz(k));
            end
        end
    end
    % Update the value of Hz using equation 20 in paper, start from
    % the bottom boundary and the top boundary of the model where Hz=0:
    HZ(1:XI,1:YJ,1,1) = 0;
    HZ(1:XI,1:YJ,ZK+1,1) = 0;
    for k = ZK+1:-1:k0+2
        for j = 1:YJ
            for i = 1:XI
                HZ(i,j,k-1,2) = HZ(i,j,k,2) ...
                    +dz(k-1)/dx(i)*(HX(i+1,j,k-1,2)-HX(i,j,k-1,2)) ...
                    +dz(k-1)/dy(j)*(HY(i,j+1,k-1,2)-HY(i,j,k-1,2));
            end
        end
    end
    for k = 2:k0
        for j = 1:YJ
            for i = 1:XI
                HZ(i,j,k,2) = HZ(i,j,k-1,2) ...
                    -dz(k-1)/dx(i)*(HX(i+1,j,k-1,2)-HX(i,j,k-1,2)) ...
                    -dz(k-1)/dy(j)*(HY(i,j+1,k-1,2)-HY(i,j,k-1,2));
            end
        end
    end
    % Calculate dBz/dt at receivers using equation 21 in paper:
    for p = Rx_st:1:Rx_end
        DBZ_Rx(iter_n,p-Rx_st+1) = ((EX(i0,p+1,k0+1,2)-EX(i0,p,k0+1,2))/dy(p) ...
            -(EY(i0+1,p,k0+1,2)-EY(i0,p,k0+1,2))/dx(i0));
    end
    % The electromagnetic field are used for the next iteration:
    EX(:,:,:,1) = EX(:,:,:,2);
    EY(:,:,:,1) = EY(:,:,:,2);
    EZ(:,:,:,1) = EZ(:,:,:,2);
    HX(:,:,:,1) = HX(:,:,:,2);
    HY(:,:,:,1) = HY(:,:,:,2);
    HZ(:,:,:,1) = HZ(:,:,:,2);
    dt(1) = dt(2);
    iter_n = iter_n+1
end
fprintf('Computation finished.\n')
%--------------------------------------------------------------------------
% [4.0] Output results (time and dBz/dt at receivers):
%--------------------------------------------------------------------------
% Output time:
save('.\Data\Result_time.txt','t_iteration_H','-ascii')
% Output dBz/dt:
save('.\Data\Result_dBz.txt','DBZ_Rx','-ascii')

