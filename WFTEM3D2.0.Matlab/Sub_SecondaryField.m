%+ Iteration and output results.
%
function [t_iteration_H,DBZ_Rx] = Sub_SecondaryField(XI,YJ,ZK,dx,dy,dz,...
    Rx,Tx,iter_n_sec,alpha_sec,model_EC,EX,EY,EZ,HX,HY,HZ)
%
%--------------------------------------------------------------------------
% In: XI,YJ,ZK,dx,dy,dz,Rx,Tx,iter_n_sec,alpha_sec,model_EC, 
%     EX,EY,EZ,HX,HY,HZ   
% Out: t_iteration_H,DBZ_Rx
% Description:
% Calculatie secondary field and output results.
%
% Method:
% Step Maxwell¡¯s equations in time using a staggered grid and
% a modified DuFort-Frankel method. See 3D Finite-difference Transient
% Electromagnetic Modeling with a Whole-space Initial Field for detail.
%
% Current Code Owner: <Fei Li and Jiulong Cheng>
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. Fei Li
%
% Declarations:
mu_0 = 4*pi*10^(-7); % Magnetic permeability of vacuum.
% XI, YJ, ZK      % Number of cells in the x-, y-, and z-directions.
% dx, dy, dz      % Grid size in the x-, y-, and z-directions.
% Tx              % Tx loop diagonal points / ground wire endpoints.
% Rx              % Rx points: start/end/interval.
% iter_n_sec     % Iteration number of secondary field.
% alpha_sec      % Coefficient of time step size for secondary field.
% model_EC        % Model conductivity.
% EX(:,:,:,:)     % Ex.
% EY(:,:,:,:)     % Ey.
% EZ(:,:,:,:)     % Ez.
% HX(:,:,:,:)     % Hx.
% HY(:,:,:,:)     % Hy.
% HZ(:,:,:,:)     % Hz.
% DBZ_Rx          % dBz/dt at receivers.
t_iteration_E = zeros(iter_n_sec,1);% Time of electric field.
t_iteration_H = zeros(iter_n_sec,1);% Time of magnetic field.
EC_x = zeros(XI,YJ,ZK);             % Conductivity of nodes.
EC_y = zeros(XI,YJ,ZK);             % ditto.
EC_z = zeros(XI,YJ,ZK);             % ditto.
dt = zeros(2);                      % Time step size.
% i,j,k,iter_n                      % Temporary loop control variables.
% Tx_z                              % Transmitter depth.
% gamma                             % Artificial dielectric constant.
% model_EC_min                      % Minimum conductivity in model.
% Rx_start_x, Rx_end_x, Rx_step_x   % Rx along x-axis: start/end/interval.
% Rx_start_y, Rx_end_y, Rx_step_y   % Rx along y-axis: start/end/interval.
% Rx_start_z, Rx_end_z, Rx_step_z   % Rx along z-axis: start/end/interval.
% R                                 % Position of each receiver.
%- End of header ----------------------------------------------------------

fprintf('Calculating secondary field...\n')
%--------------------------------------------------------------------------
% [1.0] Model and configuration initialization:
%--------------------------------------------------------------------------
Tx_z=Tx(3);                     % Transmitter depth
model_EC_min=min(model_EC(:));  % Minimum conductivity in model
% Initial time of electric field:
t_iteration_E(1) = 1.13*mu_0*model_EC_min*dz(Tx_z)^2;
% Calculate the time step size using equation 28 in paper:
dt(1)=alpha_sec*dz(Tx_z)*sqrt(mu_0*model_EC_min*t_iteration_E(1)/6);
% Initial time of magnetic field:
t_iteration_H(1) = t_iteration_E(1)+dt(1)/2;      
% Receiver positions:
Rx_start_x=Rx(1);Rx_end_x=Rx(2);Rx_step_x=Rx(3);
Rx_start_y=Rx(4);Rx_end_y=Rx(5);Rx_step_y=Rx(6);
Rx_start_z=Rx(7);Rx_end_z=Rx(8);Rx_step_z=Rx(9);
iter_n = 1;
for k=Rx_start_z:Rx_step_z:Rx_end_z
    for j=Rx_start_y:Rx_step_y:Rx_end_y
        for i=Rx_start_x:Rx_step_x:Rx_end_x
            R(iter_n,:) = [i,j,k]; % Position of each receiver.
            iter_n = iter_n+1;
        end
    end
end
% Calculate dBz/dt at receivers using equation 21 in paper:
for i=1:1:length(R(:,1))
    DBZ_Rx(1,i) = ((EX(R(i,1),R(i,2)+1,R(i,3)+1,1) ...
        -EX(R(i,1),R(i,2),R(i,3)+1,1))/dy(R(i,2))  ...
        -(EY(R(i,1)+1,R(i,2),R(i,3)+1,1)           ...
        -EY(R(i,1),R(i,2),R(i,3)+1,1))/dx(R(i,1)));
end
%--------------------------------------------------------------------------
% [2.0] Calculate conductivity of grid nodes:
%--------------------------------------------------------------------------
% Equation 13 in paper:
for k = 2:ZK % 1st & (ZK+1)th are on the boundary face, Dirichlet BC (u=0).
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
% [3.0] Calculate secondary field:
%--------------------------------------------------------------------------
iter_n = 2;
while iter_n <= iter_n_sec
    % Calculate the electric field time of the iter_nth iteration:
    t_iteration_E(iter_n) = t_iteration_E(iter_n-1)+dt(1);
    % Calculate the time step size using equation 23 in paper:
    dt(2) = alpha_sec*dz(Tx_z)*sqrt(mu_0 ...
        *model_EC_min*t_iteration_E(iter_n)/6);
    % Calculate the magnetic field time of the iter_nth iteration:
    t_iteration_H(iter_n) = t_iteration_H(iter_n-1)+(dt(1)+dt(2))/2;
    % Calculate artificial dielectric constant using equation 22 in paper:
    gamma = (4/mu_0)*(dt(1)/dz(Tx_z))^2;
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
    EX(1:XI,1:YJ+1,1,2) = 0;
    EX(1:XI,1:YJ+1,ZK+1,2) = 0;
    EX(1:XI,1,1:ZK+1,2) = 0;
    EX(1:XI,YJ+1,1:ZK+1,2) = 0;
    EY(1:XI+1,1:YJ,1,2) = 0 ;
    EY(1:XI+1,1:YJ,ZK+1,2) = 0 ;
    EY(1,1:YJ,1:ZK+1,2) = 0;
    EY(XI+1,1:YJ,1:ZK+1,2) = 0;
    EZ(1:XI+1,1,1:ZK,2) = 0;
    EZ(1:XI+1,YJ+1,1:ZK,2) = 0;
    EZ(1,1:YJ+1,1:ZK,2) = 0;
    EZ(XI+1,1:YJ+1,1:ZK,2) = 0;
    % Update the value of Hx using equation 18 in paper:
    for k = 1:ZK
        for j = 1:YJ
            for i = 1:XI+1
                HX(i,j,k,2) = HX(i,j,k,1)                 ...
                    +(dt(1)+dt(2))/(2*mu_0)       ...
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
                    +(dt(1)+dt(2))/(2*mu_0)       ...
                    *((EZ(i+1,j,k,2)-EZ(i,j,k,2))         ...
                    /dx(i)-(EX(i,j,k+1,2)-EX(i,j,k,2))/dz(k));
            end
        end
    end
    % Update the value of Hz using equation 20 in paper, start from
    % the bottom boundary and the top boundary of the model where Hz=0:
    HZ(1:XI,1:YJ,1,1) = 0;
    HZ(1:XI,1:YJ,ZK+1,1) = 0;
    for k = ZK+1:-1:Tx_z+2
        for j = 1:YJ
            for i = 1:XI
                HZ(i,j,k-1,2) = HZ(i,j,k,2) ...
                    +dz(k-1)/dx(i)*(HX(i+1,j,k-1,2)-HX(i,j,k-1,2)) ...
                    +dz(k-1)/dy(j)*(HY(i,j+1,k-1,2)-HY(i,j,k-1,2));
            end
        end
    end
    for k = 2:Tx_z
        for j = 1:YJ
            for i = 1:XI
                HZ(i,j,k,2) = HZ(i,j,k-1,2) ...
                    -dz(k-1)/dx(i)*(HX(i+1,j,k-1,2)-HX(i,j,k-1,2)) ...
                    -dz(k-1)/dy(j)*(HY(i,j+1,k-1,2)-HY(i,j,k-1,2));
            end
        end
    end
    % Calculate dBz/dt at receivers using equation 21 in paper:  
    for i=1:1:length(R(:,1))
        DBZ_Rx(iter_n,i) = ((EX(R(i,1),R(i,2)+1,R(i,3)+1,1) ...
            -EX(R(i,1),R(i,2),R(i,3)+1,1))/dy(R(i,2))  ...
            -(EY(R(i,1)+1,R(i,2),R(i,3)+1,1)           ...
            -EY(R(i,1),R(i,2),R(i,3)+1,1))/dx(R(i,1)));
    end
    % The electromagnetic field are used for the next iteration:
    EX(:,:,:,1) = EX(:,:,:,2);
    EY(:,:,:,1) = EY(:,:,:,2);
    EZ(:,:,:,1) = EZ(:,:,:,2);
    HX(:,:,:,1) = HX(:,:,:,2);
    HY(:,:,:,1) = HY(:,:,:,2);
    HZ(:,:,:,1) = HZ(:,:,:,2);
    dt(1) = dt(2);
    iter_n = iter_n+1;
    % Print progress message every 500 iterations:
    if mod(iter_n, 500) == 0
        fprintf('%6d (%.1f%%)\n', iter_n, 100.0*iter_n/iter_n_sec);
    end
end
fprintf('Secondary field calculation finished.\n')
%--------------------------------------------------------------------------
% [4.0] Output results (time and dBz/dt at receivers):
%--------------------------------------------------------------------------
% Output time:
save('.\Data\Result_time.txt','t_iteration_H','-ascii')
% Output dBz/dt:
save('.\Data\Result_dBz.txt','DBZ_Rx','-ascii')
