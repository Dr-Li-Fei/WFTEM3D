%+ Iteration and output results.
%
function [EX,EY,EZ,HX,HY,HZ] = Sub_PrimaryField(XI,YJ,ZK,dx,dy,dz, ...
    Tx,iter_n_pri,alpha_pri)
%
%--------------------------------------------------------------------------
% In: XI,YJ,ZK,dx,dy,dz,Tx,iter_n_pri,alpha_pri
% Out: EX,EY,EZ,HX,HY,HZ
% Description:
% Calculate primary field.
%
% Method:
% Wire is modeled as volume current.First calculate the primary field using
% a whole-space homogeneous model, then calculate the secondary field using 
% the true model. (Relevant paper will be provided in subsequent updates.)
%
% Current Code Owner: <Fei Li and Jiulong Cheng>
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. Fei Li
% 2.0      01/03/25  Magnetic dipole-based -> volume current-based. Fei Li
%
% Declarations:
% XI, YJ, ZK      % Number of cells in the x-, y-, and z-directions.
% dx, dy, dz      % Grid size in the x-, y-, and z-directions.
% Tx              % Tx loop diagonal points / ground wire endpoints.
% iter_n_pri      % Iteration number of primary field.
% alpha_pri       % Coefficient of time step size for primary field
EX = zeros(XI,YJ+1,ZK+1,2);  % Ex.
EY = zeros(XI+1,YJ,ZK+1,2);  % Ey.
EZ = zeros(XI+1,YJ+1,ZK,2);  % Ez.
HX = zeros(XI+1,YJ,ZK,2);    % Hx
HY = zeros(XI,YJ+1,ZK,2);    % Hy.
HZ = zeros(XI,YJ,ZK+1,2);    % Hz.
mu_0 = 4*pi*10^(-7);         % Magnetic permeability of vacuum.
t_iteration_E = zeros(iter_n_pri,1);% Time of electric field.
t_iteration_H = zeros(iter_n_pri,1);% Time of magnetic field.
dt = zeros(2);                      % Time step size.
% i,j,k,iter_n                      % Temporary loop control variables.
% Tx_x1,Tx_y1,Tx_x2,Tx_y2,Tx_z      % Position of Tx.
% EC                                % Conductivity of homogeneous model.
% current                           % Current.
% gamma                             % Artificial dielectric constant.
%- End of header ----------------------------------------------------------

fprintf('Calculating primary field...\n')
%--------------------------------------------------------------------------
% [1.0] Model and configuration initialization:
%--------------------------------------------------------------------------
EC = 10;                            % Conductivity of homogeneous model.
current=1;                          % Current is 1 A by default
Tx_x1=Tx(1);Tx_y1=Tx(2);Tx_z=Tx(3); % Position of Tx.
Tx_x2=Tx(4);Tx_y2=Tx(5);            % Position of Tx.
t_iteration_E(1) = 0;               % Initial time of electric field.
% Calculate the time step size using equation 28 in paper:
dt(1) =2.26*mu_0*EC*dz(Tx_z)^2;
%Initial time of magnetic field:
t_iteration_H(1) = t_iteration_E(1)+dt(1)/2;
%--------------------------------------------------------------------------
% [2.0] Calculate primary field:
%--------------------------------------------------------------------------
iter_n = 2;
while iter_n <= iter_n_pri
    % Calculate the electric field time of the iter_nth iteration:
    t_iteration_E(iter_n) = t_iteration_E(iter_n-1)+dt(1);
    % Calculate the time step size using equation 23 in paper:
    dt(2) = alpha_pri*dz(Tx_z)*sqrt(mu_0 ...
        *EC*t_iteration_E(iter_n)/6);
    % Calculate the magnetic field time of the iter_nth iteration:
    t_iteration_H(iter_n) = t_iteration_H(iter_n-1)+(dt(1)+dt(2))/2;
    % Calculate artificial dielectric constant using equation 22 in paper:
    gamma = (4/mu_0)*(dt(1)/dz(Tx_z))^2;
    % Update the value of Ex using equation 12 in paper:
    for k = 2:ZK
        for j = 2:YJ
            for i = 1:XI
                EX(i,j,k,2) = (2*gamma-EC*dt(1))         ...
                    /(2*gamma+EC*dt(1))*EX(i,j,k,1)      ...
                    +4*dt(1)/(2*gamma+EC*dt(1))          ...
                    *((HZ(i,j,k,1)-HZ(i,j-1,k,1))/(dy(j-1)+dy(j)) ...
                    -(HY(i,j,k,1)-HY(i,j,k-1,1))/(dz(k-1)+dz(k)));
            end
        end
    end
    % Wire along x-direction if Tx_x1~=Tx_x2:
    if Tx_x1~=Tx_x2 
        for i = Tx_x1:Tx_x2
            j = Tx_y1;
            k = Tx_z+1;
            EX(i,j,k,2) = (2*gamma-EC*dt(1))         ...
                /(2*gamma+EC*dt(1))*EX(i,j,k,1)      ...
                +2*dt(1)/(2*gamma+EC*dt(1))          ...
                *((HZ(i,j,k,1)-HZ(i,j-1,k,1))/((dy(j-1)+dy(j))/2) ...
                -(HY(i,j,k,1)-HY(i,j,k-1,1))/((dz(k-1)+dz(k))/2)  ...
                -current/((dy(j-1)+dy(j))*(dz(k-1)+dz(k))/4));
        end
    end
    % Anothor wire along x-direction if Tx_x1~=Tx_x2 and Tx_y1~=Tx_y2:
    if Tx_x1~=Tx_x2 && Tx_y1~=Tx_y2 
        for i = Tx_x1:Tx_x2
            j = Tx_y2+1;
            k = Tx_z+1;
            EX(i,j,k,2) = (2*gamma-EC*dt(1))         ...
                /(2*gamma+EC*dt(1))*EX(i,j,k,1)      ...
                +2*dt(1)/(2*gamma+EC*dt(1))          ...
                *((HZ(i,j,k,1)-HZ(i,j-1,k,1))/((dy(j-1)+dy(j))/2) ...
                -(HY(i,j,k,1)-HY(i,j,k-1,1))/((dz(k-1)+dz(k))/2)  ...
                -current*(-1)/((dy(j-1)+dy(j))*(dz(k-1)+dz(k))/4));
        end
    end
    % Update the value of Ey using equation 14 in paper:
    for k = 2:ZK
        for j = 1:YJ
            for i = 2:XI
                EY(i,j,k,2) = (2*gamma-EC*dt(1))         ...
                    /(2*gamma+EC*dt(1))*EY(i,j,k,1)      ...
                    +4*dt(1)/(2*gamma+EC*dt(1))          ...
                    *((HX(i,j,k,1)-HX(i,j,k-1,1))/(dz(k-1)+dz(k)) ...
                    -(HZ(i,j,k,1)-HZ(i-1,j,k,1))/(dx(i-1)+dx(i)));
            end
        end
    end
    % Wire along y-direction if Tx_y1~=Tx_y2:
    if Tx_y1~=Tx_y2
        for j = Tx_y1:Tx_y2
            i = Tx_x1;
            k = Tx_z+1;
            EY(i,j,k,2) = (2*gamma-EC*dt(1))         ...
                /(2*gamma+EC*dt(1))*EY(i,j,k,1)      ...
                +2*dt(1)/(2*gamma+EC*dt(1))          ...
                *((HX(i,j,k,1)-HX(i,j,k-1,1))/((dz(k-1)+dz(k))/2) ...
                -(HZ(i,j,k,1)-HZ(i-1,j,k,1))/((dx(i-1)+dx(i))/2)  ...
                -current*(-1)/((dz(k-1)+dz(k))*(dx(i-1)+dx(i))/4));
        end
    end
    % Anothor wire along y-direction if Tx_x1~=Tx_x2 and Tx_y1~=Tx_y2:
    if Tx_x1~=Tx_x2 && Tx_y1~=Tx_y2
        for j = Tx_y1:Tx_y2
            i = Tx_x2+1;
            k = Tx_z+1;
            EY(i,j,k,2) = (2*gamma-EC*dt(1))         ...
                /(2*gamma+EC*dt(1))*EY(i,j,k,1)      ...
                +2*dt(1)/(2*gamma+EC*dt(1))          ...
                *((HX(i,j,k,1)-HX(i,j,k-1,1))/((dz(k-1)+dz(k))/2) ...
                -(HZ(i,j,k,1)-HZ(i-1,j,k,1))/((dx(i-1)+dx(i))/2)  ...
                -current/((dz(k-1)+dz(k))*(dx(i-1)+dx(i))/4));
        end
    end
    % Update the value of Ez using equation 16 in paper:
    for k = 1:ZK
        for j = 2:YJ
            for i = 2:XI
                EZ(i,j,k,2) = (2*gamma-EC*dt(1))         ...
                    /(2*gamma+EC*dt(1))*EZ(i,j,k,1)      ...
                    +4*dt(1)/(2*gamma+EC*dt(1))          ...
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
        fprintf('%6d (%.1f%%)\n', iter_n, 100.0*iter_n/iter_n_pri);
    end
end
fprintf('Primary field calculation finished.\n')