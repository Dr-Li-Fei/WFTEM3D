%+ Calculate the subloop initial field.
%
function [EX_subloop,EY_subloop,EZ_subloop,HX_subloop,HY_subloop, ...
    HZ_subloop,t1_E,t1_H] = Sub_InitialField_SubLoop(i0,j0,k0,XI,YJ,ZK, ...
    model_EC,dx,dy,dz,subloop_x,subloop_y,L_subloop,Alpha)
%
% In: i0,j0,k0,XI,YJ,ZK,model_EC,dx,dy,dz,subloop_x,subloop_y,L_subloop,
%     Alpha (NOTE: model_EC is not necessary)
% Out: EX_subloop,EY_subloop,EZ_subloop,HX_subloop,HY_subloop,HZ_subloop,
%      t1_E,t1_H
% Description:
% Calculate the subloop initial field.
%
% Method:
% Calculate the whole-space initial field excited by magnetic dipole
% sources at an initial time after the current is switched off.
% See 3D Finite-difference Transient Electromagnetic Modeling with
% Whole-space Initial Field for detail.
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. <Fei Li and Jiulong Cheng>
%
% Declarations:
Permeability_vac = 4*pi*10^(-7); % Magnetic permeability of vacuum.
% i0, j0, k0           % Tx loop center is at the bottom face center
% of the cell (i0,j0,k0).
% XI, YJ, ZK           % Number of cells in the x-, y-, and z-directions.
% model_EC             % Model conductivity.
% dx, dy, dz           % Grid size in the x-, y-, and z-directions.
% subloop_x,subloop_y  % Coordinate of subloop center.
% L_subloop            % Length of subloops.
% Alpha                % Coefficient in the time step size equation.
EX_subloop = zeros(XI,YJ+1,ZK+1,2); % Subloop initial field Ex.
EY_subloop = zeros(XI+1,YJ,ZK+1,2); % Subloop initial field Ey.
EZ_subloop = zeros(XI+1,YJ+1,ZK,2); % Subloop initial field Ez.
HX_subloop = zeros(XI+1,YJ,ZK,2);   % Subloop initial field Hx.
HY_subloop = zeros(XI,YJ+1,ZK,2);   % Subloop initial field Hy.
HZ_subloop = zeros(XI,YJ,ZK+1,2);   % Subloop initial field Hz.
% t0_E                % Initial time of electric field.
% t0_H                % Initial time of magnetic field.
dt = zeros(2);        % Time step size.
% integer::i,j,k      % Temporary loop variables.
% real(8)::coord_x,coord_y,coord_z % Coordinates of the node.
% real(8)::coord_r    % Distance between the grid node and Tx center.
% real(8)::u          % Temporary variable used to calculate initial field.
% real(8)::initialfield_EC % The conductivity used to 
                           % calculate initial field.
% real(8)::Txmm            % Tx magnetic moment.
%- End of header ----------------------------------------------------------

%--------------------------------------------------------------------------
% [1.0] Set the conductivity used to calculate initial field:
%--------------------------------------------------------------------------
% NOTE: The initial field is independent of the model, and can be calculated
%       using an artificial high conductivity (e.g. initialfield_EC = 100),
%       here the conductivity of cell (i0,j0,k0+1) is used by default:
initialfield_EC = model_EC(i0,j0,k0+1);
% Tx magnetic moment (current is 1 A by default):
Txmm = L_subloop^2;
%--------------------------------------------------------------------------
% [2.0] Calculate initial time:
%--------------------------------------------------------------------------
% Calculate the initial time of electric field using equation 26 in paper:
t1_E = 1.13*Permeability_vac*model_EC(i0,j0,k0)*dz(k0)^2;
% Calculate the time step size using equation 28 in paper:
dt(1) = Alpha*dz(k0)*sqrt(Permeability_vac*model_EC(i0,j0,k0)*t1_E/6);
% Calculate the initial time of magnetic field using equation 27 in paper:
t1_H = t1_E+dt(1)/2;
%--------------------------------------------------------------------------
% [3.0] Calculate subloop initial field using equations 25 in paper:
%--------------------------------------------------------------------------
% Calculate subloop initial field Ex:
for k = 2:ZK
    for j = 2:YJ
        for i = 1:XI
            [coord_x] = coord_x_calc(i,1,i0,dx)+subloop_x;
            [coord_y] = coord_y_calc(j,1,j0,dy)+subloop_y;
            [coord_z] = coord_z_calc(k,1,k0,dz);
            coord_r = sqrt(coord_y^2+coord_z^2+coord_x^2);
            u = coord_r*sqrt(2*pi*initialfield_EC/(1e7*t1_E));
            EX_subloop(i,j,k,1) = sqrt(2/pi)*Txmm*coord_y*(u/coord_r)^5 ...
                /(4*pi*initialfield_EC*exp(u^2/2));
        end
    end
end
% Calculate subloop initial field Ey:
for k = 2:ZK
    for j = 1:YJ
        for i = 2:XI
            [coord_x] = coord_x_calc(i,2,i0,dx)+subloop_x;
            [coord_y] = coord_y_calc(j,2,j0,dy)+subloop_y;
            [coord_z] = coord_z_calc(k,2,k0,dz);
            coord_r = sqrt(coord_y^2+coord_z^2+coord_x^2);
            u = coord_r*sqrt(2*pi*initialfield_EC/(1e7*t1_E));
            EY_subloop(i,j,k,1) = sqrt(2/pi)*Txmm*coord_x*(u/coord_r)^5 ...
                /(4*pi*initialfield_EC*exp(u^2/2));
        end
    end
end
% Calculate subloop initial field Ez:
EZ_subloop(2:XI,2:YJ,1:ZK,1) = 0;
% Calculate subloop initial field Hx:
for k = 1:ZK
    for j = 1:YJ
        for i = 1:XI+1
            [coord_x] = coord_x_calc(i,4,i0,dx)+subloop_x;
            [coord_y] = coord_y_calc(j,4,j0,dy)+subloop_y;
            [coord_z] = coord_z_calc(k,4,k0,dz);
            coord_r = sqrt(coord_y^2+coord_z^2+coord_x^2);
            u = coord_r*sqrt(2*pi*initialfield_EC/(1e7*t1_H));
            HX_subloop(i,j,k,1) = (Txmm*coord_z*coord_x) ...
                /(4*pi*coord_r^5)*(3*erf(u)-sqrt(2/pi) ...
                *u*(3+u^2)/exp(u^2/2));
        end
    end
end
% Calculate subloop initial field Hy:
for k = 1:ZK
    for j = 1:YJ+1
        for i = 1:XI
            [coord_x] = coord_x_calc(i,5,i0,dx)+subloop_x;
            [coord_y] = coord_y_calc(j,5,j0,dy)+subloop_y;
            [coord_z] = coord_z_calc(k,5,k0,dz);
            coord_r = sqrt(coord_y^2+coord_z^2+coord_x^2);
            u = coord_r*sqrt(2*pi*initialfield_EC/(1e7*t1_H));
            HY_subloop(i,j,k,1) = (Txmm*coord_z*coord_y) ...
                /(4*pi*coord_r^5)*(3*erf(u)-sqrt(2/pi) ...
                *u*(3+u^2)/exp(u^2/2));
        end
    end
end
% Calculate subloop initial field Hz:
for k = 1:ZK+1
    for j = 1:YJ
        for i = 1:XI
            [coord_x] = coord_x_calc(i,6,i0,dx)+subloop_x;
            [coord_y] = coord_y_calc(j,6,j0,dy)+subloop_y;
            [coord_z] = coord_z_calc(k,6,k0,dz);
            coord_r = sqrt(coord_y^2+coord_z^2+coord_x^2);
            u = coord_r*sqrt(2*pi*initialfield_EC/(1e7*t1_H));
            if coord_r < 0.001*dz(k0)
               coord_r = 0.001*dz(k0); % Avoid singularity.
            end
            HZ_subloop(i,j,k,1) = Txmm/(4*pi*coord_r^5)             ...
                *((2*coord_z^2-sqrt(coord_x^2+coord_y^2))*erf(u)    ...
                -( 2*coord_z^2-sqrt(coord_x^2+coord_y^2)*(1+u^2)) ...
                *sqrt(2/pi)/exp(u^2/2));
        end
    end
end
%--------------------------------------------------------------------------
% [4.0] SubFunctions used to calculate coordinates in equations 25 in paper:
%--------------------------------------------------------------------------
% Calculate the value of x coordinate:
function [coord_x] = coord_x_calc(x,m,i0,dx)
coord_x = 0;
if ( m == 2 || m == 3 || m == 4)
    if (x <= i0)
        for i = x:i0
            coord_x = coord_x+dx(i);
        end
        coord_x = coord_x-0.5*dx(i0);
        coord_x = -coord_x;
    elseif(x > i0)
        for i = i0:x-1
            coord_x = coord_x+dx(i);
        end
        coord_x = coord_x-0.5*dx(i0);
    end
else
    if (x <= i0)
        for i = x:i0
            coord_x = coord_x+dx(i);
        end
        coord_x = coord_x-0.5*(dx(x)+dx(i0));
        coord_x = -coord_x;
    elseif(x > i0)
        for i = i0:x
            coord_x = coord_x+dx(i);
        end
        coord_x = coord_x-0.5*(dx(i0)+dx(x));
    end
end
% Calculate the value of y coordinate:
function [coord_y] = coord_y_calc(y,m,j0,dy)
coord_y = 0;
if (m == 1 || m == 3 || m == 5 )
    if (y <= j0)
        for i = y:j0
            coord_y = coord_y+dy(i);
        end
        coord_y = coord_y-0.5*dy(j0);
        coord_y = -coord_y;
    elseif(y>j0)
        for i = j0:y-1
            coord_y = coord_y+dy(i);
        end
        coord_y = coord_y-0.5*dy(j0);
    end
else
    if (y <= j0)
        for i = y:j0
            coord_y = coord_y+dy(i);
        end
        coord_y = coord_y-0.5*(dy(y)+dy(j0));
        coord_y = -coord_y;
    elseif(y > j0)
        for i = j0:y
            coord_y = coord_y+dy(i);
        end
        coord_y = coord_y-0.5*(dy(j0)+dy(y));
    end
end
% Calculate the value of z coordinate:
function [coord_z] = coord_z_calc(z,m,k0,dz)
coord_z = 0;
if (m == 3 || m == 4 || m == 5 )
    if (z <= k0)
        for i = z:k0
            coord_z = coord_z+dz(i);
        end
        coord_z = coord_z-0.5*dz(z);
        coord_z = -coord_z;
    elseif(z>k0)
        for i = k0+1:z
            coord_z = coord_z+dz(i);
        end
        coord_z = coord_z-0.5*dz(z);
    end
else
    if (z <= k0)
        for i = z:k0
            coord_z = coord_z+dz(i);
        end
        coord_z = -coord_z;
    elseif(z > k0)
        for i = k0+1:z-1
            coord_z = coord_z+dz(i);
        end
    end
end
