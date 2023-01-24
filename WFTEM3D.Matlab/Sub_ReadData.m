%+ Read model and configuration parameters from the input file.
%
function [L_loop,XI,YJ,ZK,dx,dy,dz,i0,j0,k0,Rx_st,Rx_end,iter_n_max, ...
    n_subloop,Alpha,model_EC]=Sub_ReadData(inputfile)
%
% In: inputfile
% Out: L_loop, XI, YJ, ZK, dx, dy, dz, i0, j0, k0, Rx_st, Rx_end,
%      iter_n_max, n_subloop, Alpha, model_EC
% Description:
% Read model and configuration parameters from the input file.
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. <Fei Li and Jiulong Cheng>
%
% Declarations:
% inputfile       % Name of the input file.
% L_loop          % Length of Tx loop. Current is 1 A by default.
% XI, YJ, ZK      % Number of cells in the x-, y-, and z-directions.
% dx, dy, dz      % Grid size in the x-, y-, and z-directions.
% i0, j0, k0      % Tx loop center is at the bottom face center
% of the cell (i0,j0,k0).
% Rx_st, Rx_end   % Receivers are at the bottom face center of cells
% from (i0,Rx_st,k0) to (i0,Rx_end,k0).
% iter_n_max      % Maximum number of iterations.
% n_subloop       % Number of subloops.
% Alpha           % Coefficient in the time step size calculation equation.
% model_EC        % Model conductivity.
%integer::i             % Temporary loop variables.
%integer::n_submodel    % Number of model subareas.
%real(8)::submodel(0:6) % Temporary array used to read model parameters.
%- End of header ----------------------------------------------------------

data=importdata(inputfile);  % Open the input file.
L_loop=data(1);   % Read length of Tx loop (m). Current is 1 A by default.
XI=data(2);       % Read number of cells in the x direction.
YJ=data(3);       % Read number of cells in the y direction.
ZK=data(4);       % Read number of cells in the z direction.
dx(1:XI)=data(5:4+XI); % Read grid size in the x direction (m).
dy(1:YJ)=data(5+XI:4+XI+YJ); % Read grid size in the y direction (m).
dz(1:ZK)=data(5+XI+YJ:4+XI+YJ+ZK); % Read grid size in the z direction (m).
o=XI+YJ+ZK;
i0=data(o+5);           % Read position of Tx loop.
j0=data(o+6);           % ditto.
k0=data(o+7);           % ditto. 
                        % Tx loop center is at the bottom face center
                        % of the cell (i0,j0,k0).
Rx_st=data(o+8);        % Read positions of the start receivers.
Rx_end=data(o+9);       % Read positions of the end receivers.
                        % Receivers are at the bottom face center of cells
                        % from (i0,Rx_st,k0) to (i0,Rx_end,k0).
iter_n_max=data(o+10);  % Read maximum number of iterations.
n_subloop=data(o+11);   % Read number of subloops.
Alpha=data(o+12);       % Read coefficient in the time step size equation.
p=data(o+13);           % Read number of model subareas.
for i=1:1:p             % Read ranges and conductivities of subareas (S/m).
    submodel(1:7)=data(o+7+i*7:o+13+i*7);
    model_EC(submodel(1):submodel(2),submodel(3):submodel(4), ...
    submodel(5):submodel(6))=submodel(7);
end
