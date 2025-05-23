%+ Read model and configuration parameters from the input file.
%
function [XI,YJ,ZK,dx,dy,dz,Tx,Rx,iter_n_pri,iter_n_sec,alpha_pri, ...
    alpha_sec,model_EC]=Sub_ReadData(inputfile)
%
% In: inputfile
% Out: XI,YJ,ZK,dx,dy,dz,Tx,Rx,iter_n_pri,iter_n_sec,alpha_pri,
%      alpha_sec,model_EC
% Description:
% Read model and configuration parameters from the input file.
%
% Current Code Owner: <Fei Li and Jiulong Cheng>
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. Fei Li
% 2.0      01/03/25  Add support for annotated input file. Fei Li
%
% Declarations:
% inputfile      % Name of the input file.
% XI, YJ, ZK     % Number of cells in the x-, y-, and z-directions.
% dx, dy, dz     % Grid size in the x-, y-, and z-directions.
% Tx             % Tx loop diagonal points / ground wire endpoints.
% Rx             % Rx points: start/end/interval.
% iter_n_pri     % Iteration number of primary field.
% iter_n_sec     % Iteration number of secondary field.
% alpha_pri      % Coefficient of time step size for primary field.
% alpha_sec      % Coefficient of time step size for secondary field.
% model_EC       % Model conductivity.
% integer::i     % Temporary variable.
% integer::n_submodel    % Number of subdomains in the model.
% real(8)::submodel(0:6) % Temporary array used to read model parameters.
%- End of header ----------------------------------------------------------

fid = fopen(inputfile, 'r'); % Open the input file.
if ~exist(inputfile, 'file') % Check file existence.
    error('ERROR: Input file "%s" not found!', inputfile);
end
%Comment lines start with #, data is separated by spaces, ',', or ';': 
dataCell = textscan(fid, '%f', 'CommentStyle', '#', 'Delimiter', ' ,;');
fclose(fid);data = [dataCell{:}];% Convert cell array to numeric values.

XI=data(1);                  % Read number of cells in the x direction.
YJ=data(2);                  % Read number of cells in the y direction.
ZK=data(3);                  % Read number of cells in the z direction.
dx(1:XI)=data(4:3+XI);       % Read grid size in the x direction (m).
dy(1:YJ)=data(4+XI:3+XI+YJ); % Read grid size in the y direction (m).
dz(1:ZK)=data(4+XI+YJ:3+XI+YJ+ZK); % Read grid size in the z direction (m).
o=XI+YJ+ZK;
Tx=data(o+4:o+9);            % Read position of Tx.
Rx=data(o+10:o+18);          % Read position of Rx.
iter_n_pri=data(o+19);       % Read iteration number of primary field.
iter_n_sec=data(o+20);       % Read iteration number of secondary field.
alpha_pri=data(o+21);  % Read coefficient of time step size for pri field.
alpha_sec=data(o+22);  % Read coefficient of time step size for sec field.
n_submodel=data(o+23); % Read number of model subdomains.
for i=1:1:n_submodel   % Read ranges and conductivities of subdomains (S/m).
    submodel(1:6)=data(o+17+i*7:o+22+i*7);
    submodel_EC=data(o+23+i*7);
    model_EC(submodel(1):submodel(2),submodel(3):submodel(4), ...
        submodel(5):submodel(6))=submodel_EC;
end
