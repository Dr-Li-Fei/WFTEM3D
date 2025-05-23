 %+ A quick and flexible 3D transient electromagnetic modeling software.
%
% WFTEM3D2.0
%
% Description:
% A quick, flexible, and extensible 3D TEM modeling open-source software. 
% Supports wire/loop sources, half/whole-space, and ground/airborne/marine/
% tunnel/borehole scenarios.
%
% Method:
% Wire is modeled as volume current.First calculate the primary field using
% a whole-space homogeneous model, then calculate the secondary field using 
% the true model. (Relevant paper will be provided in subsequent updates.)
%
% Input files:
% ".txt" file or ".dat" file,
% e.g., "Example_conductive_brick_in_a_half-space.txt"
% In which routine they are read: Sub_ReadData.m.
%
% Output files:
% "Result_time.txt" and "Result_dBz.txt".
% In which routine they are written: Sub_SecondaryField.m.
%
% Current Code Owner: <Fei Li and Jiulong Cheng>
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. Fei Li
% 2.0      01/03/25  Add wire source support(volume current-based). Fei Li
%
% Declarations:
% character::inputfile                     % Name of the input file.
%- End of header ----------------------------------------------------------

clear all;close all;clc;tic
%--------------------------------------------------------------------------
% [1.0] Read model and configuration parameters from the input file:
%--------------------------------------------------------------------------
addpath(genpath('.\Data'))
%inputfile = 'Example_half-space_homogeneous_model_with_wire_source.txt';
%inputfile = 'Example_half-space_homogeneous_model_with_loop_source.txt';
%inputfile = 'Example_whole-space_homogeneous_model_with_loop_source.txt';
inputfile = 'Example_conductive_brick_in_a_half-space.txt';
%inputfile = 'Example_airborne_TEM_model.txt';
%inputfile = 'Example_hill_model.txt';
%inputfile = 'Example_tunnel_TEM_model.txt';
[XI,YJ,ZK,dx,dy,dz,Tx,Rx,iter_n_pri,iter_n_sec,alpha_pri, ...
    alpha_sec,model_EC]=Sub_ReadData(inputfile);
%--------------------------------------------------------------------------
% [2.0] Calculate primary field:
%--------------------------------------------------------------------------
[EX,EY,EZ,HX,HY,HZ] = Sub_PrimaryField(XI,YJ,ZK,dx,dy,dz, ...
    Tx,iter_n_pri, alpha_pri);
%--------------------------------------------------------------------------
% [3.0] Calculate secondary field:
%--------------------------------------------------------------------------
[t_iteration_H1,DBZ_Rx] = Sub_SecondaryField(XI,YJ,ZK,dx,dy,dz, ...
    Rx,Tx,iter_n_sec, alpha_sec,model_EC,EX,EY,EZ,HX,HY,HZ);    
toc
%--------------------------------------------------------------------------
% [4.0] Plot results:
%--------------------------------------------------------------------------
time = load ('.\Data\Result_time.txt');
dBz = load ('.\Data\Result_dBz.txt');
figure
loglog(time,abs(dBz))
xlabel('\it{t}\rm/s','FontSize',15)
ylabel('\it{\partial}Bz/{\partial}t\rm (V/Am^2)','FontSize',15)
set(gca,'FontSize',15)
