%+ A quick and flexible 3D transient electromagnetic modeling software.
%
% WFTEM3D
%
% Description:
% WFTEM3D is a quick, flexible and extensible 3D finite-difference 
% transient electromagnetic modeling open-source software. It is 
% applicable to the tunnel model, the non-flat topography model, 
% the multi-scale model, etc.
%
% Method:
% First the scheme calculates the approximate initial field excited by
% whole-space magnetic dipole sources at an initial time after the current 
% is switched off. Then the scheme steps Maxwell¡¯s equations in time 
% using a staggered grid and a modified DuFort-Frankel method. See the 
% paper "3D Finite-difference Transient Electromagnetic Modeling with 
% Whole-space Initial Field" for detail.
%
% Input files:
% ".dat" file,
% e.g., "Example1_Conductive brick in a half-space.dat"
%     , "Example2_Complex conductor at a vertical contact.dat".
% In which routine they are read: Sub_ReadData.
%
% Output files:
% "Result_time.txt", "Result_dBz.txt" and "Run_time.txt".
% In which routine they are written: Sub_Iteration.
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. <Fei Li and Jiulong Cheng>
%
% Declarations:
% character::inputfile                     % Name of the input file.
%- End of header ----------------------------------------------------------

clear all;close all;clc;
tic
%--------------------------------------------------------------------------
% [1.0] Read model and configuration parameters from the input file:
%--------------------------------------------------------------------------
addpath(genpath('..\Data'))
inputfile = 'Example1_Conductive brick in a half-space.dat';
%inputfile = 'Example2_Complex conductor at a vertical contact.dat';
[L_loop,XI,YJ,ZK,dx,dy,dz,i0,j0,k0,Rx_st,Rx_end,iter_n_max,n_subloop, ...
    Alpha,model_EC] = Sub_ReadData(inputfile);
%--------------------------------------------------------------------------
% [2.0] Calculate initial field:
%--------------------------------------------------------------------------
[EX,EY,EZ,HX,HY,HZ,t1_E,t1_H] = Sub_InitialField(L_loop,n_subloop, ...
    i0,j0,k0,XI,YJ,ZK,dx,dy,dz,Alpha,model_EC);
%--------------------------------------------------------------------------
% [3.0] Iteration and output results:
%--------------------------------------------------------------------------
[t_iteration_H,DBZ_Rx] = Sub_Iteration(i0,j0,k0,XI,YJ,ZK,dx,dy,dz, ...
    Rx_st,Rx_end,iter_n_max, Alpha,model_EC,EX,EY,EZ,HX,HY,HZ,t1_E,t1_H);
% Output run-time:
dlmwrite('..\Data\Run_time.txt',['Computation finished. Run-time is ', ...
    num2str(toc),' s.'], 'delimiter', '');
%--------------------------------------------------------------------------
% [4.0] Plot results:
%--------------------------------------------------------------------------
time = load ('..\Data\Result_time.txt');
dBz = load ('..\Data\Result_dBz.txt');
if (strcmp(inputfile,'Example1_Conductive brick in a half-space.dat') == 1)
    figure
    loglog(time,abs(dBz(:,1)))
    legend('FDTD (This study)')
    xlabel('Time (s)','FontSize',14)
    ylabel('\it{\partial}Bz/{\partial}t\rm (V/Am^2)','FontSize',14)
    set(gcf, 'position', [200 200 480 460])
    set(gca,'xminortick','on')
    set(gca,'ticklength',[0.02 0.01])
    set(gca,'FontSize',14)
    xlim([9e-5 1.1e-2])
end
if (strcmp(inputfile, ...
        'Example2_Complex conductor at a vertical contact.dat') ==1 )
    figure
    loglog(time,abs(dBz(:,1)),time,abs(dBz(:,6)), ...
        time,abs((dBz(:,17)+dBz(:,18))/2),time,abs(dBz(:,34)))
    legend('(0,50,0)','(0,150,0)','(0,450,0)','(0,1050,0)')
    xlabel('Time (s)','FontSize',14)
    ylabel('|\it{\partial}Bz/{\partial}t\rm| (V/Am^2)','FontSize',14)
    set(gcf, 'position', [200 200 480 460]);
    set(gca,'xminortick','on')
    set(gca,'ticklength',[0.02 0.01])
    set(gca,'FontSize',14)
    xlim([9e-5 3.1e-2])
end