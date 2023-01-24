%+ Calculate initial field.
%
function [EX,EY,EZ,HX,HY,HZ,t1_E,t1_H] = Sub_InitialField(L_loop, ...
    n_subloop,i0,j0,k0,XI,YJ,ZK,dx,dy,dz,Alpha,model_EC)
%
% In: L_loop,n_subloop,i0,j0,k0,XI,YJ,ZK,dx,dy,dz,Alpha,model_EC
% Out: EX,EY,EZ,HX,HY,HZ,t1_E,t1_H
% Description:
% Calculate initial field by superimpose subloop initial fields.
%
% Method:
% See 3D Finite-difference Transient Electromagnetic Modeling with
% Whole-space Initial Field for detail.
%
% Current Code Owner: <Fei Li and Jiulong Cheng>
%
% History:
% Version    Date    Comment
% -------    ----    -------
% 1.0      01/10/21  Original code. <Fei Li and Jiulong Cheng>
%
% Code Description:
% Language: Matlab.
% Software Standards: "European Standards for Writing and
% Documenting Exchangeable Fortran 90 Code".
%
% Declarations:
% L_loop            % Length of Tx loop (m). Current is 1 A by default.
% n_subloop                  % Number of subloops.
% L_subloop                  % Length of subloops.
% subloop_x,subloop_y        % Coordinate of the subloop center.
EX = zeros(XI,YJ+1,ZK+1,2);  % Ex.
EY = zeros(XI+1,YJ,ZK+1,2);  % Ey.
EZ = zeros(XI+1,YJ+1,ZK,2);  % Ez.
HX = zeros(XI+1,YJ,ZK,2);    % Hx
HY = zeros(XI,YJ+1,ZK,2);    % Hy.
HZ = zeros(XI,YJ,ZK+1,2);    % Hz.
% EX_subloop(:,:,:,:)        % Subloop initial field Ex.
% EY_subloop(:,:,:,:)        % Subloop initial field Ey.
% EZ_subloop(:,:,:,:)        % Subloop initial field Ez.
% HX_subloop(:,:,:,:)        % Subloop initial field Hx.
% HY_subloop(:,:,:,:)        % Subloop initial field Hy.
% HZ_subloop(:,:,:,:)        % Subloop initial field Hz.
% i0,j0,k0,XI,YJ,ZK,dx,dy,dz % See Sub_InitialField_SubLoop.
% Alpha,model_EC,t0_E,t0_H   % See Sub_InitialField_SubLoop.
% integer::i,j               % Temporary loop variables.
%- End of header ----------------------------------------------------------


fprintf('Initial field is calculating...\n')
% Calculate the length of subloop:
L_subloop = L_loop/sqrt(n_subloop);
for i = 1:sqrt(n_subloop)
    for j = 1:sqrt(n_subloop)
        % Calculate the coordinates of subloop centers:
        subloop_x = (L_loop-L_subloop)/2*(-1)+L_subloop*(i-1);
        subloop_y = (L_loop-L_subloop)/2*(-1)+L_subloop*(j-1);
        % Calculate the subloop initial field:
        [EX_subloop,EY_subloop,EZ_subloop,HX_subloop,HY_subloop, ...
            HZ_subloop,t1_E,t1_H] = Sub_InitialField_SubLoop(i0,j0,k0, ...
            XI,YJ,ZK,model_EC,dx,dy,dz,subloop_x,subloop_y,L_subloop,Alpha);
        % Superimpose the subloop initial field:
        EX = EX+EX_subloop;
        EY = EY+EY_subloop;
        EZ = EZ+EZ_subloop;
        HX = HX+HX_subloop;
        HY = HY+HY_subloop;
        HZ = HZ+HZ_subloop;
    end
end
fprintf('Initial field calculation finished.\n')