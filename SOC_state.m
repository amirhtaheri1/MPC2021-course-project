function dxdt = SOC_state(x,u,Ts)
% Continuous-time state equations for the exothermic CSTR
%
% global PBpb;
% The states of the CSTR model are:
%
%   x(1) = T        Reactor temperature [K]
%   x(2) = CA       Concentration of A in reactor tank [kgmol/m^3]
%   x(3) = Dist     State of unmeasured output disturbance
%
% The inputs of the CSTR model are:
%
%   u(1) = CA_i     Concentration of A in inlet feed stream [kgmol/m^3]
%   u(2) = T_i      Inlet feed stream temperature [K]
%   u(3) = T_c      Jacket coolant temperature [K]
%   u(4) = WN       White noise

% Copyright 1990-2018 The MathWorks, Inc.

% states
SOC_pb = x(1);
SOC_li = x(2);

% inputs
P_grid = u(1);
PB_li = u(2);
P_load = u(3);
P_gen=u(4);
P_noise=u(5);
% parameters
% A = eye(2);
% B = [0.0936 0.0936 0;0 0.0752 0];
% D = [0.0936;0];
A = eye(2);
B = [-0.04352 -0.04352 0.04352;0 0.15299 0];
D = [-0.04352;0];
dxdt = zeros(2,1);
% state equations
% dxdt = A*([SOC_pb;SOC_li]) + B*([P_grid;PB_li;P_load]) + D*P_gen;
% dxdt = A*([SOC_pb;SOC_li]) + B*([P_grid;PB_li;P_load]) + D*P_gen;
% dxdt=A*([SOC_pb;SOC_li])+[-0.04352,0;-0.15299,0]*[PBpb;PB_li];
global PBpb
PBpb=u(1)+u(4)-u(3)-u(2);
dxdt(1)=SOC_pb+0.04352*PBpb;
dxdt(2)=SOC_li+0.15299*PB_li;
