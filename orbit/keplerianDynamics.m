function [state_dot] = keplerianDynamics(~,~,kep_state)
%%% Equations of W & w due to the effect of the Earth's oblatenes
%
% Inputs:
%   kep_state = 
%
% Inputs:
%   kep_state = State vector of Keplerian elements: [a,e,i,w,W,f]
%
% Reference:
%   Howard D. Curtis - Orbital Mechanics For Engineering Students
%
% 2020/6/4

% Unpack Keplerian parameters
a =  kep_state(1);  
e =  kep_state(2);  
i =  kep_state(3);  

% Import physical parameters
physicalParams

fac     = -3/2*sqrt(u)*J2*Re^2/(1-e^2)^2/a^(7/2); % Common factor 
W_dot   = fac*cos(i);                               
w_dot   = fac*(5/2*sin(i)^2 - 2);

state_dot = [W_dot, w_dot]';
end