%%% Orbit of LEO satellite using Kepler's equation: E - e*sin(E) = M  
% Takes into account the effects of the Earth's oblatenes
% 2020/6/4

clc; 
clear; 
close all

%%% Startup:
%     1) Run setup.m
%     2) Initialize in initSimulation.m

tic;
% Import parameters
initSimulationParams;  
physicalParams    
i = init_kepler(3); % Not to confuse with complex i

% E and M change with time 
E   = 2*atan(tan(f/2)*sqrt((1-e)/(1+e)));   % Initial eccentric anomaly
M   = E - e*sin(E);                         % Initial mean anomaly 

% Simulation timing parameters 
T           = 2*pi*a^(3/2)/sqrt(u);         % Time period of a orbit(s) 
start_time  = M*T/(2*pi);                   % Initial time since perigee(s)
stop_time   = start_time + sim_time*86400;  % (s)
time        = start_time:dt:stop_time;  

%%% States that are stored
dynamic_state = zeros(2,length(time));  % [w, W] vector for RK4 integration
keplerians    = init_kepler;            % Keplerian state in timestep t
alpha = zeros(1,length(time)-1);        % Right ascension/Longitude (rad)
delta = zeros(1,length(time)-1);        % Declination/Latitude (rad)
height = zeros(1,length(time)-1);       % Height of satellite
R = zeros(3,length(time)-1);            % Position vector in ECI or ECEF

% Initial states
dynamic_state(:,1)  = [w,W];
theta = 0; % Angle between ECI and ECEI (theta is 0 at t=0)

fprintf('Simulation in progress...');
for i_iters = 1:length(time)-1
    t = time(i_iters);
    
    % Update W and w using RK4 integration
    fn = @(t,y)keplerianDynamics(t, y, keplerians);
    dynamic_state(:,i_iters+1) = RK4(fn, dynamic_state(:,i_iters), dt, t);
   
    % Update M, E and f
    M = 2*pi*t/T;                           
    E = getEccentricAnomaly(e, M);
    f = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
    
    % Compute position vector in ECI frame 
    % theta(t) =  we*(t - start_time) 
    [R_eci,~]   = transformKeplerians2ECI(keplerians);
    theta       = we*(t - start_time);  
    R_ecef      = transformECI2ECEF(R_eci,theta);
    
    R(:,i_iters)   = R_ecef; % Change R_ecef -> R_eci to see the magic
    
    % Compute  right ascension, diclination and height of satellite.
    [alpha(i_iters), delta(i_iters)] = transformCartesian2Spherical(R_ecef);
    
    % Update Keplerian elements
    w = dynamic_state(1,i_iters+1);
    W = dynamic_state(2,i_iters+1);
    keplerians = [a;e;i;w;W;f];
end

%%% Plots ground track and 3d orbit   
plotSatelliteOrbit(alpha,delta,R)

%%% Print simualtion information
% Initial Keplerian elements
n_orbits = length(time)*dt/T; clc;
fprintf('\nInitial Keplerian elements');
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~');
fprintf('\na) Perigee              ...      = %g km', rP);
fprintf('\nb) Apogee               ...      = %g km', rA);
fprintf('\n1) Semimajor axis       ...      = %g km', init_kepler(1));
fprintf('\n2) Eccentricity         ...      = %g ', init_kepler(2));
fprintf('\n3) Inclination          ...      = %g rad', init_kepler(3));
fprintf('\n4) Periapsis argument   ...      = %g rad', init_kepler(4));
fprintf('\n5) Longitude of ascending node   = %g rad', init_kepler(5));
fprintf('\n6) True anomaly         ...      = %g rad', init_kepler(6));

% Simulation parameters
fprintf('\n\nSimulation parameters');
fprintf('\n~~~~~~~~~~~~~~~~~~~~~');
fprintf('\nTimestep         = %g secs', dt);
fprintf('\nTime period      = %g secs', T);
fprintf('\nNo. of orbits    = %g', n_orbits);
fprintf('\nSimulation time  = %g days', sim_time);
fprintf('\n\n')
toc;



