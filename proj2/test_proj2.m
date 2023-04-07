% Author: Nathaniel Ruhl
% MEAM 5460 PROJ2

% global variables
global theta_tw vtip sigma cl_alpha gamma;

% Rotor properties
chord = 0.12192; % m
R = 1.524; % m
Nb = 3; % number of blades
A = pi*R^2; % m^2
sigma = Nb*chord/(pi*R); % solidity
theta_tw = -8*pi/180; % rad
cl_alpha = 2*pi;
vtip = 198.12; % m/s
gamma = 7; % enforce gamma
rho = 1.225e3; % kg/m^3
% Ib = rho*cl_alpha*chord*R^4/gamma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters that I will change
% alpha = 5*pi/180;   % rotor tilt angle (rad)
v_inf = 1*42.1952;  % m/s
nu_b = 1;

alpha_list = (pi/180).*linspace(0.5,9.0,10);
[lambda, CT, theta_0, theta_1c, theta_1s, CX] = deal(zeros(length(alpha_list),1));

for i=1:length(alpha_list)
    alpha=alpha_list(i);
    [lambda(i), CT(i), theta_0(i), theta_1c(i), theta_1s(i)] = fCT_part1(alpha, v_inf, nu_b);
    CX = CT.*vtip^2.*sin(alpha)/(0.5*v_inf^2);
end

% This function returns the relevant parameters when beta_1c=beta_1s=0 in
% problem 1
function [lambda, CT, theta_0, theta_1c, theta_1s] = fCT_part1(alpha, v_inf, nu_b)
global theta_tw vtip sigma cl_alpha gamma;

syms theta_0 theta_1c theta_1s lambda CT eqn

beta_0 = 0;
beta_1c = 0;
beta_1s = 0;

mu = v_inf*cos(alpha)/vtip;

% System Ab=c given in Homework assignment
A = sym(zeros(3));
c = sym(zeros(3,1));

% First row of A
A(1,1) = 8*nu_b^2/gamma;
A(1,2) = 0;
A(1,3) = 0;

% Second row of A
A(2,1) = (4/3)*mu;
A(2,2) = 8*(nu_b^2-1)/gamma;
A(2,2) = 1+(mu^2/2);

% Third row of A
A(3,1) = 0;
A(3,2) = -(1-(mu^2/2));
A(3,3) = 8*(nu_b^2-1)/gamma;

% LHS b vector
b = [beta_0; beta_1c; beta_1s];

% RHS c vector
c(1) = (1+mu^2)*theta_0 + ((4/5)+(2/3)*mu^2)*theta_tw + (4/3)*mu*theta_1s - (4/3)*lambda;
c(2) = (1+(mu^2/2))*theta_1c;
c(3) = (8/3)*mu*theta_0+2*mu*theta_tw+(1+(3/2)*mu^2)*theta_1s-2*mu*lambda;

% Assemble lhs of eqns 1 through 5
eqn_list={0,0,0,0,0};
for i=1:3
    eqn = 0;
    for j=1:3
        eqn = eqn + A(i,j)*b(j);
    end
    eqn_list{i} = eqn - c(i);
end
eqn_list{4} = lambda - mu*tan(alpha)+CT/(2*sqrt(mu^2+lambda^2));
eqn_list{5} = CT - (sigma*cl_alpha/2)*((theta_0/3)*(1+1.5*mu^2)+0.25*theta_tw*(1+mu^2)+0.5*mu*theta_1s-0.5*lambda);
vars_to_solve = [CT lambda theta_0 theta_1c theta_1s];

sol = solve(eqn_list{1}==0,eqn_list{2}==0,eqn_list{3}==0,eqn_list{4}==0,eqn_list{5}==0,vars_to_solve);
disp(sol.lambda);

% return the average of the two correct-looking solutions
lambda = 0.5*(sol.lambda(3)+sol.lambda(4));
CT = 0.5*(sol.CT(3)+sol.CT(4));
theta_0 = 0.5*(sol.theta_0(3)+sol.theta_0(4));
theta_1c = 0.5*(sol.theta_1c(3)+sol.theta_1c(4));
theta_1s = 0.5*(sol.theta_1s(3)+sol.theta_1s(4));
end


