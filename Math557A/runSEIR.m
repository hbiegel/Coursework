%% SEIR model 
% Implemenation of the SEIR model (Schwartz and Smith, 1983)
% with measles-specific disease parameters
function runSEIR
% 
% s0 = 0.99;
% e0 = 0;
% i0 = 0.01;
% r0 = 0;


% a point as close as possible to the periodic orbit of the 
% seasonally-forced system
x0 = 0.0024;
y0 = 0.0027; 
z0 = -6.2317e-04; 

% Steady state of the basic SEIR model corresponding to the 
% periodic solution of the measles-specific seasonally-forced 
% SEIR model
s0 = 0.0635*(1+x0);
e0 = 5.2225e-04*(1+y0);
i0 = 1.8715e-04*(1+z0);
r0 = 1-s0-e0-i0; %0.9358;

X0 = [s0,e0,i0,r0];

tspan = [0,10];

[t,Xt] = ode45(@(t,X)SEIRfunc(t,X),tspan,X0);

St = Xt(:,1);
Et = Xt(:,2);
It = Xt(:,3);
Rt = Xt(:,4);



figure();
hold on
plot(t,St,'LineWidth',3)
plot(t,Et,'LineWidth',3)
plot(t,It,'LineWidth',3)
plot(t,Rt,'LineWidth',3)

legend('S','E','I','R')
xlabel('Time (years)')
ylabel('Fraction of Population')
%ylim([0,1])
set(gca,'FontSize',14)
hold off



end




function dXdt = SEIRfunc(t,x)
% measles parameters
mu = 0.02;
alpha = 1/0.0279;
gamma = 1/0.01;
beta0 = 1575.0;
%Q = 15.73807;
%epsilon = 0.29476137;
%eta = Q/(Q-1);

dXdt = zeros(4,1);
S = x(1);
E = x(2);
I = x(3);
R = x(4);


dXdt(1) = mu - mu*S - beta0*I*S;
dXdt(2) = beta0*I*S - (mu + alpha)*E;
dXdt(3) = alpha*E - (mu + gamma)*I;
dXdt(4) = gamma*I - mu*R;


end
