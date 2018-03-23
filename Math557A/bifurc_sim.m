function bifurc_sim
%% Uses the model and analysis of Schwartz and Smith (1983)
% to plot solutions to the seasonally-forced SEIR epidemic model
% near the period-doubling bifurction 
% for measles-specific parameter values

% initial state near the periodic solution
x0 = 0.0243;
y0 = -0.0364;
z0 = -0.0669;
X0 = [x0, y0, z0];

deltavec = [0.11,0.114856];
Xmax = zeros(length(deltavec),1);

tlast = 100; % end of solution
tspan = [0,tlast];

%% Run the period 1 solution
delta1 = 0.11;
  
[t,Xt] = ode45(@(t,X)Fx(t,X,delta1),tspan,X0);

N = find(t>= 96,1); % almost on the periodic solution 

% identify period-1 solution 
xt1 = Xt(N:end,1);
yt1 = Xt(N:end,2);
zt1 = Xt(N:end,3);
T1 = t(N:end)-t(N);



%% Run the period 2 solution 
% there is a bifurcation near delta = 0.114856

delta2 = 0.114856;


[t,Xt] = ode45(@(t,X)Fx(t,X,delta2),tspan,X0);

N = find(t>= 96,1); %almost on the periodic solution

xt2 = Xt(N:end,1);
yt2 = Xt(N:end,2);
zt2 = Xt(N:end,3);
T2 = t(N:end)-t(N);

    
    
figure();
hold on
plot(T1,xt1,'LineWidth',2)
plot(T2,xt2,'LineWidth',2,'LineStyle','--')
legend('\delta = 0.11','\delta = 0.114856')
set(gca,'FontSize',12)
xlabel('time (years)')
ylabel('X')
title('X over time')
hold off

    
figure();
hold on
plot(xt1,yt1,'LineWidth',2)
plot(xt2,yt2,'LineWidth',2,'LineStyle','--')
legend('\delta = 0.11','\delta = 0.114856')
set(gca,'FontSize',12)
xlabel('X')
ylabel('Y')
title('X-Y Phase Plane')
hold off
    
figure();
hold on; grid on;
plot3(xt1,yt1,zt1,'LineWidth',2)
plot3(xt2,yt2,zt2,'LineWidth',2,'LineStyle','--')
legend('\delta = 0.11','\delta = 0.114856')
xlabel('X')
ylabel('Y')
zlabel('Z')
view(45,26)
hold off
end






function dXdt = Fx(t,X,delta)
dXdt = zeros(3,1);

% measles parameters
mu = 0.02;
alpha = 1/0.0279;
gamma = 1/0.01;
beta0 = 1575.0;
Q = 15.73807;
epsilon = 0.29476137;
eta = Q/(Q-1);

% % flu parameters
% mu = 0.02/365;
% alpha = 1/(0.39); %0.00702
% gamma = 0.133; %0.0206
% beta0 = 0.388; %141.62
% Q = (beta0*alpha)/(mu+gamma)/(mu+alpha); % 906.739
% epsilon = 0.29;
% eta = Q/(Q-1);

x = X(1);
y = X(2);
z = X(3);

dXdt(1) = -epsilon*((eta + delta*cos(2*pi*t))*x + (1+ delta*cos(2*pi*t))*z...
            + delta*cos(2*pi*t) + x*z*(1+ delta*cos(2*pi*t)));
dXdt(2) = (mu + alpha)*(delta*cos(2*pi*t)+x*(1+delta*cos(2*pi*t)) + ...
            z*(1+ delta*cos(2*pi*t)) - y + x*z*(1+ delta*cos(2*pi*t)));
dXdt(3) = (mu + gamma)*(y-z);        

end




%% If desired, you can compare the solutions to the reduced model
function dXdt = redFx(t,X,v)

dXdt = zeros(2,1);

x = X(1);
y = X(2);

dXdt(1) = -v*y;
dXdt(2) = v*x*(1+y);

end


