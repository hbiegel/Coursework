function simModel

% delta = 0.01; %0.1;
% x0 = 1;
% y0 = 1;
% z0 = 1;

deltavec = [0.01, 0.05, 0.1];

x0vec = [0.0024, 0.0122, 0.0243];
y0vec = [0.0027, -0.0011, -0.0364];
z0vec = [-6.2317e-04, -0.0172, -0.0669];
% x0 = 0.0024;%[-2.7619e-04];
% y0 = 0.0027; %[0.0158];
% z0 = -6.2317e-04; %[-8.6072e-04];



for i = 10:10
    
    currFig = figure; set(currFig, 'Visible', 'on');
        hold on;
    
    for j = 1:3
        
        delta = deltavec(j);
        x0 = x0vec(j);
        y0 = y0vec(j);
        z0 = z0vec(j);
        
        
        tlast = i*0.1;
        %tlast = 1;

        tspan = [0,tlast];
        X0 = [x0, y0, z0];
        [t,Xt] = ode45(@(t,X)Fx(t,X,delta),tspan,X0);

        N = 1;
        xt = Xt(N:end,1);
        yt = Xt(N:end,2);
        zt = Xt(N:end,3);
        T = t(N:end)-t(N);

        plot(zt,xt,'k.','MarkerSize',12);

    end
    
        xlim([-0.6,0.6])
        ylim([-.03,.03])
        xlabel('z')
        ylabel('x')
        legend('delta = 0.01','delta = 0.05','delta = 0.1')
        set(gca,'FontSize',12)
        title('X-Z phase plane for period1 solutions')
        hold off;
        %  
    
end


% 
% figure();
% hold on; box on;
% plot(T,xt,'LineWidth',2)
% plot(T,yt,'LineWidth',2)
% plot(T,zt,'LineWidth',2)
% legend('x','y','z')
% xlabel('Time')
% ylabel('x,y,z')
% hold off;
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


function dXdt = redFx(t,X,v)

dXdt = zeros(2,1);

x = X(1);
y = X(2);

dXdt(1) = -v*y;
dXdt(2) = v*x*(1+y);

end



