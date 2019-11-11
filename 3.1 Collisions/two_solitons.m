% Animation showing the propagation of two solitons at different speeds
% Finds out whether mass and momentum of the waves are both conserved

h = 0.1;            % Spacial step size
dt = 0.001;         % Time step size
S = 300;            % Number of discrete steps along x 
xmax = S*h;         % Maximum x value
L = 80;             % Number of steps
x = (-xmax:h:xmax); % Discretise x values
a = 1;              % Constant alpha in kdeV equation
t = 0;              % Initial time
j = 1;              % Counter for for loop
tnum = (4*((2*xmax)/(4*a^2)))/dt;  % Number of time steps 

% u1 and u2 represent initial u values for 2 solitons 
% at different speeds and initial positions 
u1 = 12*(a^2)*(sech(a*(x-(4*(a^2)*t))).^2); 
a1 = a+0.3;
x1 = x+4;
u2 = 12*(a1^2)*(sech(a1*(x1-(4*(a1^2)*t))).^2);   
u = (u1+u2)'; % Stores transpose of initial u for 2 solitons

% Creates a matrix of zeros to store mass and momentum values
mass = zeros(round(tnum/L),1); 
momen = zeros(round(tnum/L),1); 
t = zeros(round(tnum/L),1);    

plot(x,u);

for i=1:tnum

    u = rk4(h,dt,u); % rk4 calculates the next approx u

    if mod(i,L) == 0 
        plot(x,u,'LineSmoothing','on');
        set(gca,'fontsize',15, 'FontWeight', 'bold');
        axis([-xmax, xmax, -1, (12*a^2+10)]);
        % sprintf converts floating num to string
        title(['Soliton collision for a = ', num2str(a), ' and a = ', sprintf('%1.2f',a1)]);
        xlabel('x') % x-axis label
        ylabel('u') % y-axis label  
        grid on;
        drawnow;
        
        % Finds area under u using trapezium rule
        mass(j,1) = (h/2)*(2*u(1)+2*u(end)+2*sum(u(2:end-1)));
        % u^2 is positive so can use usqr to find u^2 dx
        momen(j,1) = h*sum(u.^2); 
        t(j) = L*dt*j;
        j = j+1;
    end
    
end

scatter(t,mass,4,'filled');  % Mass plot
hold all;
scatter(t,momen,4,'filled'); % Momentum plot
set(gca,'fontsize',15, 'FontWeight', 'bold');
title(['Variation of mass and momentum over time for a = ', num2str(a), ' and a = ', sprintf('%1.2f',a1)]);
xlabel('Time')               % x-axis label
ylabel('Refer to legend')    % y-axis label
grid on;

