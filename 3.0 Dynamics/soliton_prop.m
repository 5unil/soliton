% Animation showing the propagation of a single soliton pulse for various a

h = 0.1;            % Spacial step size
dt = 0.001;         % Time step size
S = 300;            % Number of discrete steps along x 
xmax = S*h;         % Maximum x value
L = 100;            % Number of steps
x = (-xmax:h:xmax); % Discretise x values
t = 0;              % Initial time


for a=1:5           % Repeats for various a (=constant alpha)
    tnum = (4*((2*xmax)/(4*a^2)))/dt;  % Number of time steps 
    u = 12*(a^2)*(sech(a*(x-(4*(a^2)*t))).^2); % Works out initial u
    u = u';         % Transposes u 
    for i=1:tnum
        u = rk4(h,dt,u); % rk4 calculates the next approx u
        if mod(i,L)==0
            plot(x,u);
            set(gca,'fontsize',15, 'FontWeight', 'bold');
            % Axis are autoscaled
            axis([-xmax, xmax, -5, (12*a^2+2)]);
            % Title changes according to alpha value being used
            title(['Soliton propagation for a = ', int2str(a)]); 
            xlabel('x'); % x-axis label
            ylabel('u'); % y-axis label
            grid on;
            drawnow;
        end
    end
end

