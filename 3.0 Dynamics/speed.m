% Measures the soliton speed for different heights (proportional to a)

h = 0.1;            % Spacial step size
dt = 0.001;         % Time step size
tnum = 3/dt;        % Number of time steps 
S = 300;            % Number of discrete steps along x 
xmax = S*h*10;      % Maximum x value
L = 2;              % Number of steps
umax(1,1) = 0;      % Set initial maximum of u at x = 0 
t(1,1) = 0;         % Set initial maximum of u at t = 0 


for a = 1:5         % Repeats for various a (=constant alpha)
  
    t = 0;          % Initial time
    x = (-xmax:h:xmax);  % Discretise x values
    u = 12*(a^2)*(sech(a*(x-(4*(a^2)*t))).^2); % Works out initial u
    u = u';         % Transposes u 
    t = zeros(S/L,1);    % Creates an S/M by 1 matrix of zeros
    umax = zeros(S/L,1); % Creates an S/M by 1 matrix of zeros
    j = 1;
    
    for i=1:tnum
        u = rk4(h,dt,u); % rk4 calculates the next approx u     
        if mod(i,L) == 0
            j = j+1;
            [q,p] = max(u);   % Locate the maximum values of u
            umax(j,1) = x(p); % Store value of x at max u
            t(j,1) = j*L*dt;  % Store value of t at max u
        end
    end    
    
    scatter(t,umax,4,'filled');
    axis([0 2.5 0 50]);
    set(gca,'fontsize',15, 'FontWeight', 'bold');
    title('Soliton speed for a = 1 to 5'); 
    xlabel('Time');                    % x-axis label
    ylabel('Value of x at maximum u'); % y-axis label
    grid on; 
    hold all;  
    
end








