% Finds the unstable pairs of values of h and dt

S = 100;
a = 20;
j = 1;
t = 0;

% Varying dt between minimum and maximum value
maxdt = 0.001;
deldt = maxdt/50;
mindt = maxdt/50;

% Varying h between minimum and maximum value
maxh = 0.15;
dh = maxh/50; 
minh = maxh/50;

% Creates a 10000 by 1 matrix of zeros to store unstable pairs
unstabh = zeros(10000,1); 
unstabdt = zeros(10000,1); 


for dt = mindt:deldt:maxdt

    for h = minh:dh:maxh

        xmax = S*h;
        T = 4*(2*xmax)/(4*a*a); % Set for four time periods
        x = (-xmax:h:xmax);     % Discretise x values
        u = 12*(a^2)*(sech(a*(x-(4*(a^2)*t))).^2); % Works out initial u
        u = u';                 % Transposes u 
        % Finds initial area under u using trapezium rule
        iArea = (h/2)*(2*u(1)+2*u(end)+2*sum(u(2:end-1)));
        tnum = T/dt;            % Number of time steps 

        for i=1:tnum
            u = rk4(h,dt,u);    % Calls rk4 function
        end
        
        % Finds final area under u using trapezium rule
        fArea = (h/2)*(2*u(1)+2*u(end)+2*sum(u(2:end-1)));

        % Pairs of h and dt are only stable if change in area < 10%
        if fArea < (0.9*iArea) || fArea > (1.1*iArea) || isnan(fArea)
            unstabh(j,1) = h;
            unstabdt(j,1) = dt;
            j = j+1;
        end

    end
end

scatter(unstabh,unstabdt,4,'filled');
set(gca,'fontsize',15, 'FontWeight', 'bold');
axis([0, maxh, 0, maxdt]);
title(['Stability Test for a = ', int2str(a)]);
xlabel('Unstable h')  % x-axis label
ylabel('Unstable dt') % y-axis label











