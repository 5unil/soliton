% Soliton wave breaking up into a train of solitons of different
% sizes for an initial gaussian waveform (not normal mode).
% Also shows that mass and momentum are still conserved

h = 0.1;            % Spacial step size
dt = 0.001;         % Time step size
S = 300;            % Number of discrete steps along x 
xmax = S*h;         % Maximum x value
L = 2;              % Number of steps
x = (-xmax:h:xmax); % Discretise x values
a = 2;              % Constant alpha in kdeV equation
tnum = (20*L)/dt;   % Number of time steps
j = 1;              % Counter for for loop

u   = 12*a^2*gaussmf(x,[1 0]);  % Uses initial gaussian function 
u = u';                         % Transposes u

mass = zeros(round(tsteps/M)+10,1);    % Area vector trapezium
momen = zeros(round(tsteps/M)+10,1);   % u squared terms
t = zeros(round(tsteps/M)+10,1);       % time vector

figure

subplot(2,3,1)
plot(x,u,'LineSmoothing','on');
set(gca,'fontsize',15, 'FontWeight', 'bold');
axis([-L, L, -2, 100]);
xlabel('x');
ylabel('u');
title('t = 0');
grid on;

for i=1:tsteps

    u = rk4(h,dt,u);
    % Finds area under u using trapezium rule
    mass(j,1) = (h/2)*(2*u(1)+2*u(end)+2*sum(u(2:end-1)));
    % u^2 is positive so can use usqr to find u^2 dx
    momen(j,1) = h*sum(u.^2);   
    t(j) = M*j*dt;
    j = j+1;

    
    if i == round(tsteps/5) 
        b = 0.1;
        c = 2;

    elseif i == round((2*tsteps)/5) 
        b = 0.2;
        c = 3;

    elseif i == round((3*tsteps)/5)
        b = 0.3;
        c = 4;

    elseif i == round((4*tsteps)/5)           
        b = 0.4;
        c = 5;

    elseif i == round((5*tsteps)/5)
        b = 0.5;
        c = 6;  

    else
        c = 0;
    end
    
    
    if (c >= 2) && (c <= 6)
        subplot(2,3,c)
        plot(x,u,'LineSmoothing','on');
        hold all;
        set(gca,'fontsize',15, 'FontWeight', 'bold');
        axis([-L, L, -2, 100]);
        title(['t = ', sprintf('%1.2f',b)]);
        xlabel('x');    % x-axis label
        ylabel('u');    % y-axis label
        grid on;
    end
end

figure

scatter(t,mass,4,'filled');  % Mass plot
hold all;
scatter(t,momen,4,'filled'); % Momentum plot
set(gca,'fontsize',15, 'FontWeight', 'bold');
title(['Variation of mass and momentum over time for a = ', num2str(a), ' and a = ', sprintf('%1.2f',a1)]);
xlabel('Time')               % x-axis label
ylabel('Refer to legend')    % y-axis label
grid on;

u   = 12*a^2*gaussmf(x,[1 0]);  % Uses initial gaussian function 
u = u';                         % Transposes u

figure

plot(x,u);

for i=1:tnum

    u = rk4(h,dt,u);      % rk4 calculates the next approx u
    [peak,l] = findpeaks(u,'MinPeakHeight',5); % Finds peak height
    
    if length(peak) == 1; % Looking for only one peak
        pk(i) = peak(1);               
        lcl(i) = local(1);
        times(i) = i*dt;
        
        if i>8000                       
            if (pk(i) >= 47) && (pk(i) <= 49) 
                break     % Ends when gaussian is reformed
            end        
        end
        
        times = [];
        peak = [];
        pk = [];
        lcl = [];
    end
    
    if mod(i,L) == 0
        plot(x,u,'LineSmoothing','on');
        set(gca,'fontsize',15, 'FontWeight', 'bold');
        axis([-xmax, xmax, -1, 100]);
        title('Wave breaking');
        xlabel('x')  % x-axis label
        ylabel('u')  % y-axis label
        drawnow;
    end
    
end





