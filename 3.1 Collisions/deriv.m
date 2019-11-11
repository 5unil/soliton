% Function calculates derivative at u 
function [der] = deriv(u,h)

    % u makes use of periodic boundary conditions
    u = [u(end-1);u(end);u;u(1);u(2)];
    % nlin is the non linear term
    nlin = -(1/(4*h))*((u(4:end-1)).^2-(u(2:end-3)).^2);
    % disp is the dispersive term
    disp = -(1/(2*h*h*h))*(u(5:end)-(2*u(4:end-1))+(2*u(2:end-3))-u(1:end-4));
    der = nlin + disp;
 
end

