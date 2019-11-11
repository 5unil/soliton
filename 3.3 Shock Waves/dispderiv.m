% Function calculates derivative at u without dispersive term 
function [der] = dispderiv(u,h)

    % u makes use of periodic boundary conditions
    u = [u(end-1);u(end);u;u(1);u(2)];
    % nlin is the non linear term
    nlin = -(1/(4*h))*((u(4:end-1)).^2-(u(2:end-3)).^2);
    der = nlin;
 
end

