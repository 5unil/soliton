% Function calculates derivative at u with diffusive term replacing
% dispersive term
function [der] = diffderiv(u,h)

    % u makes use of periodic boundary conditions
    u = [u(end-1);u(end);u;u(1);u(2)];
    % nlin is the non linear term
    nlin = -(1/(4*h))*((u(4:end-1)).^2-(u(2:end-3)).^2);
    % diff is the diffusive term
    diff = (1/h*h)*(u(4:end-1)-(2*u(3:end-2))+u(2:end-3));
    der = nlin + diff;
 
end

