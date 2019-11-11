% Function executes RK4 with diffusive term replacing dispersive term
function [uj1] = rk4diff(h,dt,uj)

    k1 = dt * diffderiv(uj,h);
    k2 = dt * diffderiv((uj+(0.5*k1)),h);
    k3 = dt * diffderiv((uj+(0.5*k2)),h);
    k4 = dt * diffderiv((uj+k3),h);
    
    % uj1 is calculated and returned using previous value uj
    uj1 = uj+((1/6)*(k1+(2*k2)+(2*k3)+k4));

end

