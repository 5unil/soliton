% Function executes RK4 method without dispersion term
function [uj1] = rk4disp(h,dt,uj)

    k1 = dt * dispderiv(uj,h);
    k2 = dt * dispderiv((uj+(0.5*k1)),h);
    k3 = dt * dispderiv((uj+(0.5*k2)),h);
    k4 = dt * dispderiv((uj+k3),h);
    
    % uj1 is calculated and returned using previous value uj
    uj1 = uj+((1/6)*(k1+(2*k2)+(2*k3)+k4));

end

