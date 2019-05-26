function [u_n] = rk4(func, u, t, dt)
        % Compute intermediary values for k
        k1 = func(dt, u);
        k2 = func(t + dt/2, u + dt/2*k1);
        k3 = func(t + dt/2, u + dt/2*k2);
        k4 = func(t + dt, u + dt*k3);
        % Compute updated values for u and t
        u_n = u + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
