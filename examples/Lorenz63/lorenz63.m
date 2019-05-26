function du = lorenz63(t,u,sigma,rho,beta)
	du1 = -sigma*u(:,1) + sigma*u(:,2);
	du2 = -u(:,1).*u(:,3) + rho*u(:,1) - u(:,2);
	du3 = u(:,1).*u(:,2) - beta*u(:,3);
	du  = [du1, du2, du3];
end
