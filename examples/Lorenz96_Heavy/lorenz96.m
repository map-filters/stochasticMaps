function du = lorenz96(t,u,F)
    u_ip1 = circshift(u,-1,2);
    u_im1 = circshift(u,1,2);
    u_im2 = circshift(u,2,2);
    du = (u_ip1 - u_im2).*u_im1 - u + F;
end