function [rho,u,v,w,T,p,e,Et] = cons2prim(U,R,cv)
rho=squeeze(U(1,:,:,:));
u=squeeze(U(2,:,:,:))./rho;
v=squeeze(U(3,:,:,:))./rho;
w=squeeze(U(4,:,:,:))./rho;
temp = squeeze(U(5,:,:,:))./rho;
temp = temp - 1/2.*(u.^2+v.^2+w.^2);
T=temp./cv;

e=cv.*T;
p=rho.*R.*T;
Et=squeeze(U(5,:,:,:));

end