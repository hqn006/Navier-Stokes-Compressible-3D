function U = prim2cons(rho,u,v,w,T,cv)
U=zeros(5,size(u,1),size(u,2),size(u,3));
U(1,:,:,:)=rho;
U(2,:,:,:)=rho.*u;
U(3,:,:,:)=rho.*v;
U(4,:,:,:)=rho.*w;
U(5,:,:,:)=rho.*(cv.*T +1/2 *(u.^2+v.^2+w.^2));
end