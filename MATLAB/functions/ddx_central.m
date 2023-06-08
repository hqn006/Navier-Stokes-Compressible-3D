function df = ddx_central(f,dx)
df=zeros([size(f)]);

for k=1:size(f,3)
    for j=1:size(f,2)
        df(1,j,k)=(-3*f(1,j,k)+4*f(2,j,k)-1*f(3,j,k))/(2*dx); %Boundary condition forward
        for i=2:size(f,1)-1
            df(i,j,k)=(f(i+1,j,k)-f(i-1,j,k))/(2*dx); %backward
        end
        df(size(f,1))=-1*(-3*f(size(f,1),j,k)+4*f(size(f,1)-1,j,k)-1*f(size(f,1)-2,j,k))/(2*dx); %Boundary condition backwards  
    end
end
end