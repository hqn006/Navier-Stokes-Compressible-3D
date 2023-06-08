function df = ddz_central(f,dz)
df=zeros([size(f)]);

for i=1:size(f,1)
    for j=1:size(f,2)
        df(i,j,1)=(-3*f(i,j,1)+4*f(i,j,2)-1*f(i,j,3))/(2*dz); %Boundary condition forward
        for k=2:size(f,3)-1
            df(i,j,k)=(f(i,j,k+1)-f(i,j,k-1))/(2*dz); %backward
        end
        df(size(f,3))=-1*(-3*f(i,j,size(f,3))+4*f(i,j,size(f,3)-1)-1*f(i,j,size(f,3)-2))/(2*dz); %Boundary condition backwards  
    end
end
end