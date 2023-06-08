function df=ddy_bwd(f,dy)

df=zeros([size(f)]);

for i=1:size(f,1)
    for k=1:size(f,3)
        %df(i,1)=(f(i,2)-f(i,1))/dy; %BC backward
        %periodic
       % df(i,1)=(f(i,1)-f(i,end))/dy;
        %Normal
        df(i,1,k)=(f(i,2,k)-f(i,1,k))/dy;
        for j=2:size(f,2)
            df(i,j,k)=(f(i,j,k)-f(i,j-1,k))/dy;
        end
    end
end


end
