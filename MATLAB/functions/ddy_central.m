function df = ddy_central(f,dy)
df=zeros([size(f)]);


for i=1:size(f,1)
    for k=1:size(f,3)
        df(i,1,k)=((-3*f(i,1,k)+4*f(i,2,k)-1*f(i,3,k))/(2*dy));
        for j=2:size(f,2)-1
            df(i,j,k)=(f(i,j+1,k)-f(i,j-1,k))/(2*dy); %backward
        end
        df(size(f,2))=-1*(-3*f(i,size(f,2),k )+4*f(i, size(f,2)-1,k)-1*f(i,size(f,2)-2,k))/(2*dy);
    df;
    end
end
