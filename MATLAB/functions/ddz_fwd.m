function df = ddz_fwd(f,dz)

df = zeros([size(f)]);
for i=1:size(f,1)
    for j=1:size(f,2)
        for k=1:size(f,3)-1
            df(i,j,k)=(f(i,j,k+1)-f(i,j,k))/dz;
        end
        %Normal
        % BC df(size(f,3))=(f(i,j,size(f,3)) -f(i,j,size(f,3)-1))/dz;
    %Periodic
       % df(i,end)=(f(i,1)-f(i,end))/dy;
    end
end
