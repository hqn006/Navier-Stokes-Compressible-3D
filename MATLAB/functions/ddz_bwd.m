function dfdz=ddz_bwd(f,dz)

    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdz        = zeros(nx,ny,nz);
    
    % backward difference
    for i=1:nx
        for j=1:ny
            for k=2:nz
                dfdz(i,j,k) = (f(i,j,k)-f(i,j,k-1))/dz;
            end
        end
    end

    % forward difference for first point

        k = 1;
        for i=1:nx
            for j=1:ny
                dfdz(i,j,k) = (f(i,j,k+1)-f(i,j,k))/dz;
            end
        end
% df=zeros([size(f)]);
% 
% for i=2:size(f,1)-1
%     for j=2:size(f,2)-1
%         %df(i,1)=(f(i,2)-f(i,1))/dy; %BC backward
%         %periodic
%        % df(i,1)=(f(i,1)-f(i,end))/dy;
%         %Normal
%         df(i,j,1)=(f(i,j,2)-f(i,j,1))/dz;
%         for k=2:size(f,3)
%             df(i,j,k)=(f(i,j,k)-f(i,j,k-1))/dz;
%         end
%     end
% end


end
