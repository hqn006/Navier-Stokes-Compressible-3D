function dfdy=ddy_bwd(f,dy)

    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdy       = zeros(nx,ny,nz);
    
    % backward difference
    for i=1:nx
        for k=1:nz
            for j=2:ny
                dfdy(i,j,k) = (f(i,j,k)-f(i,j-1,k))/dy;
            end
        end
    end

    % forward difference for first point

        j = 1;
        for i=1:nx
            for k=1:nz
                dfdy(i,j,k) = (f(i,j+1,k)-f(i,j,k))/dy;
            end
        end
% df=zeros([size(f)]);
% 
% for i=2:size(f,1)-1
%     for k=2:size(f,3)-1
%         %df(i,1)=(f(i,2)-f(i,1))/dy; %BC backward
%         %periodic
%        % df(i,1)=(f(i,1)-f(i,end))/dy;
%         %Normal
%         df(i,1,k)=(f(i,2,k)-f(i,1,k))/dy;
%         for j=2:size(f,2)
%             df(i,j,k)=(f(i,j,k)-f(i,j-1,k))/dy;
%         end
%     end
% end


end
