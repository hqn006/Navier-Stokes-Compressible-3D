function dfdy = ddy_fwd(f,dy)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny,nz);
    
    % forward difference
    for i=2:nx-1
        for k=2:nz-1
            for j=1:ny-1
                dfdy(i,j,k) = (f(i,j+1,k)-f(i,j,k))/dy;
            end
        end
    end
    
    % backward difference for last point
    % j = ny;
    % for i=2:nx-1
    %     for k=2:nz-1
    %         dfdy(i,j,k) = (f(i,j,k)-f(i,j-1,k))/dy;
    %     end
    % end
end

% function df = ddy_fwd(f,dy)
% 
% df = zeros([size(f)]);
% for i=2:size(f,1)-1
%     for k=2:size(f,3)-1
%         for j=1:size(f,2)-1
%             df(i,j,k)=(f(i,j+1,k)-f(i,j,k))/dy;
%         end
%         %Normal
%         df(i,end,k)=(f(i,end,k)-f(i,end-1,k))/dy;
%     end
% end
