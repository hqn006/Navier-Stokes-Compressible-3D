% function df = ddy_central(f,dy)
% df=zeros([size(f)]);
% 
% 
% for i=2:size(f,1)-1
%     for k=2:size(f,3)-1
%         df(i,1,k)=((-3*f(i,1,k)+4*f(i,2,k)-1*f(i,3,k))/(2*dy));
%         for j=2:size(f,2)-1
%             df(i,j,k)=(f(i,j+1,k)-f(i,j-1,k))/(2*dy); %backward
%         end
%         df(size(f,2))=-1*(-3*f(i,size(f,2),k )+4*f(i, size(f,2)-1,k)-1*f(i,size(f,2)-2,k))/(2*dy);
%     df;
%     end
% end

function dfdy = ddy_central(f,dy)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny,nz);
    
    % central difference
    for i=2:nx-1
        for k=2:nz-1
            for j=2:ny-1
                dfdy(i,j,k) = (f(i,j+1,k)-f(i,j-1,k))/2/dy;
            end
        end
    end
    
    % % forward difference for first point
    % j = 2;
    % for i=2:nx-1
    %     for k=2:nz-1
    %         dfdy(i,j,k) = (-3*f(i,j,k)+4*f(i,j+1,k)-f(i,j+2,k))/2/dy;
    %     end
    % end
    % 
    % % backward difference for last point
    % j = ny-1;
    % for i=2:nx-1
    %     for k=2:nz-1
    %         dfdy(i,j,k) = (3*f(i,j,k)-4*f(i,j-1,k)+f(i,j-2,k))/2/dy;
    %     end
    % end
end