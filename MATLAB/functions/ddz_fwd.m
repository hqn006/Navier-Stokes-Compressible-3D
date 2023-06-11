% function df = ddz_fwd(f,dz)
% 
% df = zeros([size(f)]);
% for i=2:size(f,1)-1
%     for j=2:size(f,2)-1
%         for k=1:size(f,3)-1
%             df(i,j,k)=(f(i,j,k+1)-f(i,j,k))/dz;
%         end
%         %Normal
%         df(size(f,3))=(f(i,j,size(f,3)) -f(i,j,size(f,3)-1))/dz;
%     %Periodic
%        % df(i,end)=(f(i,1)-f(i,end))/dy;
%     end
% end
function dfdz = ddz_fwd(f,dz)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdz        = zeros(nx,ny,nz);
    
    % forward difference
    for i=2:nx-1
        for j=2:ny-1
            for k=1:nz-1
                dfdz(i,j,k) = (f(i,j,k+1)-f(i,j,k))/dz;
            end
        end
    end
    
    % backward difference for last point
    % k = nz;
    % for i=2:nx-1
    %     for j=2:ny-1
    %         dfdz(i,j,k) = (f(i,j,k)-f(i,j,k-1))/dz;
    %     end
    % end
end