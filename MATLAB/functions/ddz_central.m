% function df = ddz_central(f,dz)
% df=zeros([size(f)]);
% 
% for i=2:size(f,1)-1
%     for j=2:size(f,2)-1
%         df(i,j,1)=(-3*f(i,j,1)+4*f(i,j,2)-1*f(i,j,3))/(2*dz); %Boundary condition forward
%         for k=2:size(f,3)-1
%             df(i,j,k)=(f(i,j,k+1)-f(i,j,k-1))/(2*dz); %backward
%         end
%         df(size(f,3))=-1*(-3*f(i,j,size(f,3))+4*f(i,j,size(f,3)-1)-1*f(i,j,size(f,3)-2))/(2*dz); %Boundary condition backwards  
%     end
% end
% end

function dfdz = ddz_central(f,dz)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdz        = zeros(nx,ny,nz);
    
    % central difference
    for i=2:nx-1
        for j=2:ny-1
            for k=2:nz-1
                dfdz(i,j,k) = (f(i,j,k+1)-f(i,j,k-1))/2/dz;
            end
        end
    end
    
    % % forward difference for first point
    % k = 2;
    % for i=2:nx-1
    %     for j=2:ny-1
    %         dfdz(i,j,k) = (-3*f(i,j,k)+4*f(i,j,k+1)-f(i,j,k+2))/2/dz;
    %     end
    % end
    % 
    % % backward difference for last point
    % k = nz-1;
    % for i=2:nx-1
    %     for j=2:ny-1
    %         dfdz(i,j,k) = (3*f(i,j,k)-4*f(i,j,k-1)+f(i,j,k-2))/2/dz;
    %     end
    % end
end