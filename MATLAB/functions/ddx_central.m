% function df = ddx_central(f,dx)
% df=zeros([size(f)]);
% 
% for k=2:size(f,3)-1
%     for j=2:size(f,2)-1
%         df(1,j,k)=(-3*f(1,j,k)+4*f(2,j,k)-1*f(3,j,k))/(2*dx); %Boundary condition forward
%         for i=2:size(f,1)-1
%             df(i,j,k)=(f(i+1,j,k)-f(i-1,j,k))/(2*dx); %backward
%         end
%         df(size(f,1))=-1*(-3*f(size(f,1),j,k)+4*f(size(f,1)-1,j,k)-1*f(size(f,1)-2,j,k))/(2*dx); %Boundary condition backwards  
%     end
% end
% end

function dfdx = ddx_central(f,dx)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdx        = zeros(nx,ny,nz);
    
    % central difference
    for j=2:ny-1
        for k=2:nz-1
            for i=2:nx-1
                dfdx(i,j,k) = (f(i+1,j,k)-f(i-1,j,k))/2/dx;
            end
        end
    end
    
    % % forward difference for first point
    % i = 2;
    % for j=2:ny-1
    %     for k=2:nz-1
    %         dfdx(i,j,k) = (-3*f(i,j,k)+4*f(i+1,j,k)-f(i+2,j,k))/2/dx;
    %     end
    % end
    % 
    % % backward difference for last point
    % i = nx-1;
    % for j=2:ny-1
    %     for k=2:nz-1
    %         dfdx(i,j,k) = (3*f(i,j,k)-4*f(i-1,j,k)+f(i-2,j,k))/2/dx;
    %     end
    % end
end