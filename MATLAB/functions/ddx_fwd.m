%1.6
% function df = ddx_fwd(f,dx)
% 
% df=zeros([size(f)]);
% for k=2:size(f,3)-1
%     for j=2:size(f,2)-1
%         for i=1:size(f,1)-1
%             df(i,j,k)=(f(i+1,j,k)-f(i,j,k))/dx; %forward
%         end
%         %Periodic BC
%         %df(end,j)=(f(1,j)-f(end,j))/dx;
%         %Normal
%         df(size(f,1))=(f(size(f,1),j,k)-f(size(f,1)-1,j,k))/dx; %BOundary condition backwards  
%     end
% end
% 
% end

function dfdx = ddx_fwd(f,dx)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdx        = zeros(nx,ny,nz);
    
    % forward difference
    for j=1:ny
        for k=1:nz               
            for i=1:nx-1
                dfdx(i,j,k) = (f(i+1,j,k)-f(i,j,k))/dx;
            end
        end
    end
    
    % backward difference for last point
    i = nx;
    for j=1:ny
        for k=1:nz
            dfdx(i,j,k) = (f(i,j,k)-f(i-1,j,k))/dx;
        end
    end
end