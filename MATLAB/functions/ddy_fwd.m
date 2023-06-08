% function dfdy = ddy_fwd(f,dy)
%     
%     % determine field size
%     [nx,ny]     = size(f);
% 
%     % allocate return field
%     dfdy        = zeros(nx,ny);
%     
%     % forward difference
%     for i=1:nx
%         for j=1:ny-1
%             dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
%         end
%     end
%     
%     % backward difference for last point
%     j = ny;
%     for i=1:nx
%         dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
%     end
% end

function df = ddy_fwd(f,dy)

df = zeros([size(f)]);
for i=1:size(f,1)
    for k=1:size(f,3)
        for j=1:size(f,2)-1
            df(i,j,k)=(f(i,j+1,k)-f(i,j,k))/dy;
        end
        %Normal
        df(i,end,k)=(f(i,end,k)-f(i,end-1,k));
    %Periodic
       % df(i,end)=(f(i,1)-f(i,end))/dy;
    end
end


% function df = ddy_fwd(f,dy)
% 
% df = zeros([size(f)]);
% for i=1:size(f,1)
%     for j=1:size(f,2)-1
%         df(i,j)=(f(i,j+1)-f(i,j))/dy;
%     end
%     %Normal
%     df(i,end)=(f(i,end)-f(i,end-1));
% %Periodic
%    % df(i,end)=(f(i,1)-f(i,end))/dy;
% end
