function dfdx = ddx_bwd(f,dx)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdx        = zeros(nx,ny,nz);
    
    % backward difference
    for i=2:nx
        for j=1:ny 
            for k=1:nz
                dfdx(i,j,k) = (f(i,j,k)-f(i-1,j,k))/dx;
            end
        end
    end

    % forward difference for first point
    %for k=1:nz
    %    i = 1;
    %    for j=1:ny
    %        dfdx(i,j,k) = (f(i+1,j,k)-f(i,j,k))/dx;
    %    end
    %end
    
    %     % assuming periodicity  (left boudary)
    %     i = 1;
    %     for j=1:ny
    %         dfdx(i,j) = (f(i,j)-f(end,j))/dx;
    %     end
end
