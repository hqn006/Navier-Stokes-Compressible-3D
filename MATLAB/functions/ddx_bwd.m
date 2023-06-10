function dfdx = ddx_bwd(f,dx)
    
    % determine field size
    [nx,ny,nz]     = size(f);

    % allocate return field
    dfdx        = zeros(nx,ny,nz);
    
    % backward difference
    for j=1:ny
        for k=1:nz
            for i=2:nx
                dfdx(i,j,k) = (f(i,j,k)-f(i-1,j,k))/dx;
            end
        end
    end

    % forward difference for first point

        i = 1;
        for k=1:nz
            for j=1:ny
                dfdx(i,j,k) = (f(i+1,j,k)-f(i,j,k))/dx;
            end
        end
end