function surf_cs(A)

[nx,~] = size(A);

B = nan*zeros(4*nx,3*nx);

B(nx+1:2*nx,nx+1:2*nx) = A(:,:,3);
B(2*nx+1:3*nx,nx+1:2*nx) = A(:,:,2);
B(1:nx,nx+1:2*nx)  = A(:,:,4);
B(nx+1:2*nx,1:nx) = A(:,:,5);
B(nx+1:2*nx,2*nx+1:3*nx) = A(:,:,6);
B(3*nx+1:4*nx,nx+1:2*nx) = A(:,:,1);
surf(B'); shading flat

end
