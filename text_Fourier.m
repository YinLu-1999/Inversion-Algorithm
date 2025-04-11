function [dBxdx, dBxdy, dBxdz, dBydy, dBydz] = text_Fourier(B,dx,dy)
[nx, ny] = size(B);
%观测面的数据一定要把源点信息全覆盖
wx = (-nx/2:nx/2-1)/(nx*dx); 
kx = 2*pi*fftshift(wx);
wy = (-ny/2:ny/2-1)/(ny*dy); 
ky = 2*pi*fftshift(wy);
% K=sqrt(kx.^2 + ky.^2) ;
[KX, KY] = meshgrid(kx, ky);

K2 =sqrt(KX.^2 + KY.^2) ;

K2(K2==0) = eps; % Avoid division by zero

% Take 2D FFT of measured Bz
F_Bz = fft2(B);

% Calculate Fourier transforms of field components using Maxwell equations
% From equation (16) in the paper
F_Bx = -1i * KY .* F_Bz ./ K2;
F_By = -1i * KX .* F_Bz ./ K2;

% Calculate Fourier transforms of field derivatives
% For Bx derivatives
F_dBx_dx = 1 * KX.^2 .* F_Bz./K2;
F_dBx_dy = 1 * KX .*KY.* F_Bz./K2;
F_dBx_dz = 1i * KX .* F_Bz;

% For By derivatives

F_dBy_dy = 1 * KY.^2 .* F_Bz./K2;
F_dBy_dz = 1i * KY .* F_Bz;
% Inverse FFT for derivatives
dBx_dx = real(ifft2(F_dBx_dx));
dBx_dy = real(ifft2(F_dBx_dy));
dBx_dz = real(ifft2(F_dBx_dz));
dBy_dy = real(ifft2(F_dBy_dy));
dBy_dz = real(ifft2(F_dBy_dz));

dBxdx = dBx_dx';
dBxdy = dBx_dy';
dBxdz = dBx_dz';
dBydy = dBy_dy';
dBydz = dBy_dz';