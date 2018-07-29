%% REFERENCE
% https://en.wikipedia.org/wiki/Newton%27s_method

%% 
% COST FUNCTION
% x^* = argmin_x { 1/2 * || A(X) - Y ||_2^2 + lambda/2 * ( || D_x(X) ||_2^2 + || D_y(X) ||_2^2 ) }
% 
% Newton Method
% x^(k+1) = x^k - f(x^k) / f'(x^k)
%
% s.t.  f(x) = 0;
%       f'(x) = a( f(x) )/ax

%%
clear ;
close all;
home;

%% GPU Processing
% If there is GPU device on your board, 
% then isgpu is true. Otherwise, it is false.
bgpu    = false;
bfig    = true;

%%  SYSTEM SETTING
N       = 512;
VIEW    = 360;
THETA   = linspace(0, 180, VIEW + 1);   THETA(end) = [];

R       = @(x) radon(x, THETA);
RT      = @(y) iradon(y, THETA, 'none', N)/(pi/(2*length(THETA)));
RINV    = @(y) iradon(y, THETA, N);

%% DATA GENERATION
load('XCAT512.mat');
x       = imresize(double(XCAT512), [N, N]);
p       = R(x);
x_full  = RINV(p);

%% LOW-DOSE SINOGRAM GENERATION
i0     	= 5e4;
pn     	= exp(-p);
pn     	= i0.*pn;
pn     	= poissrnd(pn);
pn      = max(-log(max(pn,1)./i0),0);

y       = pn;
x_low   = RINV(y);

%% NEWTON METHOD INITIALIZATION
LAMBDA  = 1e2;
A0      = @(x)	(RT(R(x) - y) + LAMBDA*Dxt(Dx(x)) + LAMBDA*Dyt(Dy(x)));
A1      = @(x)  (RT(R(ones(size(x)))));

CG_A    = @(x) 	(LAMBDA*(Dxt(Dx(x)) + Dyt(Dy(x))) + RHO*x);

x0      = zeros(size(x));
% x0      = x_low;
niter   = 2.5e2;

L2              = @(x) power(norm(x, 'fro'), 2);
COST.equation   = '1/2 * || A(X) - Y ||_2^2 + lambda/2 * ( || D_x(X) ||_2^2 + || D_y(X) ||_2^2 )';
COST.function	= @(x) 1/2 * L2(R(x) - y) + LAMBDA/2 * (L2(Dx(x)) + L2(Dy(x)));


%% RUN NEWTON METHOD
if bgpu
    y  = gpuArray(y);
    x0 = gpuArray(x0);
end

[x_newton, obj]	= Newton(A0, A1, x0, niter, COST, bfig);

%% CALCUATE QUANTIFICATION FACTOR 
x_low           = max(x_low, 0);
x_newton        = max(x_newton, 0);
nor             = max(x(:));

mse_x_low       = immse(x_low./nor, x./nor);
mse_x_newton    = immse(x_newton./nor, x./nor);

psnr_x_low      = psnr(x_low./nor, x./nor);
psnr_x_newton   = psnr(x_newton./nor, x./nor);

ssim_x_low      = ssim(x_low./nor, x./nor);
ssim_x_newton	= ssim(x_newton./nor, x./nor);

%% DISPLAY
wndImg  = [0, 0.03];

figure(1); 
colormap(gray(256));

suptitle('Newton Method');
subplot(231);   imagesc(x,          wndImg);	axis image off;     title('ground truth');
subplot(232);   imagesc(x_full,     wndImg);  	axis image off;     title(['full-dose_{FBP, view : ', num2str(VIEW) '}']);
subplot(234);   imagesc(x_low,      wndImg);  	axis image off;     title({['low-dose_{FBP, view : ', num2str(VIEW) '}'], ['MSE : ' num2str(mse_x_low, '%.4e')], ['PSNR : ' num2str(psnr_x_low, '%.4f')], ['SSIM : ' num2str(ssim_x_low, '%.4f')]});
subplot(235);   imagesc(x_newton,   wndImg);  	axis image off;     title({['recon_{newton}'], ['MSE : ' num2str(mse_x_newton, '%.4e')], ['PSNR : ' num2str(psnr_x_newton, '%.4f')], ['SSIM : ' num2str(ssim_x_newton, '%.4f')]});

subplot(2,3,[3,6]); semilogy(obj, '*-');    title(COST.equation);  xlabel('# of iteration');   ylabel('Objective'); 
                                            xlim([1, niter]);   grid on; grid minor;
