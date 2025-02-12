%script to generate and analyze Freidel oscillation
%cos(k0r+phi)/r^2
clear
N=20; % number of impurities
npx = 512; % number of pixels for simulation
r = 1:1:npx;
[rx,ry] = meshgrid(r,r);
%create the windowing function
x = linspace(0,pi,npx);
y = linspace(0,pi,npx);
w = sin(x)'*sin(y);
%define the center of the grid
if mod(npx, 2) == 0
 cpx = npx/2 + 1;
else
 cpx = (npx-1)/2 + 1;
end
%define k-space
kx = 1/npx*(rx-cpx);
ky = 1/npx*(ry-cpx);
k0 = [1,0]; % set k0
kappa=0;
%set what limits you want to randomize the impurity to,
%usually it is good to avoid boundaries,
%hence the floor(npx/10) term, set it as needed
minX = min(r)+floor(npx/10);
maxX = max(r)-floor(npx/10);
%set phi
phi = 0;
%get random coordinates for the defect
%within min and max limits set in
%lines 23-24
x0 = floor((maxX-minX)*rand(N,1)+minX);
y0 = floor((maxX-minX)*rand(N,1)+minX);
dnr = zeros(npx,npx);
dnqC = zeros(npx,npx);
dnrC = zeros(npx,npx);
sigma = 9;
%%%Generate signal
for n=1:N
 ys = cos(norm(k0)*sqrt((rx-x0(n)).^2+(ry-y0(n)).^2)+phi)./ ...
 (sqrt((rx-x0(n)).^2+(ry-y0(n)).^2)).^2 ...
 .*exp(-kappa*sqrt((rx-x0(n)).^2+(ry-y0(n)).^2)) + 0.05*exp(-((rx-x0(n)).^2+(ry-y0(n)).^2)/(2*sigma^2)); ...
 %set infinity values at origin to 1, that's how sinc does it
 ys(isinf(ys)) = cos(phi)*1;
 dnr = dnr + ys;
end
dnr = dnr/N; %noramlize by number of defects
 for n=1:N
 sftr = exp(1i*2*pi*(kx*(x0(n)-cpx)+ky*(y0(n)-cpx)));
 yq = fftshift(fft2(ifftshift(dnr)));
 sft_yq = yq.*sftr; %shifted in q-space
 %bring back to real space by ifft and multiple by sine window
 wind_sft_yq = ifftshift(fftshift(ifft2(ifftshift(sft_yq))).*w);
 %corrected and shifted dnr summed
 dnrC = dnrC + real(fftshift(wind_sft_yq));
 %take fft again of the windowed data
 ft_wind_sft_yq = fftshift(fft2(wind_sft_yq));
 dnqC = dnqC + ft_wind_sft_yq;
 end
dnq = fftshift(fft2(dnr));
%perform azimuthal integration
for i=cpx:-1:2
 mask = circlematrix([npx,npx],i,cpx,cpx)...
 -circlematrix([npx,npx],i-1,cpx,cpx);
 real_dnq_curve(i) = sum(sum(mask.*real(dnq)));
 imag_dnq_curve(i) = sum(sum(mask.*imag(dnq)));
 real_dnqC_curve(i) = sum(sum(mask.*real(dnqC)));
 imag_dnqC_curve(i) = sum(sum(mask.*imag(dnqC)));
end
figure;
%plot the dnr before correction
subplot(1,2,1)
imagesc(dnr);axis equal; axis off;colormap('gray');
caxis([-0.002,0.002])
% caxis([min(dnr(:)),max(dnr(:))]/100);
title(gca,'\delta N(r)')
%plot the dnr after correction, this makes sense only for N=1, for
%multiple defects, this will show sum of more than one images
subplot(1,2,2)
imagesc(dnrC);axis equal; axis off;colormap('gray');
caxis([-0.02,0.02])
% caxis([min(dnrC(:)),max(dnrC(:))]/100);
title(gca,'Shift Corrected \delta N(r)')
% define a three color colormap with 0 in the middle
%to be able to see sign
% we choose to go from Red [1,0,0] to White [1,1,1] to Blue [0,0,1]
rwb = [linspace(1,1,128),linspace(1,0,128);
 linspace(0,1,128),linspace(1,0,128);
 linspace(0,1,128),linspace(1,1,128);
 ];
rwb = rwb';

figure;
%Plot Fourier transform of Friedel oscillations
subplot(2,2,1)
imagesc(real(dnq));axis equal; axis off; colormap(rwb);
caxis([-3,3])
title(gca,'Re\delta N(q)')
subplot(2,2,2)
imagesc(imag(dnq));axis equal; axis off; colormap(rwb);
caxis([-3,3])
title(gca,'Im\delta N(q)')
%Plot Multiple Atom Fourier transform of Friedel oscillations
subplot(2,2,3)
imagesc(real(dnqC));axis equal; axis off; colormap(rwb);
caxis([-10,10])
title(gca,'Re\delta N_{MA}(q)')
subplot(2,2,4)
imagesc(imag(dnqC));axis equal; axis off; colormap(rwb);
caxis([-1,1])
title(gca,'Im\delta N_{MA}(q)')
%plot azimuthally integrated signal
figure
subplot(1,2,1)
hold on
plot([86 86],[min(real_dnq_curve),max(real_dnq_curve)],'-k','Linewidth',2)
plot(real_dnq_curve,'r','Linewidth',1)
xlabel('q (a.u.)')
ylabel('Re\delta N(q) (I_0)')
set(gca,'Xtick',[])
axis square
box on
axis tight
subplot(1,2,2)
hold on
plot([86 86],[min(real_dnqC_curve),max(real_dnqC_curve)],'-k','Linewidth',2)
plot(real_dnqC_curve,'b','Linewidth',1)
xlabel('q (a.u.)')
ylabel('Re\delta N_{MA}(q) (I_0)')
set(gca,'Xtick',[])
axis square
box on
axis tight
function cm=circlematrix(imageSize,r,x,y)
ci = [x, y, r]; % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
cm = double(uint8((xx.^2 + yy.^2)<=ci(3)^2));
end
