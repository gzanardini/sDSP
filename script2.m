close all;
clear;

%% load images 
stored = load('img_restoration.mat');
imageArray = stored.I2;                    %%%%%%%%%%%%%%% change image here %%%%%%%%%%%%%%%
imageArray = uint8(imageArray);
imageArray = im2double(imageArray);

figure
imshow(imageArray, []);
title('Original Image')

%% parameters

T = 1;
a = 0.01;
b = 0.01;

%%%%%%%%%%%% I blurring %%%%%%%%%%%%
[M, N] = size(imageArray);


for u = 1:M
    for v = 1:N
        H(u,v) = T./(pi*(u*a + v*b)).*sin(pi*(u*a + v*b)).*exp(-1j*pi*(u*a + v*b));
    end
end
Im_fft = fft2(imageArray);
blurred_fft =Im_fft.*H;
blurred_ifft = abs(ifft2(blurred_fft));
%figure,imshow(blurred_ifft,[]),title('Motion blurred image')

%% adding gaussian noise
std = 10;
mean = 0;

im_gaus_noise = imnoise(imageArray, 'gaussian',0,0.003);
%figure,imshow(im_gaus_noise,[]),title('gaussian Noise added - NO BLUR')

im_blur_gaus_noise = imnoise(real(blurred_ifft),'gaussian',0,0.003);
%figure,imshow(im_blur_gaus_noise,[]),title('gaussian Noise added - WITH BLUR')

noise=real(blurred_ifft)-im_blur_gaus_noise;
 figure
 imshow(noise,[])
 title('noise')

% % % power spectrum calculation

power_spec = (abs(Im_fft)).^2;
noise_fft = fft2(noise)-Im_fft;
power_spec_noise = (abs(noise_fft)).^2;

estimated_nsr = power_spec_noise./power_spec;

%% Wiener Filter with Gaussian NOISE

%estimated_nsr = 0.2; %std^2 / var(imageArray(:));


gaus_blurred_fft = fft2(im_blur_gaus_noise);
f_gaus_wiener = (1./H).*(((abs(H)).^2))./((estimated_nsr+(abs(H)).^2)).*gaus_blurred_fft;
gaus_wiened = abs(ifft2(f_gaus_wiener));

%figure,imshow(gaus_wiened, []),title('gaus wiened')


%% comparison of our wiener filter with matlab wiener filter
figure
subplot(1,3,1)
imshow(imageArray, [])
title('Original Image')

subplot(1,3,2)
imshow(gaus_wiened, [])
title('Restoration with our Wiener Filter')

LEN = 25;
THETA = 45;
PSF = fspecial('motion', LEN, THETA);
blurred = imfilter(imageArray, PSF, 'conv', 'circular');
% figure, imshow(blurred,[])
% title('Simulate Blur')

matlab_blurred_noisy = imnoise(blurred,'gaussian',0,0.003);

wnr3 = deconvwnr(matlab_blurred_noisy, PSF, estimated_nsr);
%figure
subplot(1,3,3)
imshow(wnr3,[])
title('Restoration with MATLAB Wiener Filter');

%% extra filters and restorations 

%% arithmetic mean 

window = 5; % nr of windows
for i = 1:M-window+1
    for j = 1:N-window+1
        g = im_blur_gaus_noise(i:i+window-1,j:j+window-1);
        sumsum = sum(sum(g));
        im_meaned(i+(window-1)/2,j+(window-1)/2) = sumsum/(window*window);
    end
end

figure
subplot(1,3,1),imshow(imageArray,[]);
title('Original Image')
subplot(1,3,2),imshow(im_blur_gaus_noise,[]);
title('Noisy and blurry image')
subplot(1,3,3),imshow(im_meaned,[]);
title('Restored Image with the Arithmetic Mean Filter')


%% median filter

medianed_noisy_blurry = medfilt2(im_blur_gaus_noise,[5 5]);
figure
subplot(1,3,1)
imshow(imageArray)
title('Original Image')
subplot(1,3,2)
imshow(im_blur_gaus_noise,[]);
title('Noisy and blurry image')
subplot(1,3,3)
imshow(medianed_noisy_blurry,[])
title('Restored Image with the Median Filter')

%% still playing around with median filter


estimated_nsr = 0.005; %playing around

L_fft = fft2(medianed_noisy_blurry);
f_L_wiener = (1./H).*(((abs(H)).^2))./((estimated_nsr+(abs(H)).^2)).*L_fft;
L_wiened = abs(ifft2(f_L_wiener));

figure
subplot(1,3,1)
imshow(imageArray)
title('Original Image')


medianed_wiened = medfilt2(gaus_wiened,[5 5]);subplot(1,3,2)
imshow(medianed_wiened)
title('wiener then median')


subplot(1,3,3)
imshow(histeq(medianed_wiened))
title('wiener then median then histeq')