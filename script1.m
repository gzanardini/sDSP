%close all;
clear
%% load images 
stored = load('img_restoration.mat');

imageArray = stored.I1; %%%%% CHANGE NAME HERE %%%%%

figure
%subplot(1,2,1);
imshow(imageArray, []);
title('Original Image')

%% parameters for blurring

T = 1;
a = 0.01;
b = 0.01;

%%%%%%%%%%%% Image blurring %%%%%%%%%%%%

[M, N] = size(imageArray);

for u = 1:M
    for v = 1:N
        H(u,v) = T./(pi*(u*a + v*b)).*sin(pi*(u*a + v*b)).*exp(-1j*pi*(u*a + v*b));
    end
end
Im_fft = fft2(imageArray);

blurred_fft =Im_fft.*H;
blurred_ifft = abs(ifft2(blurred_fft));
%figure('Name','Image - Motion Blurred ');
% subplot(2,2,3);
% imshow(blurred_ifft,[])
% title('Image - Motion Blurred')


%% inverse filtering -  without noise

inversed = blurred_fft ./ H;
inversed_ifft = abs(ifft2(inversed));
% figure,subplot(1,2,1),imshow(inversed_ifft,[]),title('Inverse Filtered without noise')

%% adding gaussian noise to image with AWGN function 

im_noise_added_without_motion_awgn = awgn(Im_fft,10,'measured');

awgn_noisy_ifft = abs(ifft2(im_noise_added_without_motion_awgn));
% figure,subplot(1,2,1),imshow(awgn_noisy_ifft,[]),title('AWGN Noisy Image 1')

%% adding gaussian noise to motion blurred images

im_awgn_noise_added_blurred = awgn(blurred_fft,10,'measured');

blurred_awgn_noisy_ifft = abs(ifft2(im_awgn_noise_added_blurred));

%figure,subplot(2,2,1),imshow(blurred_awgn_noisy_ifft,[]),title('Motion Blurred and AWGN Noise added')

%% inverse filtering -  with gaussian noise (AWGN)

noisy_inversed = im_awgn_noise_added_blurred ./ H;
noisy_inversed_ifft = abs(ifft2(noisy_inversed));
% subplot(2,2,3),imshow(noisy_inversed_ifft,[]),title('Just Inversed')


%% adding random noise
std = 10;
mean = 0;

im_rand_noise = imageArray+std*randn(M,N)+mean;
% figure,imshow(im_rand_noise,[]),title('Random Noise added - NO BLUR')

im_blur_rand_noise = real(blurred_ifft)+std*randn(M,N);
% figure,imshow(im_blur_rand_noise,[]),title('Random Noise added - WITH BLUR')

%% Wiener Filter with RANDOM NOISE

estimated_nsr = std / sqrt(var(imageArray(:)));

rand_blurred_ifft = blurred_ifft+std*randn(M,N)+mean;

rand_blurred_fft = fft2(rand_blurred_ifft);

%f_rand_wiener = (((abs(H)).^2))./(H*(estimated_nsr+(abs(H)).^2)).*rand_blurred_fft;
f_rand_wiener = (1./H).*(((abs(H)).^2))./((estimated_nsr+(abs(H)).^2)).*rand_blurred_fft;
rand_wiened = abs(ifft2(f_rand_wiener));
% figure
% subplot(1,2,1);
% imshow(rand_blurred_ifft, []);
% title('Motion Blurred and rand Noise added')
% subplot(1,2,2);
% imshow(rand_wiened, [])
% %title(["rand wiened - K =",num2str(K)])
% title('rand wiened')

%% comparison of our wiener filter with matlab wiener filter

%K = power_spec_noise./power_spec;

figure
subplot(1,3,1)
imshow(imageArray, [])
title('Original Image')

subplot(1,3,2)
imshow(rand_wiened, [])
title('Restoration with our Wiener Filter')

LEN = 25;
THETA = 45;
PSF = fspecial('motion', LEN, THETA);
blurred = imfilter(imageArray, PSF, 'conv', 'circular');
% figure, imshow(blurred,[])
% title('Simulate Blur')


im_blur_rand_noise = blurred+std*randn(M,N)+mean;
% figure
% imshow(im1_rand_noise,[])
% title('Random Noise added')

wnr3 = deconvwnr(im_blur_rand_noise, PSF, estimated_nsr);
%figure
subplot(1,3,3)
imshow(wnr3,[])
title('Restoration of with MATLAB function');


%% plotting stuff

figure,subplot(1,3,1),imshow(imageArray, []),title('Original Image')

subplot(1,3,2),imshow(im_blur_rand_noise,[]),title('Blurred, Random Noisy Image')

subplot(1,3,3),imshow(noisy_inversed_ifft,[]),title('Inverse Filter Applied')
%
figure,subplot(1,3,1),imshow(imageArray, []),title('Original Image')

subplot(1,3,2),imshow(blurred_ifft,[]),title('Original Image - Motion Blurred')

subplot(1,3,3),imshow(inversed_ifft,[]),title('Inverse Filtered without noise')





%% not relevant
% %% Wiener AWGN
% K = 0.1;
% 
%
% %f_wiener = (1./H) .* (abs(H)^2 ./ (abs(H)^2 + K)) .* im_noise_added;
% f_wiener = (conj(H) ./ (abs(H)^2 + K)) .* im_awgn_noise_added_blurred;
% wiened = abs(ifft2(f_wiener));
% % figure,subplot(1,2,1),imshow(blurred_awgn_noisy_ifft, []),title('Motion Blurred and AWGN Noise added')
% % subplot(1,2,2),imshow(wiened, []),title(["wiened - K = ",num2str(K)])