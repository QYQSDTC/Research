clc;
clear;
A=imread('little_fairy.jpeg');
if ndims(A)==3
   A=rgb2gray(A);% 3D to 2D
end
I=fft2(double(A));% change data to type double
B=fftshift(I);% move the zero zone to the center
C=log(1+abs(B));% show beatifully
D = real(ifft2(I));
subplot(1,3,1),imshow(A),title('oringinal image');
subplot(1,3,2),imshow(C,[]),title('image of FFT');
subplot(1,3,3),imshow(D,[]),title('image of inverse FFT');
%  figure
%  imshow(C,[])