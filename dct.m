 clc;
 clear all;
 close all;
 warning off;
 cd test
 nm=uigetfile('*.tif','select an input image');
% % dr=dir('train');
% ln=length(dr)
% for i=3:ln
%     nm=strcat(num2str(i-2),'.tif');
 im=imread(nm); 
 cd ..
 im=imresize(im,[200 200]);
 I=double(im);
figure,
imshow(I,[]),title('input image')
%I=imresize(I,[256 256]);
if size(I,3)>1
    I=rgb2gray(I);    
end
figure,imshow(i);
% title('Original image');
I=double(I);
%----------------------------------
 Q=70;                   % quality facor
 [I J]=jcomatt(I,10);
 figure,imshow(J,[]);title(['Compressed image with DCT with quality factor of ',num2str(Q)]);
 mse=sum(sum((I-J).^2))/(size(I,1)*size(I,2))
 PSNR=20*log10(255/sqrt(mse))
%----------------------------------------
% wavelet based image  compression -----
wname='haar';
[C,S] = wavedec2(I,4,wname);
A1 = wrcoef2('a',C,S,wname,4);
H1 = wrcoef2('h',C,S,wname,4); 
V1 = wrcoef2('v',C,S,wname,4); 
D1 = wrcoef2('d',C,S,wname,4); 
% Compress the image and display it. 
% To compress the original image X, use the ddencmp command to calculate the default parameters
% and the wdencmp command to perform the actual compression. Type
[thr,sorh,keepapp] = ddencmp('cmp','wv',I);
[Xcomp,CXC,LXC,PERF0,PERFL2] = wdencmp('gbl',C,S,wname,4,thr,sorh,keepapp);
figure,imshow(Xcomp,[]);title(['Compressed image with DWT with ', wname ,' wavlet' ]);
mse=sum(sum((I-Xcomp).^2))/(size(I,1)*size(I,2))
PSNR=20*log10(255/sqrt(mse))