clear all;
close all;
%% Parameters
global ph_size grid_size r n_transducers c wn f_start f_end fh_start fh_end wnh snr fov;
ph_size = 64;
grid_size = 500;
r = 55;
n_transducers = 64;
c = 1500;
wn = 64;
f_start = 100000;
f_end = 2500000;

fh_end = 600000;
fh_start = 250000;
wnh = 0;

snr = 40;
fov = 300;

%% Calculate W matrix
grid = zeros(grid_size, grid_size);
left = floor((grid_size - ph_size) / 2 );
right = floor((grid_size + ph_size) / 2 - 1);
top = floor((grid_size - ph_size) / 2 );
bottom = floor((grid_size + ph_size) / 2 - 1);

grid(top:bottom, left:right) = phantom(ph_size);

theta = linspace(0*fov*2*pi/360, 1*fov*2*pi/360, n_transducers+1);
theta = theta(1:n_transducers);
offset = floor(grid_size / 2);
x1 = r .* cos(theta) + offset;
y1 = r .* sin(theta) + offset;
% 
x1 = [x1(2:end) x1(1)];
x1 = fliplr(x1);
y1 = [y1(2:end) y1(1)];
y1 = fliplr(y1);

x2 = left:right;
y2 = top:bottom;

count = 0;
dist = zeros(1,n_transducers * ph_size*ph_size);
for i=1:length(x1)
   for j=1:length(x2)
       for k=1:length(y2)
           count = count+1;
           temp = sqrt((x1(i) - x2(k))^2 + (y1(i) - y2(j))^2);
           dist(count) = temp * 0.01;
       end
   end
end
reg_param = 2e11;

dist = reshape(dist, n_transducers*ph_size*ph_size, 1);
% w_range = linspace(f_start, f_end, wn);
w_range = [linspace(f_start, f_end, wn) linspace(fh_start, fh_end, wnh)];
wn = wn + wnh;
weight_matrix = zeros(wn, n_transducers * ph_size*ph_size);

disp("Calculating Weight matrix");
for i=1:wn
    disp(i);
    weight_matrix(i,:) = Calculate_W_matrix(w_range(i), dist);
end
% weight_matrix = normalize_img(weight_matrix);
% Purge unnecessary variables
clear temp x1 x2 y1 y2;

%% Calculate P = WX
P = zeros(wn, n_transducers);
% X = reshape(phantom(ph_size), ph_size*ph_size, 1);

% 
% X = imread("Derenzo.png");
% X = imresize(X, [ph_size,ph_size]);
% X = double(im2bw(X));
% gt = X;

% load("Phantoms/07.mat");
% % sliceViewer(volimage);
% X = squeeze(volimage(166,:,:));
% X = imresize(double(X), [ph_size,ph_size]);
% X = X/5;
% gt = X;

X = imread("binary_3.png");
X = imresize(X, [ph_size,ph_size]);
X = im2bw(X);
gt = X;

% X = phantom(ph_size);
% X = flipud(X);
X = double(reshape(X', [], 1));

disp("Calculating Pressure matrix");
for i=1:wn
%     temp = reshape(weight_matrix(i,:), [], n_transducers)';
    temp = weight_matrix(i,:);
    temp = reshape(temp, [], n_transducers);
    temp = temp';
    P(i,:) = temp * X; 
end
figure;
imshow(P,[])


% Purge unnecessary variables
clear temp X;

%% Calculate small WI
% ph_size = 100;
% disp("Calculating Small Weight matrix");
% for i=1:wn
%     disp(i);
%     weight_matrix(i,:) = Calculate_W_matrix(w_range(i), dist);
% end

% Arrange P, W by concatenating Real and Imaginary parts

% % Calculating PI
% disp("Calculating PI");
% 
% PI = [reshape(real(P)', [], 1); reshape(imag(P)', [], 1)];
% PI = addNoise(PI, snr, 'peak');
% 
% disp("Calculating WI");
% % Calculating WI
% WI = zeros(wn*n_transducers, ph_size*ph_size);
% for i=1:wn
% %     temp = reshape(weight_matrix(i,:), [], n_transducers)';
%     temp = weight_matrix(i,:);
%     temp = reshape(temp, [], n_transducers);
%     temp = temp';
%     WI((i-1)*n_transducers+1:i*n_transducers,:) = temp;
% end
% WI = [real(WI); imag(WI)];




% Purge unnecessary variables
% clear temp dist bottom c f_end f_start grid grid_size i j k;
clear left offset r right theta top w_range weight_matrix wn;

%% Inversion Model

disp("Calculating B");
B = WI' * PI;

disp("Calculating A");

A1 = WI' * WI;
A = A1 + reg_param* speye(ph_size*ph_size);

disp("Calculating X0");
X0 = A\B;
% X0 = (X0 - min(X0)) / (max(X0) - min(X0));

reconimage = reshape(X0, ph_size, [])';

figure;
subplot(2,1,1);imshow(reconimage, []);subplot(2,1,2);imshow(gt,[]);

% max(max(WI))
% min(min(WI))
% max(max(dist))
% min(min(dist))
% 
% imwrite(reconimage, "Images/good_results/" + num2str(snr) + "_shepp_recon_100.png");

%% Penalization
% 
v = var(reconimage(:));

disp("Quadratic Penalization");
quadratic=1/v;
Aquadratic=(WI')*(WI)+(reg_param)*quadratic*speye(ph_size*ph_size);
X0quadratic=Aquadratic\B;
reconimagequadratic = reshape(X0quadratic, ph_size, [])';
%    figure;
%    imshow(reconimagequadratic);
   
clear quadratic Aquadratic X0quadratic;

% v = var(reconimagequadratic(:));
disp("Absolute Penalization");
absolutevalue=1./(sqrt(v)*abs(X0));
absolutevaluer=absolutevalue.*speye(ph_size*ph_size);
Anew1=(WI')*(WI)+(reg_param)*absolutevaluer;
X0new1=Anew1\B;
reconimagenew1=reshape(X0new1,ph_size, [])';
%      figure;
%    imshow(reconimagenew1);
   
clear absolutevalue absolutevaluer Anew1 X0new1;

disp("Cauchy Penalization");
cauchy=1./(v+((X0).^2));
cauchyr=cauchy.*speye(ph_size*ph_size);
Acauchy=(WI')*(WI)+(reg_param)*cauchyr;
X0new2=Acauchy\B;
reconimagenew2=reshape(X0new2,ph_size, [])';
%      figure;
%    imshow(reconimagenew2);
   
clear cauchy cauchyr Acauchy X0new2;

disp("Geman Penalization");
Geman=1./((v+((X0).^2)).^2);
Gemanr=Geman.*speye(ph_size*ph_size);
Ageman=(WI')*(WI)+(reg_param)*Gemanr;
X0new3=Ageman\B;
reconimagenew3=reshape(X0new3,ph_size,ph_size)';
%   figure;
%    imshow(reconimagenew3);
   
clear Geman Gemanr Ageman X0new3;

imwrite(gt, "Images/good_results/" + num2str(snr) + "_shepp_gt_100.png");
imwrite((reconimage), "Images/good_results/" + num2str(snr) + "_shepp_recon_100.png");
imwrite((reconimagequadratic), "Images/good_results/"  + num2str(snr) + "_shepp_quad_100.png");
imwrite((reconimagenew1), "Images/good_results/"  + num2str(snr) + "_shepp_abs_100.png");
imwrite((reconimagenew2), "Images/good_results/"  + num2str(snr) + "_shepp_cauchy_100.png");
imwrite((reconimagenew3), "Images/good_results/"  + num2str(snr) + "_shepp_geman_100.png");

save("Images/good_results/" + num2str(snr) + "_shepp_mat_100.mat", "reconimage", "reconimagequadratic", "reconimagenew1", "reconimagenew2", "reconimagenew3");

figure;subplot(2,2,1);imshow(reconimagequadratic,[]);subplot(2,2,2);imshow(reconimagenew1,[]);subplot(2,2,3);imshow(reconimagenew2,[]);subplot(2,2,4);imshow(reconimagenew3,[]);
psnr(gt, reconimagequadratic)
psnr(gt, reconimagenew1)
psnr(gt, reconimagenew2)
psnr(gt, reconimagenew3)

% %% Model Resolution : Go to shepp_modres.m
%  
% 
% 
