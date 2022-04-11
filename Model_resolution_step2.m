clear all;
% close all;

%% Parameters
global ph_size grid_size r n_transducers c wn f_start f_end fh_start fh_end wnh snr fov;
ph_size = 100;
grid_size = 500;
r = 55;
n_transducers = 100;
c = 1500;
wn = 1000;
f_start = 100000;
f_end = 50000000;

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

%% Convert weight matrix to WI
disp("Calculating WI");
% Calculating WI
WI = zeros(wn*n_transducers, ph_size*ph_size);
for i=1:wn
%     temp = reshape(weight_matrix(i,:), [], n_transducers)';
    temp = weight_matrix(i,:);
    temp = reshape(temp, [], n_transducers);
    temp = temp';
    WI((i-1)*n_transducers+1:i*n_transducers,:) = temp;
end
WI = [real(WI); imag(WI)];

clear temp dist bottom c f_end f_start grid grid_size i j k;
clear left n_transducers offset P r right theta top w_range weight_matrix wn;


%% Calculate A
disp("Calculating A");

reg_param = 5e11;
A1 = WI' * WI;
A = A1 + reg_param* speye(ph_size*ph_size);

%% Model resolution part
load("Images/good_results/" + num2str(snr) + "_shepp_mat_100.mat");

% % Recon
% disp("Recon");
% disp("Calculating Blur");
% Blur=A\((WI')*(WI));
% [M, N] = size(Blur);
% 
% disp("Calculating F");
% F1 = (Blur')*(Blur); 
% 
% lambda=0.01;
% alpha = lambda / 0.01;
% 
% F = F1 + (alpha)*speye(ph_size*ph_size);
% issparse(F);
% 
% X0new3 = reshape(reconimage', [], 1);
% modrescac = (Blur')*( X0new3);
% t = zeros(size(modrescac));
% 
% Nit=100;
%  for i = 1:Nit
%      i
%      u = soft(modrescac+t,0.5*(lambda/alpha))-t;
%      modrescac= F\((Blur')*(X0new3) + alpha*u);
%      t = modrescac-u;
%      
%  end
%  
%  modresoutput_quad = reshape(modrescac, ph_size, [])';
%  imwrite(modresoutput_quad, "Images/good_results/" + num2str(snr) + "_shepp_modres_recon_0025.png");
%  clear Blur X0new3 modrescac t F u;


% Quadratic
disp("Quadratic");
disp("Calculating Blur");
Blur=A\((WI')*(WI));

[M, N] = size(Blur);

disp("Calculating F");
F1 = (Blur')*(Blur); 

lambda=0.1;
alpha = lambda / 0.01;

F = F1 + (alpha)*speye(ph_size*ph_size);
issparse(F);

X0new3 = reshape(reconimagequadratic', [], 1);
modrescac = (Blur')*( X0new3);
modrescac2 = (Blur')*( X0new3);
t = zeros(size(modrescac));

Nit=100;
 for i = 1:Nit
     i
     u = soft(modrescac+t,0.5*(lambda/alpha))-t;
     modrescac= F\(modrescac2 + alpha*u);
     t = modrescac-u;
     
 end
 
 modresoutput_quad = reshape(modrescac, ph_size, [])';
 imwrite((modresoutput_quad), "Images/good_results/" + num2str(snr) + "_shepp_modres_quad_0025.png");
%  clear Blur X0new3 modrescac t F u;
 
 
 % ABS
 disp("Absolute");
disp("Calculating Blur");
% Blur=A\((WI')*(WI));
% [M, N] = size(Blur);
% 
% disp("Calculating F");
% F1 = (Blur')*(Blur); 

% lambda=0.1;
% alpha = lambda / 0.01;

F = F1 + (alpha)*speye(ph_size*ph_size);
issparse(F);

X0new3 = reshape(reconimagenew1', [], 1);
modrescac = (Blur')*( X0new3);
t = zeros(size(modrescac));

% Nit=100;
 for i = 1:Nit
     i
     u = soft(modrescac+t,0.5*(lambda/alpha))-t;
     modrescac= F\((Blur')*(X0new3) + alpha*u);
     t = modrescac-u;
     
 end
 
 modresoutput_abs = reshape(modrescac, ph_size, [])';
 imwrite((modresoutput_abs), "Images/good_results/" + num2str(snr) + "_shepp_modres_abs_0025.png");
 % clear Blur X0new3 modrescac t F u;
 
 % cauchy
 disp("Cauchy");
% disp("Calculating Blur");
% Blur=A\((WI')*(WI));
% [M, N] = size(Blur);
% 
% disp("Calculating F");
% F1 = (Blur')*(Blur); 

% lambda=0.1;
% alpha = lambda / 0.01;

F = F1 + (alpha)*speye(ph_size*ph_size);
issparse(F);

X0new3 = reshape(reconimagenew2', [], 1);
modrescac = (Blur')*( X0new3);
t = zeros(size(modrescac));

% Nit=100;
 for i = 1:Nit
     i
     u = soft(modrescac+t,0.5*(lambda/alpha))-t;
     modrescac= F\((Blur')*(X0new3) + alpha*u);
     t = modrescac-u;
     
 end
 
 modresoutput_cauchy = reshape(modrescac, ph_size, [])';
%  figure;subplot(2,1,1);imshow(reconimagenew2);subplot(2,1,2);imshow(modresoutput_cauchy);
 
 imwrite((modresoutput_cauchy), "Images/good_results/" + num2str(snr) + "_shepp_modres_cauchy_0025.png");
 % clear Blur X0new3 modrescac t F u;
 
 % Geman
 disp("Geman");
disp("Calculating Blur");
% Blur=A\((WI')*(WI));
% [M, N] = size(Blur);
% 
% disp("Calculating F");
% F1 = (Blur')*(Blur); 

% lambda=0.1;
% alpha = lambda / 0.01;

F = F1 + (alpha)*speye(ph_size*ph_size);
issparse(F);

X0new3 = reshape(reconimagenew3', [], 1);
modrescac = (Blur')*( X0new3);
t = zeros(size(modrescac));



% Nit=100;
 for i = 1:Nit
     disp(i);
     u = soft(modrescac+t,0.5*(lambda/alpha))-t;
     modrescac= F\((Blur')*(X0new3) + alpha*u);
     t = modrescac-u;
     
 end
 
 modresoutput_geman = reshape(modrescac, ph_size, [])';
%  imshow(modresoutput_geman);
%  figure;subplot(2,1,1);imshow(reconimagenew3);subplot(2,1,2);imshow(modresoutput_geman);
 
 imwrite((modresoutput_geman), "Images/good_results/" + num2str(snr) + "_shepp_modres_geman_0025.png");
%  clear Blur X0new3 modrescac t F;
 
save("Images/good_results/" + num2str(snr) + "_shepp_modres_all_003.mat", "modresoutput_quad", "modresoutput_abs", "modresoutput_cauchy", "modresoutput_geman");
figure;subplot(2,2,1);imshow(modresoutput_quad,[]);subplot(2,2,2);imshow(modresoutput_abs,[]);subplot(2,2,3);imshow(modresoutput_cauchy,[]);subplot(2,2,4);imshow(modresoutput_geman,[]);
