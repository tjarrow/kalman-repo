clc
close all;
clear all;

%настройки для чтения изображений
prompt = {'File mask:','File format:','N:'};
title = 'Input data';
dims = [1 40];
definput = {'img','.bmp','500'};
answer = inputdlg(prompt,title,dims,definput);

file_mask = answer{1};
file_format = answer{2};

N = str2num(answer{3}); %Количество файлов

%Выбираем папку с данными и точку
folder = uigetdir('C:\Users\userg\Documents\inst\Mosquito'); 
%folder = uigetdir('A:\Mosquito'); 

preview_image = imread(sprintf('%s\\%s1%s',folder,file_mask,file_format));

%[rows, columns, channels_num] = size(preview_image); %определение размера изображения
rows = 1000; %чтобы можно было смотреть в воркспейсе значения
columns = 1000;

%figure(1),imshow(preview_image);

%f = waitbar(0,'','Name','Reading...');

%%Инициализация компонент сигнала и фильтра
K = 500; %количество отсчетов (картинок)

%параметр шага по пикселям
step_pixels = 3;

rows_short = round (rows /(step_pixels*2+1));
cols_short = round (columns /(step_pixels*2+1));

kk = zeros(1, rows_short*cols_short); %номера отсчетов для визуализации графиков 

%s = zeros(1, rows*columns); %сигнал
%s_f = zeros(1, rows*columns);

s = zeros(1,500); %сигнал
s_f = zeros(1,500);

Thetta = zeros(500, 500); %нужно будет здесь поменять на 4
H = zeros(500, 500);
P = zeros(500, 500);

Thetta2 = zeros(1, 3, 500, 500);
H2 = zeros(1, 3, 500, 500);
P2 = zeros(1, 3, 500, 500);

% for y = 1:500
%     for x = 1:500
%         %B(y,x) = 0;
%         B_est(y,x) = 0;
%         %A(y,x) = 0;
%         A_est(y,x) = 0;
%         %ph_init(y,x) = 0;
%         ph_init_est(y,x) = 0;
%         Thetta(:,:,y,x) = [B_est(1); A_est(1);ph_init_est(1) ]; %начальное значение вектора параметров
%     end;
% end;
A_est_sc = zeros(1,500); 
B_est_sc = zeros(1,500);
ph_est_sc = zeros(1,500);
%выходные изображения. Здесь стоит сделать немного по-другому - массивы в
%количестве rows/(step_pixels*2+1)
%output_imgs = zeros(500,cols_short,rows_short);

%ковариационная матрица шума системы (параметров)
Rpr1 = 0.00001;
Rpr2 = [0 0 0; 0 0.003 0; 0 0 0.000005];
%ковариационная матрица шума наблюдения
Rn1 =  0.0001;
Rn2 =  0.001;
x = rows/2; y = columns / 2;
yy = rows/2; 
xx = columns / 2;

 img_name = sprintf('%s\\%s%d%s',folder,file_mask,1,file_format);
 img = rgb2gray(imread(img_name));
 s(1) = mean(mean(img(x-step_pixels:x+step_pixels,y-step_pixels:y+step_pixels)));
 Thetta(yy,xx) = s(1); %будет в цикле по (x,y) в матричном

for k = 1 : 500 %по кол-ву изображений k=500
   img_name = sprintf('%s\\%s%d%s',folder,file_mask,k,file_format);
   img = rgb2gray(imread(img_name));
   s(k) = mean(mean(img(x-step_pixels:x+step_pixels,y-step_pixels:y+step_pixels)));
        %waitbar((y*columns + x)/(rows*columns),f,sprintf('%d from %d',(y*columns + x),rows*columns)); 
        
    %первый фильтр
     %H(:,:,yy,xx) = [1; cos(2*pi*k*0.33+Thetta(:,3,yy,xx)); -sin(2*pi*k*0.33+Thetta(:,3,yy,xx))];
     H(yy,xx) = 1;
     P(yy,xx) = Rpr1*H(yy,xx)'*inv(H(yy,xx)*Rpr1*H(yy,xx)'+Rn1);
     Thetta(yy,xx) = Thetta(yy,xx) + P(yy,xx) * (s(k) - Thetta(yy,xx));
     B_est_sc(k) = Thetta(yy,xx); 
          
     s_f(k) = s(k) - B_est_sc(k);
            
    % второй фильтр
    H2(:,:,yy,xx) = [1; cos(2*pi*k*0.33+Thetta2(:,3,yy,xx)); -sin(2*pi*k*0.33+Thetta2(:,3,yy,xx))];
    P2(:,:,yy,xx) = Rpr2*H2(:,:,yy,xx)'*inv(H2(:,:,yy,xx)*Rpr2*H2(:,:,yy,xx)'+Rn2);
    Thetta2(:,:,yy,xx) = Thetta2(:,:,yy,xx) + P2(:,:,yy,xx) * (s_f(k) - (Thetta2(:,1,yy,xx)+Thetta2(:,2,yy,xx)*cos(2*pi*x*0.33+Thetta2(:,3,yy,xx))));
            
    Thetta2(:,2,yy,xx) = abs(Thetta2(:,2,yy,xx));
    A_est_sc(k) = Thetta2(:,2,yy,xx);
    ph_est_sc(k) = Thetta2(:,3,yy,xx);
   
    idx(k) = k;
    disp(k);
end

% i = zeros(1,500);
% a = A_est_sc(:,1)';
% b = B_est_sc(:,1)';
% phi = ph_est_sc(:,1)';

%sig = signal_if - b;
A = smooth(A_est_sc, 15);
figure (1);
plot(idx,B_est_sc,idx,s);
figure (2);
plot(idx,A,idx,s_f);
figure(3);
plot(idx,ph_est_sc);