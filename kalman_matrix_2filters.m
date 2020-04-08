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

figure(1),imshow(preview_image);

f = waitbar(0,'','Name','Reading...');

%%Инициализация компонент сигнала и фильтра
K = 500; %количество отсчетов (картинок)

%параметр шага по пикселям
step_pixels = 3;

rows_short = round (rows /(step_pixels*2+1));
cols_short = round (columns /(step_pixels*2+1));

kk = zeros(1, rows_short*cols_short); %номера отсчетов для визуализации графиков 

%s = zeros(1, rows*columns); %сигнал
%s_f = zeros(1, rows*columns);

s = zeros(rows_short,cols_short); %сигнал
s_f = zeros(rows_short,cols_short);

Thetta = zeros(1, 3, rows_short, cols_short);
H = zeros(1, 3, rows_short, cols_short);
P = zeros(1, 3, rows_short, cols_short);

Thetta2 = zeros(1, 3, rows_short, cols_short);
H2 = zeros(1, 3, rows_short, cols_short);
P2 = zeros(1, 3, rows_short, cols_short);

for k = 1 : K
    kk(k) = k;
end

for y = 1:rows_short
    for x = 1:cols_short
        %B(y,x) = 0;
        B_est(y,x) = 0;
        %A(y,x) = 0;
        A_est(y,x) = 0;
        %ph_init(y,x) = 0;
        ph_init_est(y,x) = 0;
        Thetta(:,:,y,x) = [B_est(1); A_est(1);ph_init_est(1) ]; %начальное значение вектора параметров
    end;
end;


%выходные изображения. Здесь стоит сделать немного по-другому - массивы в
%количестве rows/(step_pixels*2+1)
output_imgs = zeros(500,cols_short,rows_short);

%ковариационная матрица шума системы (параметров)

Rpr1 = [0.04 0 0; 0 0 0; 0 0 0]; 
Rpr2 = [0 0 0; 0 0.5 0; 0 0 0.00005];
%ковариационная матрица шума наблюдения
Rn =  0.01;


%Сначала фильтр только по фону. И только здесь читаем сигнал

for k = 1 : 500 %по кол-ву изображений k=500
    %чтение пикселей 
    img_name = sprintf('%s\\%s%d%s',folder,file_mask,k,file_format);
    img = rgb2gray(imread(img_name));
    yy = 1; 
    
    %здесь читаем сигнал и оцениваем только фон
    for y = step_pixels+1:step_pixels*2+1:rows-step_pixels-1
        xx = 1;
        for x = step_pixels+1:step_pixels*2+1:columns-step_pixels-1
            s(yy,xx) = mean(mean(img(x-step_pixels:x+step_pixels,y-step_pixels:y+step_pixels))); % надо понять, как здесь сделать усреднение, чтобы начать с 1;1
            waitbar((y*columns + x)/(rows*columns),f,sprintf('%d from %d',(y*columns + x),rows*columns)); 
            
            %первый фильтр
            H(:,:,yy,xx) = [1; cos(2*pi*(y*columns + x)*0.33+Thetta(:,3,yy,xx)); -sin(2*pi*(y*columns + x)*0.33+Thetta(:,3,yy,xx))];
            P(:,:,yy,xx) = Rpr1*H(:,:,yy,xx)'*inv(H(:,:,yy,xx)*Rpr2*H(:,:,yy,xx)'+Rn);
            Thetta(:,:,yy,xx) = Thetta(:,:,yy,xx) + P(:,:,yy,xx) * (s(yy,xx) - (Thetta(:,1,yy,xx)+Thetta(:,2,yy,xx)*cos(2*pi*x*0.33+Thetta(:,3,yy,xx))));
            B_est(yy,xx) = Thetta(:,1,yy,xx); 
            
            s_f(yy,xx) = s(yy,xx) - B_est(yy,xx);
            
            % второй фильтр
            H2(:,:,yy,xx) = [1; cos(2*pi*(y*columns + x)*0.33+Thetta2(:,3,yy,xx)); -sin(2*pi*(y*columns + x)*0.33+Thetta2(:,3,yy,xx))];
            P2(:,:,yy,xx) = Rpr2*H2(:,:,yy,xx)'*inv(H2(:,:,yy,xx)*Rpr2*H2(:,:,yy,xx)'+Rn);
            Thetta2(:,:,yy,xx) = Thetta2(:,:,yy,xx) + P2(:,:,yy,xx) * (s_f(yy,xx) - (Thetta2(:,1,yy,xx)+Thetta2(:,2,yy,xx)*cos(2*pi*x*0.33+Thetta2(:,3,yy,xx))));
            
            Thetta2(:,2,yy,xx) = abs(Thetta2(:,2,yy,xx));
            A_est(yy,xx) = Thetta2(:,2,yy,xx);
            ph_init_est(yy,xx) = Thetta2(:,3,yy,xx);
            xx = xx+1;
        end
       yy = yy+1;
    end
    
    output_imgs(k,:,:) = A_est;
    disp(k);
end

%нормировка b-сканов
for i = 1:rows_short %проход по всем б-сканам
    bscan = output_imgs(:,:,i);
    for j = 1:cols_short %проход по всем столбцам текущего б-скана
        bscan_col = bscan(:,j);
        min_elem = min(bscan_col);
        max_elem = max(bscan_col);
    
        for k = 1:500 %проход по всем значениям текущего столбца текущего б-скана
            bscan_col(k) = (bscan_col(k) - min_elem)/(max_elem-min_elem)*255;
        end;
        bscan(:,j) = bscan_col;
    end;
    output_imgs(:,:,i) = bscan;
end;

%здесь идёт запись в картинки с Б-сканами
for i = 1:rows_short
   imwrite(output_imgs(:,:,i),sprintf('%d.png',i));
end

% A_est_sm = smooth(A_est,15);

% B_est_sm = smooth(B_est,25);
% 
% s_sub = s - B_est;
% 
% A_est_sm = A_est_sm';
% B_est_sm = B_est_sm';


% figure (3),
% plot(kk,B_est,kk, s),title('Фон');
% plot(kk,A_est, kk, s),title('Амплитуда');
% plot(kk,ph_init_est),title('Начальная фаза');