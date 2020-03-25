clc
close all;
clear all;

%��������� ��� ������ �����������
prompt = {'File mask:','File format:','N:'};
title = 'Input data';
dims = [1 40];
definput = {'img','.bmp','500'};
answer = inputdlg(prompt,title,dims,definput);

file_mask = answer{1};
file_format = answer{2};

N = str2num(answer{3}); %���������� ������

%�������� ����� � ������� � �����
folder = uigetdir('A:\Mosquito');

preview_image = imread(sprintf('%s\\%s1%s',folder,file_mask,file_format));

[rows, columns, channels_num] = size(preview_image); %����������� ������� �����������
%rows = 100;
%columns = 100;

figure(1),imshow(preview_image);

f = waitbar(0,'','Name','Reading...');

%%������������� ��������� ������� � �������
K = 500; %���������� �������� (��������)

kk = zeros(1, rows*columns); %������ �������� ��� ������������ �������� 

s = zeros(1, rows*columns); %������
s_f = zeros(1, rows*columns);

Thetta = zeros(1, 3, rows, columns);
H = zeros(1, 3, rows, columns);
P = zeros(1, 3, rows, columns);

for y = 1:rows
    for x = 1:columns
        B(y,x) = 0;
        B_est(y,x) = 0;
        A(y,x) = 0;
        A_est(y,x) = 0;
        ph_init(y,x) = 0;
        ph_init_est(y,x) = 0;
        Thetta(:,:,y,x) = [B_est(1); A_est(1);ph_init_est(1) ]; %��������� �������� ������� ����������
    end;
end;

%�������� ���� �� ��������
step_pixels = 5;
%�������� ����������� 
output_imgs = zeros(500,200,round(rows/step_pixels));

%�������������� ������� ���� ������� (����������)
Rpr = [0.073 0 0; 0 0.021 0; 0 0 0.006];
%�������������� ������� ���� ����������
Rn =  0.01;

for k = 1 : 500 %�� ���-�� �����������
    %������ ��������
    img_name = sprintf('%s\\%s%d%s',folder,file_mask,k,file_format);
    img = rgb2gray(imread(img_name));
    for y = 1:rows 
        for x = 1:columns
        s(x) = mean(mean(img(x:x+step_pixels:x-step_pixels,y:y+step_pixels:y-step_pixels))); % ���� ������, ��� ����� ������� ����������, ����� ������ � 1;1
        waitbar((y*columns + x)/(rows*columns),f,sprintf('%d from %d',(y*columns + x),rows*columns)); 
        %����������
        H(:,:,y,x) = [1; cos(2*pi*(y*columns + x)*0.33+Thetta(:,3,y,x)); -sin(2*pi*(y*columns + x)*0.33+Thetta(:,3,y,x))];
        P(:,:,y,x) = Rpr*H(:,:,y,x)'*inv(H(:,:,y,x)*Rpr*H(:,:,y,x)'+Rn);
        Thetta(:,:,y,x) = Thetta(:,:,y,x) + P(:,:,y,x) * (s(x) - (Thetta(:,1,y,x)+Thetta(:,2,y,x)*cos(2*pi*x*0.33+Thetta(:,3,y,x))));
        Thetta(:,2,y,x) = abs(Thetta(:,2,y,x));
        B_est(y,x) = Thetta(:,1,y,x);
        A_est(y,x) = Thetta(:,2,y,x);
        ph_init_est(y,x) = Thetta(:,3,y,x);
        
        if (y == columns) %����� �� ����� �����������
            %�������� � ���� ������ ������
            output_imgs(:,k, y/5) = A_est; %������, ����� ������ ���� ������
            ekf_data = fopen(sprintf('%s%d%s','ekf_data',k,'.txt'),'wt');
            fprintf(ekf_data,'%f \n', A_est);
            fclose(ekf_data);
            %�������� ��� ��������� ��� ���������� �����������
            A_est = zeros(1,rows*columns);
            B_est = zeros(1,rows*columns);
            ph_init_est = zeros(1,rows*columns);  
        end;
       end
    end
end
A_est_sm = smooth(A_est,15);
B_est_sm = smooth(B_est,25);

s_sub = s - B_est;

A_est_sm = A_est_sm';
B_est_sm = B_est_sm';


figure (3),
plot(kk,B_est_sm,kk, s_sub),title('���');
plot(kk,A_est_sm, kk, s_sub),title('���������');
plot(kk,ph_init_est),title('��������� ����');