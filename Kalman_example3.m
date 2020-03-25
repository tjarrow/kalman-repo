clc
close all;
clear all;

%%
%���������� ������
K = 15562; %���������� ��������
kk = zeros(1, K); %������ �������� ��� ������������ �������� 
B = zeros(1, K); %���
A = zeros(1, K); %���������
ph_init = zeros(1,K); %��������� ����
for k = 1 : K
    %B(k) = 0.2*sin(2*pi*k*0.002);
    %A(k) = exp(-((k-400)^2)/100^2)+0.6*exp(-((k-750)^2)/100^2);
    %s(k) = B(k) + A(k)*cos(2*pi*k*0.05 + ph_init) + n(k);
    kk(k) = k;
end
%figure(1), plot(s);
%plot(B);
T = 5;

%�������� ������
file_signal = fopen('signal_all.txt', 'r');
formatSpec = '%f';
signal = fscanf(file_signal, formatSpec);

% signal_img = signal';
% signal_f = fft(signal_img);
% signal_f2 = fft(signal_img);
% fft_idx = 100;
% %fft_idx = 150; - ��������� ����� ������ 
% fft_idx2 = 420;
% for i = 1:fft_idx;
%     signal_f(i) = 0;
% end;
% for i = 1:fft_idx2
%     signal_f2(i) = 0;
% end;
% signal_if = ifft(signal_f);
% signal_if2 = ifft(signal_f2);
% signal_sub2 = signal_if; % - signal_if2;

%%
% �������� ����������� ��������� �������
B_est = zeros(1, K); %������ ����
A_est = zeros(1, K); %������ ���������
ph_init_est = zeros(1, K); %������ ����
T_est = zeros(1, K); %������ �������


%��������� ��������
B_est(1) = 0;
A_est(1) = 0;
ph_init_est(1) = 0;
T_est(1) = T + 0.1;

Thetta = [B_est(1); A_est(1);ph_init_est(1); T_est(1) ]; %��������� �������� ������� ����������

%�������������� ������� ���� ������� (����������)
Rpr = [0.01 0 0 0; 0 0.008 0 0; 0 0 0.0001 0; 0 0 0 0.0001]; %��������� ��������� ��� ���������� ���-��
%�������������� ������� ���� ����������
Rn =  0.008;

for k = 1 : K
    H = [1; cos(2*pi*k*(1/Thetta(4))+Thetta(3)); -sin(2*pi*k*(1/Thetta(4))+Thetta(3)); (2*pi*Thetta(2)*k*(-sin(Thetta(3) + (2*pi*k)/Thetta(4))))/Thetta(4)^2]'; %������ �� ����������� �� ���������� (��. ������ � ���. 16)
    P = Rpr*H'*inv(H*Rpr*H'+Rn);
    Thetta = Thetta + P * (signal(k) - (Thetta(1)+Thetta(2)*cos(2*pi*k*(1/Thetta(4))+Thetta(3))));
    Thetta(2) = abs(Thetta(2));
    B_est(k) = Thetta(1);
    A_est(k) = Thetta(2)^2;
    ph_init_est(k) = Thetta(3);
    T_est(k) = Thetta(4);

   % ph_init_est(k) = mod( ph_init_est(k), 2*3.14); 
end
A_est_sm = smooth(A_est,15);
B_est_sm = smooth(B_est,25);

%signal_sub3 = signal_sub2 - B_est;
%signal_sub3 = signal_if - signal_if2;

%plot(B, k);
subplot(4,1,1),plot(kk,B_est,kk, signal),title('���');
subplot(4,1,2),plot(kk,A_est, kk, signal),title('���������');
subplot(4,1,3),plot(kk,ph_init_est),title('��������� ����');
subplot(4,1,4),plot(kk,T_est),title('������');