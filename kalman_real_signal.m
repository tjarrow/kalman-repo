clc
close all;
clear all;

%%
%Моделируем сигнал
K = 500; %количество отсчетов
kk = zeros(1, K); %номера отсчетов для визуализации графиков 
B = zeros(1, K); %Фон
A = zeros(1, K); %Амплитудаa
ph_init = zeros(1,K); %начальная фаза
for k = 1 : K
    %B(k) = 0.2*sin(2*pi*k*0.002);
    %A(k) = exp(-((k-400)^2)/100^2)+0.6*exp(-((k-750)^2)/100^2);
    %s(k) = B(k) + A(k)*cos(2*pi*k*0.05 + ph_init) + n(k);
    kk(k) = k;
end
%figure(1), plot(s);
%plot(B);


%реальный сигнал
file_signal = fopen('signal.txt', 'r');
formatSpec = '%f';
signal = fscanf(file_signal, formatSpec);
signal_img = signal';

signal_f = fft(signal_img);
signal_f2 = fft(signal_img);
fft_idx = 100;
fft_idx2 = 420;
for i = 1:fft_idx;
    signal_f(i) = 0;
end;
for i = 1:fft_idx2
    signal_f2(i) = 0;
end;
signal_if = ifft(signal_f);
signal_if2 = ifft(signal_f2);
signal_sub2 = signal_if; % - signal_if2;

 devA = zeros(1,50);
 devB = zeros(1,50);
 devPhi = zeros(1,50);
 devT = zeros(1,50);
 
 
 for i = 1:50
 
 deviationA = zeros(1,K);
    deviationB = zeros(1,K);
    deviationPhi = zeros(1,K);
    deviationT = zeros(1,K);
    
    dev_sumA = 0;
    dev_sumB = 0;
    dev_sumPhi = 0;
    dev_sumT = 0;

%%
% начинаем фильтровать параметры сигнала
B_est = zeros(1, K); %Оценка фона
A_est = zeros(1, K); %Оценка амплитуды
ph_init_est = zeros(1, K); %оценка фазы

%начальные значения
B_est(1) = 0;
A_est(1) = 0;
ph_init_est(1) = 0;
Thetta = [B_est(1); A_est(1);ph_init_est(1) ]; %начальное значение вектора параметров

%ковариационная матрица шума системы (параметров)
Rpr = [0.073 0 0; 0 0.021 0; 0 0 0.006];
%ковариационная матрица шума наблюдения
Rn =  0.01;

for k = 1 : K
    H = [1; cos(2*pi*k*0.33+Thetta(3)); -sin(2*pi*k*0.33+Thetta(3))]'; %вектор из производных по параметрам (см. модель в стр. 16)
    P = Rpr*H'*inv(H*Rpr*H'+Rn);
    Thetta = Thetta + P * (signal_sub2(k) - (Thetta(1)+Thetta(2)*cos(2*pi*k*0.33+Thetta(3))));
    Thetta(2) = abs(Thetta(2));
    B_est(k) = Thetta(1);
    A_est(k) = Thetta(2);
    ph_init_est(k) = Thetta(3);
   % ph_init_est(k) = mod( ph_init_est(k), 2*3.14); 
end
A_est_sm = smooth(A_est,15);
B_est_sm = smooth(B_est,25);

signal_sub3 = signal_sub2 - B_est;

for j = 1:K
      deviationA(j) = ((A_est_sm(j) - A(j))^2);
    end

    dev_sumA = sum(deviationA)/K;
    devA(i) = dev_sumA;

    for j = 1:K
      deviationB(j) = (B_est_sm(j)^2);
    end

    dev_sumB = sum(deviationB)/K;
    devB(i) = dev_sumB;

    %for j = 1:K
      %deviationPhi(j) = ((ph_init_est(j) - (2*3.14*(1/T)*x)^2));
    %end

    %dev_sumPhi = sum(deviationPhi)/K;
    devPhi(i) = std(ph_init_est);
    
    end;
    
    
devA_t = devA';
devB_t = devB';
devPhi_t = devPhi';
devT_t = devT';



%plot(B, k);
figure (2),
subplot(3,1,1),plot(kk,B_est_sm,kk,signal_sub2 ),title('Фон');
subplot(3,1,2),plot(kk,A_est_sm, kk,signal_sub3),title('Амплитуда');
subplot(3,1,3),plot(kk,ph_init_est),title('Начальная фаза');