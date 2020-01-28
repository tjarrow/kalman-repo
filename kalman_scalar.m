clc
close all;
clear all;

%%
%Моделируем сигнал

T = 6;
K = round(T*50);%количество отсчетов

DEV_A = zeros(50, 10);
DEV_B = zeros(50, 10);
DEV_PHI = zeros(50, 10);
DEV_T = zeros(50, 10);

for k = 1 : K
    %B(k) = 0.2*sin(2*pi*k*0.002);
    %A(k) = exp(-((k-400)^2)/100^2)+0.6*exp(-((k-750)^2)/100^2);
    %s(k) = B(k) + A(k)*cos(2*pi*k*0.05 + ph_init) + n(k);
    kk(k) = k;
end


sigmaPsi=0.05; %реальные погрешности (ошибка модели)
sigmaEtaModel=0.1; %ошибка измерений прибора
sigmaEtaFilter = 0.8;
kk = zeros(1, K); %номера отсчетов для визуализации графиков 

%ковариационная матрица шума системы (параметров)
Rpr = [0.073 0 0 0; 0 0.021 0 0; 0 0 0.006 0; 0 0 0 0.0001]; %подобрать дисперсию для последнего пар-ра
%ковариационная матрица шума наблюдения
Rn =  0.01;

for n = 1:10


%figure(1), plot(s);
%plot(B);

 %B = zeros(1, K); %Фон
 A = zeros(1, K); %Амплитуда
 y = zeros(1, K);
 z = zeros(1, K);
 %ph_init = zeros(1,K); %начальная фаза
 
 devA = zeros(1,50);
 devB = zeros(1,50);
 devPhi = zeros(1,50);
 devT = zeros(1,50);


for i = 1:50
    
    for x = 1:K
        A(x) = exp((-(x-45)^2)/150);
        y(x) = A(x)*cos(2*3.14*(1/T)*x)+normrnd(0,sigmaPsi); %T - кол-во отчётов на период
        z(x)=y(x)+normrnd(0,sigmaEtaModel); 
    end
   
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
    T_est = zeros(1, K); %оценка периода

    %начальные значения
    B_est(1) = 0;
    A_est(1) = 0;
    ph_init_est(1) = 0;
    T_est(1) = 5;
    Thetta = [B_est(1); A_est(1); ph_init_est(1); T_est(1) ]; %начальное значение вектора параметров

    for k = 1 : K
        H = [1; cos(2*pi*k*(1/Thetta(4))+Thetta(3)); -sin(2*pi*k*(1/Thetta(4))+Thetta(3)); (2*pi*Thetta(2)*k*sin(Thetta(3) + (2*pi*k)/Thetta(4)))/Thetta(4)^2]'; %вектор из производных по параметрам (см. модель в стр. 16)
        P = Rpr*H'*inv(H*Rpr*H'+Rn);                                                        
        Thetta = Thetta + P * (z(k) - (Thetta(1)+Thetta(2)*cos(2*pi*k*(1/Thetta(4))+Thetta(3)) ) ); %d/d(thetta4)) (Thetta(1)+Thetta(2)*cos(2*pi*k*(1/Thetta(4))+Thetta(3))
        Thetta(2) = abs(Thetta(2));
        Thetta(4) = abs(Thetta(4));
        B_est(k) = Thetta(1);
        A_est(k) = Thetta(2);
        ph_init_est(k) = Thetta(3);
        T_est(k) = Thetta(4);
       % ph_init_est(k) = mod( ph_init_est(k), 2*3.14); 
    end
    A_est_sm = smooth(A_est,15);
    B_est_sm = smooth(B_est,25);

    %signal_sub3 = signal_sub2 - B_est;
    signal_sub3 = z - B_est;

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
    
    for j = 1:K
      deviationT(j) = (T_est(j)^2);
    end
    
    dev_sumT = sum(deviationT)/K;
    devT(i) = dev_sumT;
    
   
end;

devA_t = devA';
devB_t = devB';
devPhi_t = devPhi';
devT_t = devT';

for f = 1:50
    DEV_A(f, n) = devA_t(f, :);
    DEV_B(f, n) = devB_t(f, :);
    DEV_PHI(f, n) = devPhi_t(f, :);
    DEV_T(f, n) = devT_t(f, :);
end;

T = T + 2;

end;
%boxplot(devB, T);

%plot(B, k);
figure (2),
%subplot(3,1,1),plot(kk,B_est_sm,kk, z ),title('Фон');
%subplot(3,1,2),plot(kk,A_est_sm, kk,z),title('Амплитуда');
%subplot(3,1,3),plot(kk,ph_init_est),title('Начальная фаза');