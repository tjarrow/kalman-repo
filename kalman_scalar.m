clc
close all;
clear all;

%%
%���������� ������
T_vec = zeros(1,50);
T = 1.8;
K = round(T*50);%���������� ��������

DEV_A = zeros(50, 50);
DEV_B = zeros(50, 50);
DEV_PHI = zeros(50, 50);
DEV_T = zeros(50, 50);

for k = 1 : K
    %B(k) = 0.2*sin(2*pi*k*0.002);
    %A(k) = exp(-((k-400)^2)/100^2)+0.6*exp(-((k-750)^2)/100^2);
    %s(k) = B(k) + A(k)*cos(2*pi*k*0.05 + ph_init) + n(k);
    kk(k) = k;
end

sigmaPsi=0.05; %�������� ����������� (������ ������)
sigmaEtaModel=0.01; %������ ��������� �������
kk = zeros(1, K); %������ �������� ��� ������������ �������� 

%�������������� ������� ���� ������� (����������)
%Rpr = [0.3 0 0 0; 0 0.021 0 0; 0 0 0.006 0; 0 0 0 0.1]; %��������� ��������� ��� ���������� ���-��
Rpr = [0.1 0 0 0; 0 0.5 0 0; 0 0 0.015 0; 0 0 0 0.05];

%�������������� ������� ���� ����������
Rn =  0.01;

for n = 1:50

%figure(1), plot(s);
%plot(B);

 B = zeros(1, K); %���
 A = zeros(1, K); %���������
 y = zeros(1, K);
 z = zeros(1, K);
 %ph_init = zeros(1,K); %��������� ����
 
 devA = zeros(1,50);
 devB = zeros(1,50);
 devPhi = zeros(1,50);
 devT = zeros(1,50);

for i = 1:50 % 50 ������������� (�.�. ������ ��������������)
    
    for x = 1:K
        A(x) = exp((-(x-(K/2))^2)/(K/1.5));
        y(x) =A(x)*cos(2*3.14*(1/T)*x)+normrnd(0,sigmaPsi); %T - ���-�� ������� �� ������
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
    Thetta = [B_est(1); A_est(1); ph_init_est(1); T_est(1) ]; %��������� �������� ������� ����������

    for k = 1 : K
        H = [1; cos(2*pi*k*(1/Thetta(4))+Thetta(3)); -Thetta(2)*sin(2*pi*k*(1/Thetta(4))+Thetta(3)); (2*pi*Thetta(2)*k*(-sin(Thetta(3) + (2*pi*k)/Thetta(4))))/Thetta(4)^2]'; %������ �� ����������� �� ���������� (��. ������ � ���. 16)
        P = Rpr*H'*inv(H*Rpr*H'+Rn);                                                        
        Thetta = Thetta + P * (z(k) - (Thetta(1)+Thetta(2)*cos(2*pi*k*(1/Thetta(4))+Thetta(3))) ); %d/d(thetta4)) (Thetta(1)+Thetta(2)*cos(2*pi*k*(1/Thetta(4))+Thetta(3))
        Thetta(2) = abs(Thetta(2));
        Thetta(4) = abs(Thetta(4)); % ���������: ���-�� �������� �� ����� ���� �������������
        B_est(k) = Thetta(1);
        A_est(k) = Thetta(2);
        ph_init_est(k) = Thetta(3);
        T_est(k) = Thetta(4);
       % ph_init_est(k) = mod( ph_init_est(k), 2*3.14); 
    end
    %A_est_sm = smooth(A_est,15);
    %B_est_sm = smooth(B_est,25);

    %signal_sub3 = signal_sub2 - B_est;
    signal_sub3 = z - B_est;

    for j = (K/2)- 25 :(K/2) + 25
      deviationA(j) = ((A_est(j) - A(j))^2);
    end

    dev_sumA = sum(deviationA)/K;
    devA(i) = dev_sumA;

    for j = (K/2)- 25 :(K/2) + 25
      deviationB(j) = ((B_est(j) - B(j))^2);
    end

    dev_sumB = sum(deviationB)/K;
    devB(i) = dev_sumB;

    %for j = 1:K
      %deviationPhi(j) = ((ph_init_est(j) - (2*3.14*(1/T)*x)^2));
    %end

    %dev_sumPhi = sum(deviationPhi)/K;
    devPhi(i) = var(ph_init_est); %std - ��� ����������� ����������, var - ���������
    
    for j = (K/2)- 25 :(K/2) + 25
      deviationT(j) = ((T_est(j) - T)^2);
    end
    
    dev_sumT = sum(deviationT)/K;
    devT(i) = dev_sumT;
    
end;
T_vec(n)= T;

if (T < 3)
   T = T + 0.1;
else 
   T = T + 0.5;
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

end;
figure(1);
set(gcf,'color','w');
boxplot(DEV_A, T_vec), xlabel('���������� �������� �� ������ �������'), ylabel('�������� �����������, ���.��.');
figure(2);
set(gcf,'color','w');
boxplot(DEV_B, T_vec), xlabel('���������� �������� �� ������ �������'), ylabel('�������� �����������, ���.��.');
figure(3);
set(gcf,'color','w');
boxplot(DEV_PHI, T_vec), xlabel('���������� �������� �� ������ �������'), ylabel('�������� �����������, ������');
figure(4);
set(gcf,'color','w');
boxplot(DEV_T, T_vec), xlabel('���������� �������� �� ������ �������'), ylabel('�������� �����������, ���.��.');
%plot(B, k);
%figure (2),
%subplot(3,1,1),plot(kk,B_est_sm,kk, z ),title('���');
%subplot(3,1,2),plot(kk,A_est_sm, kk,z),title('���������');
%subplot(3,1,3),plot(kk,ph_init_est),title('��������� ����');