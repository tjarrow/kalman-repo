clc
close all;
clear all;

% � ���� ������� � �� ��������

T = 5;
K = round(T*100);%���������� ��������

sigmaPsi=0.05; %�������� ����������� (������ ������)
sigmaEtaModel=0.1; %������ ��������� �������
sigmaEtaFilter = 0.8;
kk = zeros(1, K); %������ �������� ��� ������������ �������� 

%�������������� ������� ���� ������� (����������)
Rpr = [0.073 0 0 0; 0 0.021 0 0; 0 0 0.006 0; 0 0 0 0.0001]; %��������� ��������� ��� ���������� ���-��
%�������������� ������� ���� ����������
Rn =  0.01;

 %B = zeros(1, K); %���
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
        A(x) = exp((-(x-500)^2)/50000);
        y(x) = A(x)*cos(2*3.14*(1/T)*x)+normrnd(0,sigmaPsi); %T - ���-�� ������� �� ������
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
    T_est(1) = 3;
    Thetta = [B_est(1); A_est(1); ph_init_est(1); T_est(1) ]; %��������� �������� ������� ����������

    for k = 1 : K
        H = [1; cos(2*pi*k*(1/Thetta(4))+Thetta(3)); -Thetta(2)*sin(2*pi*k*(1/Thetta(4))+Thetta(3)); (2*pi*Thetta(2)*k*sin(2*pi*k*(Thetta(3) + 1/Thetta(4)))) / Thetta(4)^2]'; %������ �� ����������� �� ���������� (��. ������ � ���. 16)
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
      deviationT(j) = ((T_est(j) - T)^2);
    end
    
    dev_sumT = sum(deviationT)/K;
    devT(i) = dev_sumT;
    
   
end;

devA_t = devA';
devB_t = devB';
devPhi_t = devPhi';
devT_t = devT';
