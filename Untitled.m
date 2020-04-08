% � ���� ������� � �� ��������

T1 = 1.8; % T = 5 ( � ������ )
T2 = 5;

K = round(T1*50);%���������� ��������
kk = zeros(1, K); %������ �������� ��� ������������ �������� 

A = zeros(1, K); 
y = zeros(1, K);
z1 = zeros(1, K);
z2 = zeros(1, K);

for k = 1 : K
    kk(k) = k;
end

for x = 1:K
        B(x) = 0;
        A(x) = exp((-(x-(K/2))^2)/(K/1.5)); %��� ����� ������ �� ���� ������
        y(x) = B(x) + A(x)*cos(2*3.14*(1/T1)*x)+normrnd(0,sigmaPsi); % �� ��� ����� � - �������� ���������� ���. ����� ������ ���� ����� �, ������� ���-�� �� 0.01
        z1(x)=y(x)+normrnd(0,sigmaEtaModel); 
        idx2(x)= x;
 %      Thetta(1)+Thetta(2)*cos(2*pi*k*(1/Thetta(4))+Thetta(3))
end

T2 = 5;
K = round(T2*50);%���������� ��������
kk = zeros(1, K); %������ �������� ��� ������������ �������� 

A = zeros(1, K); 
y = zeros(1, K);
z2 = zeros(1, K);

for x = 1:K
        B(x) = 0;
        A(x) = exp((-(x-(K/2))^2)/(K/1.5)); %��� ����� ������ �� ���� ������
        y(x) = B(x) + A(x)*cos(2*3.14*(1/T2)*x)+normrnd(0,sigmaPsi); % �� ��� ����� � - �������� ���������� ���. ����� ������ ���� ����� �, ������� ���-�� �� 0.01
        z2(x)=y(x)+normrnd(0,sigmaEtaModel); 
        idx1(x)= x;
 %      Thetta(1)+Thetta(2)*cos(2*pi*k*(1/Thetta(4))+Thetta(3))
end

set(gcf,'color','w');
subplot(2,1,1),plot(idx2,z1),ylabel('�������� �������, ���.��.'), xlabel('����� ������� �������');
subplot(2,1,2),plot(idx1,z2);