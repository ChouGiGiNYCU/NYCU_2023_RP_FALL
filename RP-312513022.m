 %% parameter setting
clc;
clear;
% set Regulare Process 
noiseLength = 1000000; % length s[n],n=[1:1000000]
s=zeros([1 noiseLength+3]); % s inital value
a=[1,-1.9,1.18,-0.24]; % mean a0,a1,a2,a3
b0=2; 
% set white Gaussian noise value

meanValue = 0;     % mean value
whiteNoise = meanValue + randn(1, noiseLength);
varianceValue = var(whiteNoise);
disp('White Noise Variance :');
disp(varianceValue);

% find s[n] value
% Note matlab idx from 1 not 0
for i=4:noiseLength+3
    s(i)=b0*whiteNoise(i-3)-a(2)*s(i-1)-a(3)*s(i-2)-a(4)*s(i-3);
end

s=s(4:end); % remove init point
%% plot noise and process
figure("Name","whiteNoise");
subplot(2,2,1);
len=noiseLength/100;
plot(whiteNoise(1:len), 'o-', 'LineWidth', 1, 'DisplayName', 'whiteNoise');
hold on; 
plot(s(1:len), 'x-','LineWidth', 1, 'DisplayName', 'Process');
title('Signal');
xlabel('Index');
ylabel('Values');
legend('Location', 'best');
grid on;
hold off;

%% Problem-1

% Rss[0]+a1Rss[1]+a2Rss[2]+a3Rss[3]=b0^2
% Rss[1]+a1Rss[0]+a2Rss[1]+a3Rss[2]=0
%  -> a1Rss[0]+(1+a2)Rss[1]+a3Rss[2]=0
% Rss[2]+a1Rss[1]+a2Rss[0]+a3Rss[1]=0
%  -> a2Rss[0] + (a1+a3)Rss[1] + Rss[2]=0
% Rss[3]+a1Rss[2]+a2Rss[1]+a3Rss[0]=0
%  -> a3Rss[0] + a2Rss[1] + a1Rss[2] + Rss[3]=0
M=20; 
% slove 4-order equation at first
% Ax=b find x 、 x is Rss[0]~Rss[3]
A=[a(1),a(2),a(3),a(4);
   a(2),a(1)+a(3),a(4),0;
   a(3),a(2)+a(4),a(1),0;
   a(4),a(3),a(2),a(1)];
b=[b0^2;0;0;0];
Rss=transpose(A\b);
disp('Rss[0] ~ Rss[3]:');
disp(Rss);
Rss=[Rss,zeros([1,17])];
% slove Rss[4]~Rss[20]
for i=5:M+1
    % Rss[4(i)]=-a1*Rss[3(i-1)]-a2*Rss[2(i-2)]-a3*Rss[1(i-3)]
    Rss(i)=-a(2)*Rss(i-1)-a(3)*Rss(i-2)-a(4)*Rss(i-3);
end

% plot
subplot(2,2,2);
plot([0:M],Rss,'x-','LineWidth', 1, 'DisplayName', 'Rss');
title('Rss[m] ,  ∀ 0 ≤ m ≤ 20 ');
xlabel('Index(m)');
ylabel('Values');
legend('Location', 'best');
grid on;

%% Problem-2
n=[10000:10100];
% Count h[k] at first ∀ 1 ≤ k ≤ 20
% ^s[n]= ∑h[k]*s[n-k] ∀ 1 ≤ k ≤ 20
% R[1]= ∑ h[k]*R[1-k]   
% R[2]= ∑ h[k]*R[2-k]
%   .
%   .
% R[20]= ∑ h[k]*R[20-k]

% Ax=b、b=[R[1];....;R[20]]
b=transpose(Rss(2:21));
A=[];
for i=1:M
    tmp=flip(Rss(1:i));
    if M-i>0
        tmp=[tmp,Rss(2:M-i+1)];
    end
    A=[A;tmp];
end

if transpose(A)==A
    disp('Problem2 - A matrix is same');
end
h=transpose(A\b);

s_estimate = [];
for i=n
    tmp=h.*flip(s(i-M+1:i));
    s_estimate=[s_estimate,sum(tmp)];
end

% plot
subplot(2,2,3);
plot(n,s_estimate,'x-','LineWidth', 1, 'DisplayName', 's-estimate');
hold on; 
plot(n,s(n), 'o-','LineWidth', 1, 'DisplayName', 'Process');
title('s-estimate[n] ,  ∀ 10000 ≤ n ≤ 10100 ');
xlabel('Index(n)');
ylabel('Values');
legend('Location', 'best');
grid on;
%% Count all s_estimate
s_estimate=zeros([1 noiseLength]); % clean all
for n=1:M
    tmp=0;
    for k=1:M
        if n-k>0
            tmp=tmp+h(k)*s(n-k);
        end
    end
    s_estimate(n)=tmp;
end

for n=M+1:noiseLength
    s_estimate(n)=sum(h.*flip(s(n-M:n-1)));
end
%%
% plot
% subplot(3,2,4);
% plot(s_estimate(1:200),'x-','LineWidth', 1, 'DisplayName', 's-estimate');
% hold on; 
% plot(s(1:200), 'o-','LineWidth', 1, 'DisplayName', 'Process');
% title('s-estimate[n] ,  previous 200 point ');
% xlabel('Index(n)');
% ylabel('Values');
% legend('Location', 'best');
% grid on;
%%  Problem-3  n=10-1
p_10=0;
n_10=10;
for i=1:n_10
    tmp=(s_estimate(i)-s(i))^2;
    p_10=p_10+tmp;
end
p_10=p_10/10;
disp('p_10 :');
disp(p_10);
%%  Problem-3  n=100-2
p_100=0;
n_100=100;
for i=n_10+1:n_100
    tmp=(s_estimate(i)-s(i))^2;
    p_100=p_100+tmp;
end
p_100=((p_10*n_10) + p_100)/n_100;
disp('p_100 :');
disp(p_100);
%%  Problem-3  n=1000-3
p_1000=0;
n_1000=1000;
for i=n_100+1:n_1000
    tmp=(s_estimate(i)-s(i))^2;
    p_1000=p_1000+tmp;
end
p_1000=((p_100*n_100) + p_1000)/n_1000;
disp('p_1000 :');
disp(p_1000);
%%  Problem-3  n=10000-4
p_10000=0;
n_10000=10000;
for i=n_1000+1:n_10000
    tmp=(s_estimate(i)-s(i))^2;
    p_10000=p_10000+tmp;
end
p_10000=((p_1000*n_1000) + p_10000)/n_10000;
disp('p_10000 :');
disp(p_10000);
%%  Problem-3  n=100000-5
p_100000=0;
n_100000=100000;
for i=n_10000+1:n_100000
    tmp=(s_estimate(i)-s(i))^2;
    p_100000=p_100000+tmp;
end
p_100000=((p_10000*n_10000) + p_100000)/n_100000;
disp('p_100000 :');
disp(p_100000);
%%  Problem-3  n=1000000-6
p_1000000=0;
n_1000000=1000000;
for i=n_100000+1:n_1000000
    tmp=(s_estimate(i)-s(i))^2;
    p_1000000=p_1000000+tmp;
end
p_1000000=((p_100000*n_100000) + p_1000000)/n_1000000;
disp('p_1000000 :');
disp(p_1000000);
%% plot  p_10~p_1000000
figure;
stem([p_10 p_100 p_1000 p_10000 p_100000 p_1000000],'x-','LineWidth', 1, 'DisplayName', 'P');

title('P  ,  n=10~1000000 ');
xlabel('Index(n)');
ylabel('Values');
legend('Location', 'best');
grid on;
%%  Problem-4
% 参数设置
clear;
clc;
noiseLength=1000000;
s=zeros([1 noiseLength+3]); % s inital value
a=[1,-1.9,1.18,-0.24]; % mean a0,a1,a2,a3
b0=2;

m=1;
mu = 0; 
sigma = sqrt(2);
u = rand(m, noiseLength)-0.5;
b = sigma / sqrt(2);
laplaceRandomNoise = mu - b * sign(u).* log(1- 2* abs(u));

for i=4:noiseLength+3
    s(i)=b0*laplaceRandomNoise(i-3)-a(2)*s(i-1)-a(3)*s(i-2)-a(4)*s(i-3);
end

s=s(4:end); % remove init point
%% plot noise and process
figure("Name","LaplaceRandomNoise");
subplot(2,2,1);
len=noiseLength/100;
plot(laplaceRandomNoise(1:len), 'o-', 'LineWidth', 1, 'DisplayName', 'LaplaceRandomNoise');
hold on; 
plot(s(1:len), 'x-','LineWidth', 1, 'DisplayName', 'Process');
title('Signal');
xlabel('Index');
ylabel('Values');
legend('Location', 'best');
grid on;
hold off;
%% Count Rss
M=20; 
% Ax=b find x 、 x is Rss[0]~Rss[3]
A=[a(1),a(2),a(3),a(4);
   a(2),a(1)+a(3),a(4),0;
   a(3),a(2)+a(4),a(1),0;
   a(4),a(3),a(2),a(1)];
b=[b0^2;0;0;0];
Rss=transpose(A\b);
disp('Rss[0] ~ Rss[3]:');
disp(Rss);
Rss=[Rss,zeros([1,17])];
% slove Rss[4]~Rss[20]
for i=5:M+1
    % Rss[4(i)]=-a1*Rss[3(i-1)]-a2*Rss[2(i-2)]-a3*Rss[1(i-3)]
    Rss(i)=-a(2)*Rss(i-1)-a(3)*Rss(i-2)-a(4)*Rss(i-3);
end

% plot
subplot(2,2,2);
plot([0:M],Rss,'x-','LineWidth', 1, 'DisplayName', 'Rss');
title('Rss[m] ,  ∀ 0 ≤ m ≤ 20 ');
xlabel('Index(m)');
ylabel('Values');
legend('Location', 'best');
grid on;
%%
n=[10000:10100];
% Count h[k] at first ∀ 1 ≤ k ≤ 20
% ^s[n]= ∑h[k]*s[n-k] ∀ 1 ≤ k ≤ 20
% R[1]= ∑ h[k]*R[1-k]   
% R[2]= ∑ h[k]*R[2-k]
%   .
%   .
% R[20]= ∑ h[k]*R[20-k]

% Ax=b、b=[R[1];....;R[20]]
b=transpose(Rss(2:21));
A=[];
for i=1:M
    tmp=flip(Rss(1:i));
    if M-i>0
        tmp=[tmp,Rss(2:M-i+1)];
    end
    A=[A;tmp];
end

if transpose(A)==A
    disp('Problem2 - A matrix is same');
end
h=transpose(A\b);

s_estimate = [];
for i=n
    tmp=h.*flip(s(i-M+1:i));
    s_estimate=[s_estimate,sum(tmp)];
end

% plot
subplot(2,2,3);
plot(n,s_estimate,'x-','LineWidth', 1, 'DisplayName', 's-estimate');
hold on; 
plot(n,s(n), 'o-','LineWidth', 1, 'DisplayName', 'Process');
title('s-estimate[n] ,  ∀ 10000 ≤ n ≤ 10100 ');
xlabel('Index(n)');
ylabel('Values');
legend('Location', 'best');
grid on;
%% Count all s_estimate
s_estimate=zeros([1 noiseLength]); % clean all
for n=1:M
    tmp=0;
    for k=1:M
        if n-k>0
            tmp=tmp+h(k)*s(n-k);
        end
    end
    s_estimate(n)=tmp;
end

for n=M+1:noiseLength
    s_estimate(n)=sum(h.*flip(s(n-M:n-1)));
end
%%
% plot
% subplot(3,2,4);
% plot(s_estimate(1:200),'x-','LineWidth', 1, 'DisplayName', 's-estimate');
% hold on; 
% plot(s(1:200), 'o-','LineWidth', 1, 'DisplayName', 'Process');
% title('s-estimate[n] ,  previous 200 point ');
% xlabel('Index(n)');
% ylabel('Values');
% legend('Location', 'best');
% grid on;
%%  Problem-4  n=10-1
p_10=0;
n_10=10;
for i=1:n_10
    tmp=(s_estimate(i)-s(i))^2;
    p_10=p_10+tmp;
end
p_10=p_10/10;
disp('p_10 :');
disp(p_10);
%%  Problem-4  n=100-2
p_100=0;
n_100=100;
for i=n_10+1:n_100
    tmp=(s_estimate(i)-s(i))^2;
    p_100=p_100+tmp;
end
p_100=((p_10*n_10) + p_100)/n_100;
disp('p_100 :');
disp(p_100);
%%  Problem-4  n=1000-3
p_1000=0;
n_1000=1000;
for i=n_100+1:n_1000
    tmp=(s_estimate(i)-s(i))^2;
    p_1000=p_1000+tmp;
end
p_1000=((p_100*n_100) + p_1000)/n_1000;
disp('p_1000 :');
disp(p_1000);
%%  Problem-4  n=10000-4
p_10000=0;
n_10000=10000;
for i=n_1000+1:n_10000
    tmp=(s_estimate(i)-s(i))^2;
    p_10000=p_10000+tmp;
end
p_10000=((p_1000*n_1000) + p_10000)/n_10000;
disp('p_10000 :');
disp(p_10000);
%%  Problem-4  n=100000-5
p_100000=0;
n_100000=100000;
for i=n_10000+1:n_100000
    tmp=(s_estimate(i)-s(i))^2;
    p_100000=p_100000+tmp;
end
p_100000=((p_10000*n_10000) + p_100000)/n_100000;
disp('p_100000 :');
disp(p_100000);
%%  Problem-4  n=1000000-6
p_1000000=0;
n_1000000=1000000;
for i=n_100000+1:n_1000000
    tmp=(s_estimate(i)-s(i))^2;
    p_1000000=p_1000000+tmp;
end
p_1000000=((p_100000*n_100000) + p_1000000)/n_1000000;
disp('p_1000000 :');
disp(p_1000000);
%% plot  p_10~p_1000000
figure;
stem([p_10 p_100 p_1000 p_10000 p_100000 p_1000000],'x-','LineWidth', 1, 'DisplayName', 'P');

title('P  ,  n=10~1000000 ');
xlabel('Index(n)');
ylabel('Values');
legend('Location', 'best');
grid on;