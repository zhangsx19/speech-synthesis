clear all, close all, clc;
a1 = 1.3789;a2 = -0.9506;%为参数赋值
b = [1, 0, 0];%定义差分方程左侧系数
a = [1, -a1, -a2];%定义差分方程右侧系数
sT = 0.01/80;%采样是以每 10ms 采 80 个点进行的
sys = tf(b,a,sT,'variable','z^-1');%创造传递函数
sys
%residuez(B,A),B and A are the numerator(分子） and denominator polynomial coefficients
%合成模型和预测模型互为逆，所以a和b应该反过来
[r,p] = residuez(b,a);
f = abs(angle(p))/(2*pi*sT);
f
figure;
zplane(b,a),title('零极点分布图');%圆圈表示零点，×表示极点
figure;
freqz(b,a),title('频率响应');
figure;
subplot(2,1,1);
impz(b,a),title('impz单位样值响应');
subplot(2,1,2);
x = [1 zeros(1,399)];%x是单位样值激励
y = filter(b,a,x);
plot(y,'m'),title('filter单位样值响应'),xlabel("n(采样)"),ylabel("振幅");