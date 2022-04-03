clear all, close all, clc;
a1 = 1.3789;a2 = -0.9506;%为参数赋值
%为方便起见，把s(n)看做输出，e(n)看做输入
a = [1, -a1, -a2];%定义差分方程左侧系数
b = [1, 0, 0];%定义差分方程右侧系数
[r,p,k] = residuez(b,a);%求出原极点
angle = 150 * 2 * pi / 8000;%θ
p = p .* exp(sign(imag(p)) *1i *  angle);%已知150Hz对应的极点
[b,a] = residue(r,p,k);
a
