function signal = MultiUnitSample(last_time, Fs)
%生成基音周期随时间变化的单位样值串
% last_time: 信号持续时间
% Fs: 采样频率
L = round(Fs * last_time);%信号长度
signal = zeros(L, 1);
n = 1;
while n <= L
    signal(n) = 1;
    m = floor(n/80);%10ms的分段
    PT = 80 + 5 * mod(m, 50);
    n = n + PT;
end

end