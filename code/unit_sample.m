function signal = unit_sample(f,last_time,Fs)
%生成特定频率、持续时间和采样频率的单位样值串
% last_time: 信号持续时间
% Fs: 采样频率
% f: 单位样值串的频率
signal = zeros(round(Fs * last_time), 1);
NS = round(last_time * f);
N = round(Fs / f);
k = 0:NS - 1;
signal(k * N + 1) = 1;
end