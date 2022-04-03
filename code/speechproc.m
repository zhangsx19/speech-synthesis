function speechproc()

    % 定义常数
    FL = 80;                % 帧长
    WL = 240;               % 窗长
    P = 10;                 % 预测系数个数
    s = readspeech('voice.pcm',100000);             % 载入语音s
    L = length(s);          % 读入语音长度
    FN = floor(L/FL)-2;     % 计算帧数
    % 预测和重建滤波器
    exc = zeros(L,1);       % 激励信号（预测误差）
    zi_pre = zeros(P,1);    % 预测滤波器的状态
    s_rec = zeros(L,1);     % 重建语音
    zi_rec = zeros(P,1);
    % 合成滤波器
    pos = 2*FL+1;           % 激励信号起始位置
    zi_syn = zeros(P, 1);   %滤波器状态
    exc_syn = zeros(L,1);   % 合成的激励信号（脉冲串）
    s_syn = zeros(L,1);     % 合成语音
    % 变调不变速滤波器
    pos_t = 2 * FL + 1;
    zi_syn_t = zeros(P, 1);
    exc_syn_t = zeros(L,1);   % 合成的激励信号（脉冲串）
    s_syn_t = zeros(L,1);     % 合成语音
    % 变速不变调滤波器（假设速度减慢一倍）
    FL_v = FL * 2;
    pos_v = 2 * FL_v + 1;
    zi_syn_v = zeros(P, 1);
    exc_syn_v = zeros(2*L,1);   % 合成的激励信号（脉冲串）
    s_syn_v = zeros(2*L,1);     % 合成语音
    
    hw = hamming(WL);       % 汉明窗
    % 依次处理每帧语音
    for n = 3:FN

        % 计算预测系数（不需要掌握）
        s_w = s(n*FL-WL+1:n*FL).*hw;    %汉明窗加权后的语音
        [A E] = lpc(s_w, P);            %用线性预测法计算P个预测系数
                                        % A是预测系数，E会被用来计算合成激励的能量

        if n == 27
        % (3) 在此位置写程序，观察预测系统的零极点图
        %对预测系统，e(n)是输出，s(n)是输入，e(n)=s(n)-∑ak*s(n-k)
            figure;
            zplane(A,1),title('零极点分布图');%A对应差分方程右侧系数
        end
        
        s_f = s((n-1)*FL+1:n*FL);       % 本帧语音，下面就要对它做处理

        % (4) 在此位置写程序，用filter函数和s_f计算激励，注意保持滤波器状态
        [e_pre,zf_pre] = filter(A,1,s_f,zi_pre);
        %左边的zf_pre为本次循环得到的滤波器的最终状态，需要以其作为下一次循环的filter函数的zi初始状态
        %第一次循环的zi为0
        zi_pre = zf_pre;
        exc((n-1)*FL+1:n*FL) = e_pre;%将你计算得到的激励写在这里

        % (5) 在此位置写程序，用filter函数和exc重建语音，注意保持滤波器状态
        [srec,zf_rec] = filter(1,A,exc((n-1)*FL+1:n*FL),zi_rec);
        %由激励exc得到重建语音srec,同样要维持滤波器状态不变，用的是zi_rec和zf_rec
        zi_rec = zf_rec;
        s_rec((n-1)*FL+1:n*FL) = srec;%将你计算得到的重建语音写在这里

        % 注意下面只有在得到exc后才会计算正确
        s_Pitch = exc(n*FL-222:n*FL);
        PT = findpitch(s_Pitch);    % 计算基音周期PT（不要求掌握）
        G = sqrt(E*PT);           % 计算合成激励的能量G（不要求掌握）

        
        % (10) 在此位置写程序，生成合成激励，并用激励和filter函数产生合成语音
        %仿照1.2.2的（8），不同的是PT是在单帧中是固定的常数
        while pos <= n * FL
            exc_syn(pos) = G;
            pos = pos + PT;
        end
        %exc_syn((n-1)*FL+1:n*FL) =  将你计算得到的合成激励写在这里
        %仿照（5），由激励exc_syn得到合成语音s_syn,同样要维持滤波器状态不变，用的是zi_syn和zf_syn
        [s_syn1,zf_syn] = filter(1,A,exc_syn((n-1)*FL+1:n*FL),zi_syn);
        zi_syn = zf_syn;
        s_syn((n-1)*FL+1:n*FL) = s_syn1;%将你计算得到的合成语音写在这里

        % (11) 不改变基音周期和预测系数，将合成激励的长度增加一倍，再作为filter
        % 的输入得到新的合成语音，听一听是不是速度变慢了，但音调没有变。
        %仿照（10），不同的是FL变成了FL_v
        while pos_v <= n * FL_v
            exc_syn_v(pos_v) = G;
            pos_v = pos_v + PT;
        end
        [s_syn1_v,zf_syn_v] = filter(1,A,exc_syn_v((n-1)*FL_v+1:n*FL_v),zi_syn_v);
        zi_syn_v = zf_syn_v;
        s_syn_v((n-1)*FL_v+1:n*FL_v) = s_syn1_v;%将你计算得到的加长合成语音写在这里
        % exc_syn_v((n-1)*FL_v+1:n*FL_v) = ... 将你计算得到的加长合成激励写在这里
        
        % (13) 将基音周期减小一半，将共振峰频率增加150Hz，重新合成语音，听听是啥感受～
        while pos_t <= n * FL
            exc_syn_t(pos_t) = G;
            pos_t = pos_t + round(PT/2);
        end
        [r,p,k] = residuez(1,A);%求出原极点,生成模型中A是差分方程左侧系数
        angle = 150 * 2 * pi / 8000;%θ
        p = p .* exp(sign(imag(p)) *1i *  angle);%已知150Hz对应的极点
        [B_t,A_t] = residue(r,p,k);
        [s_syn1_t,zf_syn_t] = filter(1,A_t,exc_syn_t((n-1)*FL+1:n*FL),zi_syn_t);
        zi_syn_t = zf_syn_t;
        s_syn_t((n-1)*FL+1:n*FL) = s_syn1_t;%将你计算得到的变调合成语音写在这里
        
    end

    % (6) 在此位置写程序，听一听 s ，exc 和 s_rec 有何区别，解释这种区别
    % 后面听语音的题目也都可以在这里写，不再做特别注明
    Fs = 8000;%采样频率,采样是以每 10ms 采 80 个点进行的
    sound([exc;s;s_rec]/2^15,Fs,16);
    t = (0:L-1) / Fs; % 由采样率生成对应时间
    figure;
    subplot(3,1,1);plot(t,exc),title('预测的激励信号e(n)'),xlabel("时间/s"),ylabel("振幅");
    subplot(3,1,2);plot(t,s),title('原始语音信号s(n)'),xlabel("时间/s"),ylabel("振幅");
    subplot(3,1,3);plot(t,s_rec),title('重建的语音信号srec'),xlabel("时间/s"),ylabel("振幅");
    tpart = (9000:10000) / Fs; % 截取一小段
    figure;
    subplot(3,1,1);plot(tpart,exc(9000:10000)),title('预测的激励信号e(n)'),xlabel("时间/s"),ylabel("振幅");
    subplot(3,1,2);plot(tpart,s(9000:10000)),title('原始语音信号s(n)'),xlabel("时间/s"),ylabel("振幅");
    subplot(3,1,3);plot(tpart,s_rec(9000:10000)),title('重建的语音信号srec'),xlabel("时间/s"),ylabel("振幅");
    %试听变化的单位样值串
    pause(6);
    last_time = 1;
    s1 = unit_sample(200,last_time,Fs);
    s2 = unit_sample(300,last_time,Fs);
    sound([s1;s2],Fs);
    figure;
    subplot(2,1,1);plot(s1);title('200Hz');
    subplot(2,1,2);plot(s2);title('300Hz');
    %试听单位样值串
    pause(3);
    last_time = 1;
    s3 = MultiUnitSample(last_time,Fs);
    sound(s3,Fs);
    figure;
    plot(s3);title('MultiUnitSample');
    %把MultiUnitSample作为输入e(n)输入到ex_1_1.m的滤波器，输出为s(n)
    a1 = 1.3789;a2 = -0.9506;%为参数赋值
    b = [1, 0, 0];%定义差分方程右侧系数(e(n))
    a = [1, -a1, -a2];%定义差分方程左侧系数(s(n))
    s_multi = filter(b,a,s3);
    pause(3);
    sound([s_multi/max(abs(s_multi));s3],Fs);
    figure;
    subplot(4,1,1);plot(s_multi/max(abs(s_multi)));title('时域输出s(n)');
    subplot(4,1,2);plot(s3);title('时域输入e(n)');
    subplot(4,1,3);fft_plot(s_multi/max(abs(s_multi)),Fs);title('频域输出s(n)');
    subplot(4,1,4);fft_plot(s3,Fs);title('频域输入e(n)');
    %试听合成语音s_syn和原始语音s的区别
    pause(3);
    sound([s_syn;s]/2^15,Fs);
    figure;
    subplot(4,1,1);plot(s_syn);title('合成语音时域');
    subplot(4,1,2);plot(s);title('原始语音时域');
    subplot(4,1,3);fft_plot(s_syn,Fs);title('合成语音频域');
    subplot(4,1,4);fft_plot(s,Fs);title('原始语音频域');
    %变速不变调
    pause(3);
    sound(s_syn_v/2^15,Fs);
    figure;
    subplot(4,1,1);plot(s_syn_v);title('变速不变调语音时域');
    subplot(4,1,2);plot(s);title('原始语音时域');
    subplot(4,1,3);fft_plot(s_syn_v,Fs);title('变速不变调语音频域');
    subplot(4,1,4);fft_plot(s,Fs);title('原始语音频域');
    %变调不变速
    pause(6);
    sound(s_syn_t/2^15,Fs);
    figure;
    subplot(4,1,1);plot(s_syn_t);title('变调不变速语音时域');
    subplot(4,1,2);plot(s);title('原始语音时域');
    subplot(4,1,3);fft_plot(s_syn_t,Fs);title('变调不变速语音频域');
    subplot(4,1,4);fft_plot(s,Fs);title('原始语音频域');
    % 保存所有文件
    writespeech('exc.pcm',exc);
    writespeech('rec.pcm',s_rec);
    writespeech('exc_syn.pcm',exc_syn);
    writespeech('syn.pcm',s_syn);
    writespeech('exc_syn_t.pcm',exc_syn_t);
    writespeech('syn_t.pcm',s_syn_t);
    writespeech('exc_syn_v.pcm',exc_syn_v);
    writespeech('syn_v.pcm',s_syn_v);
return

% 从PCM文件中读入语音
function s = readspeech(filename, L)
    fid = fopen(filename, 'r');
    s = fread(fid, L, 'int16');
    fclose(fid);
return

% 写语音到PCM文件中
function writespeech(filename,s)
    fid = fopen(filename,'w');
    fwrite(fid, s, 'int16');
    fclose(fid);
return

% 计算一段语音的基音周期，不要求掌握
function PT = findpitch(s)
[B, A] = butter(5, 700/4000);
s = filter(B,A,s);
R = zeros(143,1);
for k=1:143
    R(k) = s(144:223)'*s(144-k:223-k);
end
[R1,T1] = max(R(80:143));
T1 = T1 + 79;
R1 = R1/(norm(s(144-T1:223-T1))+1);
[R2,T2] = max(R(40:79));
T2 = T2 + 39;
R2 = R2/(norm(s(144-T2:223-T2))+1);
[R3,T3] = max(R(20:39));
T3 = T3 + 19;
R3 = R3/(norm(s(144-T3:223-T3))+1);
Top = T1;
Rop = R1;
if R2 >= 0.85*Rop
    Rop = R2;
    Top = T2;
end
if R3 > 0.85*Rop
    Rop = R3;
    Top = T3;
end
PT = Top;
return