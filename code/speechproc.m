function speechproc()

    % ���峣��
    FL = 80;                % ֡��
    WL = 240;               % ����
    P = 10;                 % Ԥ��ϵ������
    s = readspeech('voice.pcm',100000);             % ��������s
    L = length(s);          % ������������
    FN = floor(L/FL)-2;     % ����֡��
    % Ԥ����ؽ��˲���
    exc = zeros(L,1);       % �����źţ�Ԥ����
    zi_pre = zeros(P,1);    % Ԥ���˲�����״̬
    s_rec = zeros(L,1);     % �ؽ�����
    zi_rec = zeros(P,1);
    % �ϳ��˲���
    pos = 2*FL+1;           % �����ź���ʼλ��
    zi_syn = zeros(P, 1);   %�˲���״̬
    exc_syn = zeros(L,1);   % �ϳɵļ����źţ����崮��
    s_syn = zeros(L,1);     % �ϳ�����
    % ����������˲���
    pos_t = 2 * FL + 1;
    zi_syn_t = zeros(P, 1);
    exc_syn_t = zeros(L,1);   % �ϳɵļ����źţ����崮��
    s_syn_t = zeros(L,1);     % �ϳ�����
    % ���ٲ�����˲����������ٶȼ���һ����
    FL_v = FL * 2;
    pos_v = 2 * FL_v + 1;
    zi_syn_v = zeros(P, 1);
    exc_syn_v = zeros(2*L,1);   % �ϳɵļ����źţ����崮��
    s_syn_v = zeros(2*L,1);     % �ϳ�����
    
    hw = hamming(WL);       % ������
    % ���δ���ÿ֡����
    for n = 3:FN

        % ����Ԥ��ϵ��������Ҫ���գ�
        s_w = s(n*FL-WL+1:n*FL).*hw;    %��������Ȩ�������
        [A E] = lpc(s_w, P);            %������Ԥ�ⷨ����P��Ԥ��ϵ��
                                        % A��Ԥ��ϵ����E�ᱻ��������ϳɼ���������

        if n == 27
        % (3) �ڴ�λ��д���򣬹۲�Ԥ��ϵͳ���㼫��ͼ
        %��Ԥ��ϵͳ��e(n)�������s(n)�����룬e(n)=s(n)-��ak*s(n-k)
            figure;
            zplane(A,1),title('�㼫��ֲ�ͼ');%A��Ӧ��ַ����Ҳ�ϵ��
        end
        
        s_f = s((n-1)*FL+1:n*FL);       % ��֡�����������Ҫ����������

        % (4) �ڴ�λ��д������filter������s_f���㼤����ע�Ᵽ���˲���״̬
        [e_pre,zf_pre] = filter(A,1,s_f,zi_pre);
        %��ߵ�zf_preΪ����ѭ���õ����˲���������״̬����Ҫ������Ϊ��һ��ѭ����filter������zi��ʼ״̬
        %��һ��ѭ����ziΪ0
        zi_pre = zf_pre;
        exc((n-1)*FL+1:n*FL) = e_pre;%�������õ��ļ���д������

        % (5) �ڴ�λ��д������filter������exc�ؽ�������ע�Ᵽ���˲���״̬
        [srec,zf_rec] = filter(1,A,exc((n-1)*FL+1:n*FL),zi_rec);
        %�ɼ���exc�õ��ؽ�����srec,ͬ��Ҫά���˲���״̬���䣬�õ���zi_rec��zf_rec
        zi_rec = zf_rec;
        s_rec((n-1)*FL+1:n*FL) = srec;%�������õ����ؽ�����д������

        % ע������ֻ���ڵõ�exc��Ż������ȷ
        s_Pitch = exc(n*FL-222:n*FL);
        PT = findpitch(s_Pitch);    % �����������PT����Ҫ�����գ�
        G = sqrt(E*PT);           % ����ϳɼ���������G����Ҫ�����գ�

        
        % (10) �ڴ�λ��д�������ɺϳɼ��������ü�����filter���������ϳ�����
        %����1.2.2�ģ�8������ͬ����PT���ڵ�֡���ǹ̶��ĳ���
        while pos <= n * FL
            exc_syn(pos) = G;
            pos = pos + PT;
        end
        %exc_syn((n-1)*FL+1:n*FL) =  �������õ��ĺϳɼ���д������
        %���գ�5�����ɼ���exc_syn�õ��ϳ�����s_syn,ͬ��Ҫά���˲���״̬���䣬�õ���zi_syn��zf_syn
        [s_syn1,zf_syn] = filter(1,A,exc_syn((n-1)*FL+1:n*FL),zi_syn);
        zi_syn = zf_syn;
        s_syn((n-1)*FL+1:n*FL) = s_syn1;%�������õ��ĺϳ�����д������

        % (11) ���ı�������ں�Ԥ��ϵ�������ϳɼ����ĳ�������һ��������Ϊfilter
        % ������õ��µĺϳ���������һ���ǲ����ٶȱ����ˣ�������û�б䡣
        %���գ�10������ͬ����FL�����FL_v
        while pos_v <= n * FL_v
            exc_syn_v(pos_v) = G;
            pos_v = pos_v + PT;
        end
        [s_syn1_v,zf_syn_v] = filter(1,A,exc_syn_v((n-1)*FL_v+1:n*FL_v),zi_syn_v);
        zi_syn_v = zf_syn_v;
        s_syn_v((n-1)*FL_v+1:n*FL_v) = s_syn1_v;%�������õ��ļӳ��ϳ�����д������
        % exc_syn_v((n-1)*FL_v+1:n*FL_v) = ... �������õ��ļӳ��ϳɼ���д������
        
        % (13) ���������ڼ�Сһ�룬�������Ƶ������150Hz�����ºϳ�������������ɶ���ܡ�
        while pos_t <= n * FL
            exc_syn_t(pos_t) = G;
            pos_t = pos_t + round(PT/2);
        end
        [r,p,k] = residuez(1,A);%���ԭ����,����ģ����A�ǲ�ַ������ϵ��
        angle = 150 * 2 * pi / 8000;%��
        p = p .* exp(sign(imag(p)) *1i *  angle);%��֪150Hz��Ӧ�ļ���
        [B_t,A_t] = residue(r,p,k);
        [s_syn1_t,zf_syn_t] = filter(1,A_t,exc_syn_t((n-1)*FL+1:n*FL),zi_syn_t);
        zi_syn_t = zf_syn_t;
        s_syn_t((n-1)*FL+1:n*FL) = s_syn1_t;%�������õ��ı���ϳ�����д������
        
    end

    % (6) �ڴ�λ��д������һ�� s ��exc �� s_rec �к����𣬽�����������
    % ��������������ĿҲ������������д���������ر�ע��
    Fs = 8000;%����Ƶ��,��������ÿ 10ms �� 80 ������е�
    sound([exc;s;s_rec]/2^15,Fs,16);
    t = (0:L-1) / Fs; % �ɲ��������ɶ�Ӧʱ��
    figure;
    subplot(3,1,1);plot(t,exc),title('Ԥ��ļ����ź�e(n)'),xlabel("ʱ��/s"),ylabel("���");
    subplot(3,1,2);plot(t,s),title('ԭʼ�����ź�s(n)'),xlabel("ʱ��/s"),ylabel("���");
    subplot(3,1,3);plot(t,s_rec),title('�ؽ��������ź�srec'),xlabel("ʱ��/s"),ylabel("���");
    tpart = (9000:10000) / Fs; % ��ȡһС��
    figure;
    subplot(3,1,1);plot(tpart,exc(9000:10000)),title('Ԥ��ļ����ź�e(n)'),xlabel("ʱ��/s"),ylabel("���");
    subplot(3,1,2);plot(tpart,s(9000:10000)),title('ԭʼ�����ź�s(n)'),xlabel("ʱ��/s"),ylabel("���");
    subplot(3,1,3);plot(tpart,s_rec(9000:10000)),title('�ؽ��������ź�srec'),xlabel("ʱ��/s"),ylabel("���");
    %�����仯�ĵ�λ��ֵ��
    pause(6);
    last_time = 1;
    s1 = unit_sample(200,last_time,Fs);
    s2 = unit_sample(300,last_time,Fs);
    sound([s1;s2],Fs);
    figure;
    subplot(2,1,1);plot(s1);title('200Hz');
    subplot(2,1,2);plot(s2);title('300Hz');
    %������λ��ֵ��
    pause(3);
    last_time = 1;
    s3 = MultiUnitSample(last_time,Fs);
    sound(s3,Fs);
    figure;
    plot(s3);title('MultiUnitSample');
    %��MultiUnitSample��Ϊ����e(n)���뵽ex_1_1.m���˲��������Ϊs(n)
    a1 = 1.3789;a2 = -0.9506;%Ϊ������ֵ
    b = [1, 0, 0];%�����ַ����Ҳ�ϵ��(e(n))
    a = [1, -a1, -a2];%�����ַ������ϵ��(s(n))
    s_multi = filter(b,a,s3);
    pause(3);
    sound([s_multi/max(abs(s_multi));s3],Fs);
    figure;
    subplot(4,1,1);plot(s_multi/max(abs(s_multi)));title('ʱ�����s(n)');
    subplot(4,1,2);plot(s3);title('ʱ������e(n)');
    subplot(4,1,3);fft_plot(s_multi/max(abs(s_multi)),Fs);title('Ƶ�����s(n)');
    subplot(4,1,4);fft_plot(s3,Fs);title('Ƶ������e(n)');
    %�����ϳ�����s_syn��ԭʼ����s������
    pause(3);
    sound([s_syn;s]/2^15,Fs);
    figure;
    subplot(4,1,1);plot(s_syn);title('�ϳ�����ʱ��');
    subplot(4,1,2);plot(s);title('ԭʼ����ʱ��');
    subplot(4,1,3);fft_plot(s_syn,Fs);title('�ϳ�����Ƶ��');
    subplot(4,1,4);fft_plot(s,Fs);title('ԭʼ����Ƶ��');
    %���ٲ����
    pause(3);
    sound(s_syn_v/2^15,Fs);
    figure;
    subplot(4,1,1);plot(s_syn_v);title('���ٲ��������ʱ��');
    subplot(4,1,2);plot(s);title('ԭʼ����ʱ��');
    subplot(4,1,3);fft_plot(s_syn_v,Fs);title('���ٲ��������Ƶ��');
    subplot(4,1,4);fft_plot(s,Fs);title('ԭʼ����Ƶ��');
    %���������
    pause(6);
    sound(s_syn_t/2^15,Fs);
    figure;
    subplot(4,1,1);plot(s_syn_t);title('�������������ʱ��');
    subplot(4,1,2);plot(s);title('ԭʼ����ʱ��');
    subplot(4,1,3);fft_plot(s_syn_t,Fs);title('�������������Ƶ��');
    subplot(4,1,4);fft_plot(s,Fs);title('ԭʼ����Ƶ��');
    % ���������ļ�
    writespeech('exc.pcm',exc);
    writespeech('rec.pcm',s_rec);
    writespeech('exc_syn.pcm',exc_syn);
    writespeech('syn.pcm',s_syn);
    writespeech('exc_syn_t.pcm',exc_syn_t);
    writespeech('syn_t.pcm',s_syn_t);
    writespeech('exc_syn_v.pcm',exc_syn_v);
    writespeech('syn_v.pcm',s_syn_v);
return

% ��PCM�ļ��ж�������
function s = readspeech(filename, L)
    fid = fopen(filename, 'r');
    s = fread(fid, L, 'int16');
    fclose(fid);
return

% д������PCM�ļ���
function writespeech(filename,s)
    fid = fopen(filename,'w');
    fwrite(fid, s, 'int16');
    fclose(fid);
return

% ����һ�������Ļ������ڣ���Ҫ������
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