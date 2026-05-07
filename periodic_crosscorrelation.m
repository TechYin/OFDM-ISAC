function [Rk_right,Rk_left,Rk_both] = periodic_crosscorrelation(x,y)
    % 计算离散周期信号的周期互相关函数（PCCF）
    %
    % 参数:
    %   x : 输入信号（支持实数和复数信号），长度为 N
    %
    % 返回:
    %   Rk_right: 循环移位时，将取共轭的拷贝信号向右依次进行循环移位，得到的PACF函数（Rk）,长度为 N, 对应位移 k=0,1,...,N-1
    %   Rk_left: 循环移位时，将取共轭的拷贝信号向左依次进行循环移位，得到的PACF函数（R-k）,长度为 N, 对应位移 k=0,1,...,N-1
    %   关系：Rk_right = conj(Rk_right)
    %   Rk_both : 双边周期自相关函数，长度为 2N+1，对应位移 k=0,1,...,N-1

    x = x(:); 
    y = y(:);
%     N = length(x);

    % 使用FFT加速计算（维纳辛钦定理：周期信号的自相关函数和功率谱密度呈傅里叶变换关系）
    X = fft(x);                     % 计算FFT
    Y = fft(y);
    power_spectrum = X .* conj(Y);  % 计算互功率谱密度
    Rk_right = ifft(power_spectrum);       % 逆FFT
%     Rk_right = real(R);                    % 取实部（消除数值误差）
%     R = circshift(R, 1);            % 将位移 k=0 对齐到第一个元素
    Rk_right = Rk_right / sqrt(sum(abs(x.^2))*sum(abs(y.^2)));                      
                                      % 能量归一化，得到相关系数
    Rk_left = conj(Rk_right);
    Rk_both = [flip(Rk_left);Rk_right(2:end)];  
        % 由于Rk_left和Rk_right中都有k=0的重复项，需要略去其中一项
        % Rk_left需要翻转一下
end