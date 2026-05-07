clc;
clear
close all
%% Calculation correlation function between different pilot-pattern detection signals.
load PilotPatterns.mat N K SymbolInPilot PilotPattern_Dictionaries DFT
Dictionary_index = 1; % Totally 5 dictionaries, Dictionary_index is an integer from 1 to 5
pop = PilotPattern_Dictionaries{Dictionary_index};


M = size(pop,2);
for i = 1:M
    F_P =  DFT(:,pop(:,i));
    P = F_P*diag(SymbolInPilot);


    % PD信号自相关ISL、PSL计算
    X_PPDS = zeros(1,K);
    X_PPDS(:,pop(:,i)) = SymbolInPilot;        % PIDS:Pilot-Patterns Detection Signal
    SignalAftNull = X_PPDS.';
    SignalIFFT = ifft(fftshift(SignalAftNull,1),K,1) * sqrt(K);  % ifft(Y,n,1) 返回每一列的 n 点逆变换。
    % 计算滑动相关
    PilotDetectSignal(i,:) = SignalIFFT;  
end

for m = 1:M
    for j = 1:M
        [PilotDetect_corr,~,~] = periodic_crosscorrelation(PilotDetectSignal(m,:),PilotDetectSignal(j,:));
        PCC_PL_single(m,j) = max(abs(PilotDetect_corr));
    end
end
PCC_Opt = [];
for x = 1:size(PCC_PL_single,1)
    for y = 1:size(PCC_PL_single,2)
        if x~=y
            PCC_Opt = [PCC_Opt PCC_PL_single(x,y)];
        end
    end
end



figure;
set(gcf, 'Position', [1100, 450, 350, 250]);  
PCC_show = PCC_PL_single(12:19,12:19);
imagesc(PCC_show); 
set(gca, 'YDir', 'normal'); % 翻转Y轴
set(gca,'FontSize',10,'Fontname', 'Times New Roman')
% colorbar; 
for x = 1:size(PCC_show,1)
    for y = 1:size(PCC_show,2)
        if x==y
            text(x, y, num2str(PCC_show(x,y),'%.3f'), 'FontSize', 9, 'FontName','Times New Roman','HorizontalAlignment','center','Color','Black'); % x, y 是文本的位置
        else
            text(x, y, num2str(PCC_show(x,y),'%.3f'), 'FontSize', 9, 'FontName','Times New Roman','HorizontalAlignment','center','Color','White'); % x, y 是文本的位置
        end
    end
end
xlabel('Index','Fontsize',12); ylabel('Index','Fontsize',12);

