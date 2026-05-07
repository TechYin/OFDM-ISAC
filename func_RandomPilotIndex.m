%% N个元素的数组里包含x个1
function [Position,PilotIndex,DataIndex] = func_RandomPilotIndex(PilotNum,AllNum)
Position = zeros(AllNum,1);
ind = randperm(AllNum, PilotNum);   % 给出x个不重复的随机值作为索引
Position(ind) = 1;
PilotIndex = find(Position==1);
DataIndex = find(Position==0);

end