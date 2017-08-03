function [ImTrans,tau] = preAlign(ImData,Imref,tau)
% for pre-alignment
numFrame = size(ImData,3);
ImTrans = ImData;
if nargin<3
tau = zeros(6,numFrame);
end
numLevel = fix(log(size(ImData,1)*size(ImData,2))/log(2)/2-3);
for i = 1:numFrame
    [ImTrans(:,:,i),tau(:,i)] = regMGNC(Imref,ImTrans(:,:,i),tau(:,i),numLevel,2);
end
end
