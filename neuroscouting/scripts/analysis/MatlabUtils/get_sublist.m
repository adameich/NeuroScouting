function [ subArray ] = get_sublist( subDir )
%UNTITLED Summary of this function goes here
%   Give it a path where all your subjects' data folder are and it will
%   return a cell with all the subject numbers
preSubs = dir(subDir);
preSubsCell = {preSubs.name};
tmp = 1;
for i = 1:size(preSubsCell,2)
	if size(preSubsCell{i}, 2)==7;
		subArray(tmp) = preSubsCell(i);
		tmp = tmp + 1;
	end
end

end

