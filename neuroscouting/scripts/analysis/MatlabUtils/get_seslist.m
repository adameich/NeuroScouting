function [ sessArray ] = get_seslist( subDir, sub )
%UNTITLED Summary of this function goes here
%   Give it a path where all your subjects' data folder are and it will
%   return a cell with all the subject numbers
preSess = dir([subDir char(sub)]);
preSessCell = {preSess.name};
tmp = 1;
for i = 1:size(preSessCell,2)
	if size(preSessCell{i}, 2)==7;
		sessArray(tmp) = preSessCell(i);
		tmp = tmp + 1;
	end
end

end

