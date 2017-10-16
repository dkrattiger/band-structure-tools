function DependencyCopy(myfunction)

% myfunction is a string (i.e. myfunction = "myfun.m") containing the
% function in question

% get a list of 
[file_list,~] = matlab.codetools.requiredFilesAndProducts(myfunction);

for i = 1:length(file_list)
    status = copyfile(file_list{i});
end