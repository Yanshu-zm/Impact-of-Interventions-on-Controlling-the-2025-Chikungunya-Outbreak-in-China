function success = createDirIfNotExist(dirPath)
% CREATEDIRIFNOTEXIST 检查目录是否存在，不存在则创建
%   success = createDirIfNotExist(dirPath)
%   输入:
%     dirPath - 要检查或创建的目录路径
%   输出:
%     success - 逻辑值，表示操作是否成功

% 检查目录是否已存在
if ~exist(dirPath, 'dir')
    try
        % 尝试创建目录
        success = mkdir(dirPath);
        if success
            fprintf('目录 "%s" 创建成功\n', dirPath);
        else
            fprintf('创建目录 "%s" 失败\n', dirPath);
        end
    catch ME
        % 捕获并显示异常信息
        fprintf('创建目录时出错: %s\n', ME.message);
        success = false;
    end
else
    % 目录已存在
    fprintf('目录 "%s" 已存在\n', dirPath);
    success = true;
end