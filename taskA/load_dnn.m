function dnn = load_dnn(folder)
% LOAD_DNN  Load DNN weights from .npy files into a struct.
%
%   dnn = load_dnn()          % looks in current directory
%   dnn = load_dnn(folder)    % looks in the given folder
%
%   The function expects six files:
%     W1.npy, b1.npy, W2.npy, b2.npy, W3.npy, b3.npy
%
%   Output struct fields:
%     dnn.W1, dnn.b1, dnn.W2, dnn.b2, dnn.W3, dnn.b3

if nargin < 1 || isempty(folder)
    folder = '';
end

pjoin = @(f) fullfile(folder, f);

dnn.W1 = readNPY(pjoin('W1.npy'));
dnn.b1 = readNPY(pjoin('b1.npy')); dnn.b1 = dnn.b1(:)';
dnn.W2 = readNPY(pjoin('W2.npy'));
dnn.b2 = readNPY(pjoin('b2.npy')); dnn.b2 = dnn.b2(:)';
dnn.W3 = readNPY(pjoin('W3.npy'));
dnn.b3 = readNPY(pjoin('b3.npy')); dnn.b3 = dnn.b3(:)';

fprintf('DNN loaded | W1:%s  W2:%s  W3:%s\n', ...
    mat2str(size(dnn.W1)), mat2str(size(dnn.W2)), mat2str(size(dnn.W3)));

end
