% measure how symmetric a matrix is
% 1 = perfectly symmetric
% -1 = perfectly asymmetric
% note that indifference point is not 0! need to do permutation / random test
%
% from https://math.stackexchange.com/questions/2048817/metric-for-how-symmetric-a-matrix-is
function s = how_sym(A)

    A_sym = 0.5 * (A + A');
    A_anti = 0.5 * (A - A');

    s = (norm(A_sym) - norm(A_anti)) / (norm(A_sym) + norm(A_anti));
