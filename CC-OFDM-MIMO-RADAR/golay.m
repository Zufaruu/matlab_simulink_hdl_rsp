function g_mat = golay(N)
g(1)=1;
g_com(1)=1;
% Ncb=256;
n=log2(N);
for i_mat=1:n
    g_mat = [g g_com;g -g_com];
    g_com_mat = [g -g_com;g g_com];
    g = g_mat;
    g_com = g_com_mat;
end

