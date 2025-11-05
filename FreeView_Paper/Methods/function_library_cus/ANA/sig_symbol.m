function s = sig_symbol(p)
% 将 p 值转为显著性符号
if p < 1e-3
    s = '***';
elseif p < 1e-2
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = 'n.s.';
end
end