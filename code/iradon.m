function f = iradon(R)

% iradon - Inv_FastSlantStack in quiet mode (intercept stdout)
%
%   f = iradon(R);
%
%   Copyright (c) 2015 Gabriel Peyre

[~,f] = evalc('real( Inv_FastSlantStack(R) )');

end