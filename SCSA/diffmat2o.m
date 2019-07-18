    %**********************************************************************
    %*********             second-order differentiation matrix           *********
    %**********************************************************************
    
    
%Author: Zineb Kaisserli
% modified by Gabriel Benigno to be consistent with notation in SCSA paper (Math. Control Signals Syst. (2013) 25:37?61);
% also got rid of Fs argument; made =1 always

function D2=diffmat2o(y, Delta_x)
    M = length(y);
    Delta = 2*pi/M;
    ex = kron( (M-1):-1:1, ones(M,1) ); % 'j-k' from eqns 35 and 36 of article
    if mod(M,2)==0
        dx = -pi^2/(3*Delta^2)-(1/6)*ones(M,1);
        test_bx = -(-1).^ex*(0.5)./(sin(ex*Delta*0.5).^2);
        test_tx =  -(-1).^(-ex)*(0.5)./(sin((-ex)*Delta*0.5).^2);
    else
        dx = -pi^2/(3*Delta^2)-(1/12)*ones(M,1);
        test_bx = -0.5*((-1).^ex).*cot(ex*Delta*0.5)./(sin(ex*Delta*0.5));
        test_tx = -0.5*((-1).^(-ex)).*cot((-ex)*Delta*0.5)./(sin((-ex)*Delta*0.5));
    end
    Ex = full(spdiags([test_bx dx test_tx],[-(M-1):0 (M-1):-1:1],M,M));
    
    D2=(Delta/Delta_x)^2*Ex;
    
end