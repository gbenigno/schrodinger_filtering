function T_vec = Tvec(T,t)

% author: Gabriel Benigno

% given Ton_shift or Toff_shift only, make binary vector with zeros padded between gaps

%%% inputs:
% T: Ton_shift or Toff_shift...ie, Ton or Toff with indices starting at 1
% t: entire time vector (also shifted so that index starts at 1)

% output: T_vec: vector of T values as ones and zeros padded between gaps

T_vec = zeros(length(t),1);
T_vec(T) = 1;

end

