function l_ei = lambda_ei(n,T_e)

%log_factor = 22.7; % not sure where this number came from
log_factor = 23.0; % from NRL and Fitzpatrick
%log_factor = 23.4633; % In order to agree with Chen


if length(n) == 1 && length(T_e) == 1

    if (sqrt(n)*T_e^(-3/2)) > 1
        l_ei = log_factor - log(sqrt(n)*T_e^(-3/2));
    else
        l_ei = log_factor;
    end
    
elseif length(n) == 1
    
    l_ei = log_factor*ones(size(T_e));
    bool = (sqrt(n)*T_e.^(-3/2)) > 1;
    l_ei(bool) = log_factor - log(sqrt(n)*T_e(bool).^(-3/2));
    
elseif length(T_e) == 1
    
    l_ei = log_factor*ones(size(n));
    bool = (sqrt(n)*T_e^(-3/2)) > 1;
    l_ei(bool) = log_factor - log(sqrt(n(bool)*T_e^(-3/2)));
    
else
    
    l_ei = log_factor*ones(size(n));
    bool = (sqrt(n).*T_e.^(-3/2)) > 1;
    
    l_ei(bool) = log_factor - log(sqrt(n(bool)).*T_e(bool).^(-3/2));
    
end