function [entropy, p_i] = presence_entropy(data, base)
%Degree Structure (Presence) Entropy
%   Calculates entropy based on the network's strength distribution. Exact
%   formulation from Wen & Jiang (2019), refferencing Wang et al., 2006.
%   Input:
%       data: a Time X StateStrength matrix OR a vector of State Strengths.
%       base: either 'bits' for log2 or 'nats' for natural log (default).
%   Output:
%       entropy: A Time X 1 matrix of the presence (or "strength") entropy 
%           of the input data.
%       p_i: A matrix of the size(data) with probabilities for entropy
%           calculation.


if iscolumn(data)
    data = data';
end

if nargin<2
    base = 'nats';
end

data = abs(data);

n_time = size(data,1);
n_states = size(data,2);

entropy = zeros(n_time,1); %initialize output
p_i = zeros(size(data));

for t = 1:n_time %loop through time

    s_vec = zeros(n_states,1);
    s_sum = sum(data(t,:));

    for i = 1:n_states %loop through space

        p_i(t,i) = data(t, i)/s_sum; %presence of state i

        if p_i(t,i) == 0
            continue;
        end

        if strcmp(base,'nats')
            s_vec(i) = p_i(t,i)*log(p_i(t,i)); % plog(p)
        elseif strcmp(base,'bits')
            s_vec(i) = p_i(t,i)*log2(p_i(t,i)); % plog(p)
        end

    end

    entropy(t) = -1 .* sum(s_vec);
    
end

end

