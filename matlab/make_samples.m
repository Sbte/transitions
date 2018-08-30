function [data, varargout] = make_samples(f, nsamples, varargin)
data.tp = [];

data.sigma = 0;
data.mu = 0;
data.median = 0;

data.tp = [];
data.mfpt = [];

data.Q1 = 0;
data.Q3 = 0;

varargout{1} = 0;
if nargout == 3
    varargout{2} = 0;
end
for i=1:nsamples
    if nargout == 3
        [a, b] = f(varargin{:});
        data.tp = [data.tp, a];
        data.mfpt = [data.mfpt, b];
        data.mu = data.mu + b / nsamples;
        varargout{1} = varargout{1} + a / nsamples;
        varargout{2} = data.mu;
    else
        [a] = f(varargin{:});
        data.tp = [data.tp, a];
        data.mu = data.mu + a / nsamples;
        varargout{1} = data.mu;
    end
end

if nsamples > 1
    if nargout == 3
        d = sort(data.mfpt);
    else
        d = sort(data.tp);
    end
    data.sigma = sqrt(1/(nsamples-1) * sum((d-data.mu).^2));
    data.median = compute_median(d);
    data.Q1 = compute_median(d(1:floor(length(d)/2)));
    data.Q3 = compute_median(d(ceil(length(d)/2):end));
end

end

function out = compute_median(d)
    if mod(length(d), 2) == 0
        out = (d(length(d) / 2) + d(length(d) / 2 + 1)) / 2;
    else
        out = d(ceil(length(d) / 2));
    end
end