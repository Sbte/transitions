function [data, varargout] = make_samples(f, nsamples, varargin)
data.tp = [];

data.sigma = 0;
data.mu = 0;
data.median = 0;
data.avg_steps = 0;

data.tp = [];
data.mfpt = [];
data.steps = [];

data.Q1 = 0;
data.Q3 = 0;

varargout{1} = 0;
varargout{2} = 0;
if nargout == 3
    varargout{3} = 0;
end
for i=1:nsamples
    % fs = functions(f);
    % fprintf('%s: sample %d\n', fs.function, i);
    if nargout == 3
        [a, b, c] = f(varargin{:});
        data.tp = [data.tp, a];
        data.mfpt = [data.mfpt, c];
        data.steps = [data.steps, b];
        data.avg_steps = data.avg_steps + b / nsamples;
        data.mu = data.mu + c / nsamples;
        varargout{1} = varargout{1} + a / nsamples;
        varargout{2} = data.mu;
    else
        [a, b] = f(varargin{:});
        data.tp = [data.tp, a];
        data.steps = [data.steps, b];
        data.avg_steps = data.avg_steps + b / nsamples;
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
    data.normalized_error = data.sigma / data.mu * sqrt(data.avg_steps);
end

end

function out = compute_median(d)
    if mod(length(d), 2) == 0
        out = (d(length(d) / 2) + d(length(d) / 2 + 1)) / 2;
    else
        out = d(ceil(length(d) / 2));
    end
end