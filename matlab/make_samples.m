function [error, varargout] = make_samples(f, nsamples, varargin)
out_list = [];
error = 0;
mu = 0;
varargout{1} = 0;
if nargout == 3
    varargout{2} = 0;
end
for i=1:nsamples
    if nargout == 3
        [a, b] = f(varargin{:});
        out_list = [out_list, b];
        mu = mu + b / nsamples;
        varargout{1} = varargout{1} + a / nsamples;
        varargout{2} = mu;
    else
        [a] = f(varargin{:});
        out_list = [out_list, a];
        mu = mu + a / nsamples;
        varargout{1} = mu;
    end
end
if nsamples > 1
    error = sqrt(1/(nsamples-1) * sum((out_list-mu).^2)) / mu;
end
end