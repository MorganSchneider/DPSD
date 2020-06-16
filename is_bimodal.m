function bimodal_inds = is_bimodal(x,dim)
% x = Sxx_mat;
% dim = 1; % dimension that time series is oriented along

ndims = length(size(x));
ndims(ndims == dim) = [];
xx = permute(x, [dim ndims]);
xx = reshape(xx, [size(x,dim), numel(x)/size(x,dim)]);

xmins = repmat(min(abs(xx),[],1), [size(x,dim) 1]);
fx = round(10*log10(abs(xx) ./ xmins));

p_vals = zeros(1,numel(x)/size(x,dim));
bimodal_inds = p_vals;
for i = 1:size(xx,2)
    ffx = fx(:,i);
    vec = [];
    for j = 1:length(ffx)
        vec = [vec, j*ones(1,ffx(j))];
    end
    [~, p_value, ~, ~] = HartigansDipSignifTest(vec, 500);
    p_vals(i) = p_value;
end

bimodal_inds(p_vals < 0.01) = 1;

sz = size(x);
sz(dim) = [];
if length(sz) == 1
    sz = [1 sz];
end
bimodal_inds = reshape(bimodal_inds, sz);


end