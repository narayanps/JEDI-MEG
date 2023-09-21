function [ind_, r_, G_] = set_part_traj(pnts, G, dist_thr, num_sources, T)

n_pnts=size(pnts,1);
ind_ = zeros(num_sources, T);
distMat = inf*eye(num_sources,num_sources);
dmin = 0;
nr = 3*num_sources;
%

while dmin < dist_thr
      idx = randperm(n_pnts);
      id = idx(1:num_sources);

    for i=1:1:num_sources
        d= pnts(id, :) - repmat(pnts(id(i), :), num_sources,1);
        for j=1:1:num_sources
            if i ~= j
                distMat(i,j) = norm(d(j,:));
            end
        end
    end
    dmin = min(distMat(:));
end

for j=1:1:num_sources
    ind_(j,:) = repmat(id(j), [1 T]);
end

r_ = repmat(reshape(pnts(ind_(:,1),:)',nr,1),[1 T]);
G_ =repmat(G(:, ind_(:,1)),[1 1 T]);


end