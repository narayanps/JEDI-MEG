function [source_inds] = select_source_inds(mesh, num_sources, dist_thr)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
msize = length(mesh.p);

distMat = inf*eye(num_sources,num_sources);
dmin = 0;
while dmin < dist_thr % just to ensure the dipoles drawn are not too close (in mm )
    idx = randperm(msize);
    source_inds = idx(1:num_sources);
    for i=1:1:num_sources
        d= mesh.p(source_inds,:) - repmat(mesh.p(source_inds(i),:), num_sources,1);
    for j=1:1:num_sources
        if i ~= j
        distMat(i,j) = norm(d(j,:));
        end
    end
    end
    dmin = min(distMat(:));
end
end

