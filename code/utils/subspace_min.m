% Computes the minimum angle between two subspaces. This function is a
% modification of the matlab function SUBSPACE, which returns the largest
% angle. SUBSPACE_MIN returns the smallest one, and can also return the
% aligments for all dimensions, sorted by decreasing angle.

%SUBSPACE Angle between subspaces.
%   SUBSPACE(A,B) finds the angle between two subspaces specified by the
%   columns of A and B.

function [theta, theta_all] = subspace_min(A,B)

A = orth(A);
B = orth(B);
%Check rank and swap
if size(A,2) < size(B,2)
   [A,B] = swap(A,B); 
end
% Compute the projection according to [1].
% If a and b are vectors, this equation substracts to b the projection of b
% into a. Geometrically, this results in a vector going from the tip of b
% to the tip of the a projection, the vector difference.
B = B - A*(A'*B);

% Make sure it's magnitude is less than 1.
% Here we apply min(svd(B)), instead of norm(B)=max(svd(B)) as in the
% original subspace function. This is because the maximum singular value of
% B gives the largest distance between B and A (for vectors a and b, the
% largest vector difference) By taking the min we obtain the minimum angle.
theta = asin(min(ones(superiorfloat(A,B)),min(svd(B))));

% compute all angles
theta_all = asin(min(ones(superiorfloat(A,B)),svd(B)));

function [B, A] = swap(A, B)
end

end