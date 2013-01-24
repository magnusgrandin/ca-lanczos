%% function [y, idx] = leja (x, which)
%
%   Given a vector of complex numbers x, compute a Leja ordering of x:
%   return a vector y and a permutation array idx such that x(idx) = y
%   is the Leja ordering of x.
% 
%   If the second argument is 'modified', compute the "modified" Leja
%   ordering instead.  This assumes that complex conjugate pairs of x 
%   are both adjacent, and ordered so that the element of the pair with
%   positive imaginary part occurs first.
% 
%   Note that Matlab (via LAPACK's DHSEQR) already returns the eigenvalues
%   of a real matrix so that complex conjugate pairs occur together,
%   ordered so that the element with positive real part occurs first.
%  
%   Note: the Leja ordering of a set of points is not unique:  see e.g.
%   Baglama, Calvetti and Reichel:  "Fast Leja points", ETNA, Vol. 7,
%   1998, pp. 124-140.
% 
%   I also think this was broken before 27 Sep 2009.  Before then, I only
%   had the non-modified Leja ordering anyway...

function [y, idx] = leja (x, which)
    if (nargin < 2)
        % FIXME (mfh 01 Oct 2009) should count the actual multiplicities,
        % rather than assuming all the given shifts are unique...
        [y, idx] = nonmodified_leja(x, length(x), ones(1,length(x)));
    else
        [y, idx] = real_leja(x);
    end
end

% -------------------------------------------------------------------------
%  Contribution by Erin Carson, University of California at Berkeley
% -------------------------------------------------------------------------
