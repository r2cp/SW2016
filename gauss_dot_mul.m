function y = gauss_dot_mul(a,b)
%Implements dot multiplication (.*) in Gauss mode
% by Timos Papadopoulos 
% Both inputs a and b have to be 2D numeric or logical arrays.
% 
% The output y is equal to y = a.*b = b.*a in the manner of Gauss dot multiplication
% 
% y = gauss_dot_mul(a,b)


validateattributes(a,{'double','logical'},{'2d'});
validateattributes(b,{'double','logical'},{'2d'});

if      isequal(size(a),size(b))                                        % both inputs NxK
    y = a.*b;
elseif  isscalar(a) || isscalar (b)                                     % one of two inputs scalar (excluding previous cases)
    y = repmat(a,size(b)).*repmat(b,size(a));
elseif  (size(a,2) == size(b,2)) && (size(a,1) == 1 || size(b,1) == 1)	% one input 1xN the other KxN (excluding previous cases)
    y = repmat(a,size(b,1),1).*repmat(b,size(a,1),1);
elseif  (size(a,1) == size(b,1)) && (size(a,2) == 1 || size(b,2) == 1)  % one input Nx1 the other NxK (excluding previous cases)
    y = repmat(a,1,size(b,2)).*repmat(b,1,size(a,2));
else
    error('gauss_dot_mul:incompatible_input_dims','This dot multiplication is not valid in Gauss')
end
