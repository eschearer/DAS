function AB = mult(A,B)
% multiply a 3D matrix A by a column vector B
    n = size(A,3);            % size of 3rd dimension of A
    m = size(B,1);            % size of B
    
    if (n~=m)
        error('Third dimension of first matrix must be the same as the length of the second matrix (column vector).');
    end
    
    AB = zeros(size(A,1), size(A,2));
    for i=1:n
        AB = AB + A(:,:,i) * B(i);
   end
end