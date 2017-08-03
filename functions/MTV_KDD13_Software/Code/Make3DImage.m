function image = Make3DImage(n,K)

image = zeros(n,n,K);
I=MakeRDSquares(n,7,40);
I(I~=0) = 1;
for k = 1:fix(K/3)
    image(:,:,k) = I;
end

I=MakeRDSquares(n,7,40);
I(I~=0) = 1;
for k = fix(K/3)+1:fix(K*2/3)
    image(:,:,k) = I;
end

I=MakeRDSquares(n,7,40);
I(I~=0) = 1;
for k = fix(K*2/3)+1:K
    image(:,:,k) = I;
end