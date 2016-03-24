function img = filterDoG(img, n1, n2)
n1 = ceil(n1);
n2 = ceil(n2);
if n1~=0,
    k1 = hamming(n1*2+1);
    k1 = k1/sum(k1);
    img = imfilter(imfilter(img, k1, 'symmetric'), k1', 'symmetric');
end
if n2 ~=0;
    k2 = hamming(n2*2+1);
    k2 = k2/sum(k2);
    img = img - imfilter(imfilter(img, k2, 'symmetric'), k2', 'symmetric');
end
