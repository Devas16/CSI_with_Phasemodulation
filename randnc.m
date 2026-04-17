function out = randnc(m,n)
out = randn(m,n)+1j*randn(m,n);
out = out/sqrt(2);
end