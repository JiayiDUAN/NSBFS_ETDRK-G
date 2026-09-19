function averxA=averx_m(A)
averxA=(A+A(:,[end 1:end-1]))/2;
end