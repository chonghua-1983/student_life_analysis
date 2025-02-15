function S = ochiia(coMat)
   D = diag(coMat);
   D = sqrt(D);
   dec = D*D';
   S = coMat./dec;
end