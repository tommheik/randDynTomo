function x = dec2waverec3(d,wSys)
%DEC2WAVEREC3 uses the decomposition coefficients d and the wavelet system
%wSys to reconstruct the output x using wavrec3
%
% T H   2023
wSys.dec = d;
x = waverec3(wSys);
end

