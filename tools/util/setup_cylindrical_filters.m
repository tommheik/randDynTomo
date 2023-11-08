function [shear_fe]=setup_cylindrical_filters(BP_sizes,decomp,dsize,level)
% The function sets up the cylinder decompostion filters for data of size
% xSz using decompositions with sizes dsize up to level specified.
%  Input: BP_sizes - Laplace pyramid sizes
%         decomp - set of exponents of 2 that determing shearing filters
%         dsize - set of dyatic lengths determing the sizes of the shearling filters 
%         level - the levels of shearlet decomp. contained in dst 
%
%  Output: shear_fe - 3D cylinder shearing filters 
%
% Original function 'setup_cylindrical_shear.m'
% Written by Glenn Easley on October 12, 2020.
% Adapted by T. Heikkilä  2022
%


shear_fe=cell(1,level); % declare cell array containing shearing filters

% reverse order of decomp and dsize to match 3D LP transform
decomp=fliplr(decomp);dsize=fliplr(dsize);

for i=1:level
     L0(3)= BP_sizes{i}(3);
    shear_f=shearing_filters_Myer(dsize(i),decomp(i)).*sqrt(dsize(i));
    %%%%%%%%%%%%%%%
    % replicate the shearing filter to be of the same z-dimension as BP{i}
    %%%%%%%%%%%%%%%
    for k=1:2^decomp(i)
        for ss0=1:L0(3)
            shear_fe{i}(:,:,ss0,k)=shear_f(:,:,k);
        end
    end
    %%%%%%%%%%%%%%%
end


