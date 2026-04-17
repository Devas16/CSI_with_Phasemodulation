function Xret=compute_phase(Distmat_2D,phase_mod,rgrid,indx_valid)
rgrid=rgrid(:);
phase_mod=phase_mod(:);
rm=min(rgrid);
rM=max(rgrid);
N=length(rgrid);
dvec=Distmat_2D(:);
dvec=dvec(indx_valid);

kvec=1+((dvec-rm)*(N-1)/(rM-rm));

% Interpolate for the points in the index set
kleft=floor(kvec);
kright=ceil(kvec);

% kfrac=(dvec-rgrid(kleft))./(rgrid(kright)-rgrid(kleft));

xvec=zeros(size(Distmat_2D(:)));
xvec(indx_valid)=phase_mod(kleft);%((1-kfrac).*phase_mod(kleft))+(kfrac.*phase_mod(kright));%
Xret=reshape(xvec,size(Distmat_2D)); 
end