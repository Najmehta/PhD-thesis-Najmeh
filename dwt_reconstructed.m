function [wt_recon]=dwt_reconstructed(newwc, wl, wfilter, scale);
% The function is used to reconstruct the signals with the selected DWT coefficients
%[wt_recon]=dwt_reconstructed(newwc, wl, wfilter, scale);
%
%syntax: [wt_recon]=dwt_reconstructed(newwc, wl, wfilter, scale);
%
%Input
%newwc:  the selected DWT coefficient Matrix in row vectors, the unselected coefficients were assigned to 0,
%              newwc(i,:)=[C(N) D(N) D(N-1) ... D(1)];
% wl:        the length of discrete wavelet coeffcients corresponding to each scale
%              wl=[length(C(N)), length(D(N)),..., length(D(1))];
%wfilter:   the wavelet filter used for DWT
%scale:     the scale used for DWT
%
%Output
%wp_recon: the reconstructed signals
%
%Edited by Da Chen, Dec 18, 2008
%version 1.0
%Modified by Da Chen, Feb 20, 2009

[calrow,calcol]=size(newwc);

for i=1:calrow
    wt_recon(i,:)=waverec(newwc(i,:), wl, wfilter);
end