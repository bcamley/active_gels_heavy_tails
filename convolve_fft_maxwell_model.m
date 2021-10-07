function u = convolve_fft_maxwell_model(t,F,params)
%beta = params.beta;

assert(~contains(params.dipole_resp_func,'saffman','IgnoreCase',true),'Simple convolution is not correct for Saffman-Delbruck tensor')
N = length(t);
T = t(end)-t(1);

om = (2*pi/T)*[0:N/2-1 0 -N/2+1:-1]; 

%try
Kom = (1-1i*params.omega0./om)./(params.Gscale);
%catch
%        Kom = 1./(params.G0 + (1i*om/1).^beta);
%        warning('Assuming omega0, Gscale = 1');
%end
Kom(abs(om)<100*eps) = 0; % remove the zero mode -- this will only be important
                          % if G_0 = 0
                          
u = real(ifft(Kom.*fft(F)));