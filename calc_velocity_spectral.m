function vr_new = calc_velocity_spectral(psd_filt, M, Ts, lambda)

km = find(psd_filt == max(psd_filt)) - M/2;
k = km - M/2 + 1: km + M/2;
P = sum(psd_filt);

sum_term = 0;
for i = 1:length(k)
    sum_term = sum_term + ((k(i) - km) * psd_filt(mod(k(i)+M/2,M)+1));
end

vr_new = -lambda/(2*M) * (km/Ts + (1/(P*Ts) * sum_term));

end