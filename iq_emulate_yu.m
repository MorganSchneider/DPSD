% Generates IQ data based on Yu et al. (2012)

function [V_h, V_v, ZDR, PHV] = iq_emulate_yu(Sv, sphv, szdr, K)

    M = max(size(Sv));
    debug_flag = true;
    
    Zv = zeros(K,M);
    
    for k = 1:K
        rand_u = rand(1,M);
        rand_phi = rand(1,M) * 2 * pi - pi;
        Zv(k,:) = sqrt(-Sv .* log(rand_u)) .* exp(1i .* rand_phi);
        
        % Requires second independent realization of Zv
        rand_u = rand(1,M);
        rand_phi = rand(1,M) * 2 * pi - pi;
        Zv2 = sqrt(-Sv .* log(rand_u)) .* exp(1i .* rand_phi);
        
        Zh(k,:) = sqrt(szdr) .* (sphv .* Zv(k,:) + sqrt(1-sphv.^2) .* Zv2);
        V_h(k,:) = ifft(fftshift(Zh(k,:)), M);
        V_v(k,:) = ifft(fftshift(Zv(k,:)), M);
        
        SHH(k) = mean(V_h(k,:) .* conj(V_h(k,:)));
        SVV(k) = mean(V_v(k,:) .* conj(V_v(k,:)));
        SXX(k) = mean(V_h(k,:) .* conj(V_v(k,:)));
        ZDR(k) = 10 * log10(SHH(k) ./ SVV(k));
        PHV(k) = abs(SXX(k)) ./ sqrt(SHH(k) .* SVV(k));
    end

    if debug_flag
        figure()
        ind = ceil(rand(1) * K);
        plot(1:M, real(V_h(ind,:)), '-k', 1:M, imag(V_h(ind,:)), '--k',...
            1:M, real(V_v(ind,:)), '-b', 1:M, imag(V_v(ind,:)), '--b');
        title(['ZDR=', num2str(ZDR(ind)), '; RHV=', num2str(PHV(ind))]);
    end
    
    V_h = V_h';
    V_v = V_v';
end