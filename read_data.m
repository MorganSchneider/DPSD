%clear all

load('koun_data_orig.mat')
fEL_desired = 0.5; %deg
pulses_per_rad = 64;
EL_persist = 100;
ptot = size(pulse,2);
cut_started = 0;
cut_ended = 0;
pind1 = 0;
pind2 = 0;
rind = 0;
while (~cut_started || (cut_started && ~cut_ended)) && pind1 < ptot
    pind1 = pind1 + 1;
    if ~cut_started
        if all(abs(header.fEL(pind1:min([pind1+EL_persist,ptot])) - fEL_desired) < 0.05)
            cut_started = 1;
        end
    end
    if cut_started && ~cut_ended
        if sum(abs(header.fEL(pind1:min([pind1+EL_persist,ptot])) - fEL_desired) > 0.05) > EL_persist/2
            cut_ended = 1;
        end
    end
    if cut_started && ~cut_ended
        if pind2 == 0
            rind = rind + 1;
            az(rind) = header.fAZ(pind1 + round(pulses_per_rad/2));
        end
        pind2 = pind2 + 1;
        % Vh and Vv dimensions are range_gate_no, radial_no, and pulse_no
        % (r, az, m)
        iqh(:, rind, pind2) = pulse(1, pind1).I + 1i*pulse(1, pind1).Q;
        iqv(:, rind, pind2) = pulse(2, pind1).I + 1i*pulse(2, pind1).Q;
        if pind2 == pulses_per_rad
            pind2 = 0;
        end
    end
    if cut_ended && pind2 ~= 0
        iqh(:, rind, :) = [];
        iqv(:, rind, :) = [];
        az(rind) = [];
        rind = rind - 1;
    end
end




       
    
        