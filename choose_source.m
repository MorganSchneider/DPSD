function [iq_source, signal_type] = choose_source()

[st,~] = dbstack('-completenames',1);

if strcmp(st.name, 'dpsd_bootstrap') || strcmp(st.name, 'dpsd_bootstrap_v2') || strcmp(st.name, 'dpsd_plotter')
    x = input(['Choose I/Q data to use: ', newline,...
        '  (1) iq_emulate_yu', newline,...
        '  (2) SimRadar', newline]);
    y = input(['Choose simulation type: ', newline,...
        '  (1) Rain only', newline,...
        '  (2) Debris only', newline,...
        '  (3) Rain and debris', newline]);
    
    if x == 1
        iq_source = 'yu';
    elseif x == 2
        iq_source = 'simradar';
    else
        disp('Try again, asshole!')
        choose_source
    end
    
    if y == 1
        signal_type = 'rain';
    elseif y == 2
        signal_type = 'debris';
    elseif y == 3
        signal_type = 'multi';
    else
        disp('Try again, asshole!')
        choose_source
    end
end