function freq_gen = measure_gen(tide,freq,tide_query)

freq_gen = interp1(tide,mean(freq),tide_query)...
    +mean(std(freq))*randn(size(tide_query));
end