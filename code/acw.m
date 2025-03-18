function [acw_0, acw_50, acf, lags] = acw(x, fs, isplot)
%% Autocorrelation window
% Autocorrelation function (ACF) is the cross-correlation of a signal with
% itself. Cross correlation is the correlation of two signals at different
% time lags. A direct example of autocorrelation with matlab convention:
% lag 0: corr(x, x)
% lag 1: corr(x(1:(end-1)), x(2:end)) (note matching the number of samples)
% lag 2: corr(x(1:(end-2)), x(3:end))
% ...
% When you plot the lags on x axis and r-values of correlation on y axis,
% you get the ACF. You would want to reduce this to a single value to do
% statistics. The conventional way to do it is to take the first lag that
% gives a r value below 0.5. This is the half-life of ACF. In other words,
% how much time do we have to wait until the similarity of the signal with
% itself reduces to half of what's possibly achievable (1). Honey et al
% picked this method in their seminal paper (2012,
% https://doi.org/10.1016%2Fj.neuron.2012.08.011). For certain processes,
% this value is analytically obtainable, so it makes perfect sense to use
% ACW-50. However, in a methodological analysis, Golesorkhi et al (2021) found
% that ACW-0: the first lag where correlation drops below 0 better differentiates 
% regions in the brain (https://doi.org/10.1038/s42003-021-01785-z). Other
% methods include taking the first inflection point and fitting an exponential
% decay function and estimating the decay rate of that.
%  
% The calculation is thus very simple. The theory is as deep as it can get.
% Whole books can be written about ACF and what it means. For starters, the
% response of a system to a small enough perturbation can be estimated by
% the autocorrelation structure of the system at equilibrium. Small is
% relative. In the example of water, increasing its temperature from 10 to
% 30 can be small enough. On the other hand, increasing the temperature
% from -1 to 1 results in dramatic changes, that is too big. In that case,
% autocorrelation is not really helpful and one should turn to nonlinear
% dynamics and bifurcation theory. A long ACW will make the relaxation time
% (the time that takes for the perturbed system to go back to equilibrium)
% longer. Intuitively, if one looks at the time series, one can think that
% the system remembers the perturbation for a longer duration. Since it
% takes a long time for the perturbation the decay, a second perturbation
% can be added on top of it. Wolff et al (2022) called this temporal
% integration and the opposite (short ACW, perturbations not getting added)
% temporal segregation (https://doi.org/10.1016/j.tics.2021.11.007). Two
% sides of the same coin. 
% 
% Another interesting aspect is the Wiener-Khinchin theorem that says
% the ACF is the inverse fourier transformation of the power spectrum.
% Software such as MATLAB and Python toolboxes use this relation when
% calculating ACF since fast fourier is, as its name suggests, extremely
% fast. Kasdin (1995) uses this relationship to simulate 1/f noise, which
% ties ACW to PLE (https://doi.org/10.1109/5.381848, I should note that
% this is my favourite paper). 
% 
% With all that being said, this function doesn't do anything fancy. It just
% calculates the ACF and extracts ACW from that. As usual, it is the 
% up to the researcher to use it in any way that they find interesting. 
% 
% Input:
%   x: time series
%   fs: sampling frequency
%   isplot: plot the ACF? (logical)
% 
% Output:
%   acw_0: ACW-0
%   acw_50: ACW-50
%   acf: ACF (y axis)
%   lags: lags (x axis)
% Example usage: [acw_0, acw_50] = acw(x, 0.5, true);
% Authored by Yasir Çatal
% catalyasir@gmail.com
%% The calculation

[acf, lags] = xcorr(x, 'coeff'); 
index = find(acf == 1); % Get rid of the left-side
acf = acf(index:end); 
lags = lags(index:end);

[~, ACW_50_i] = max(acf<=0.5);
acw_50 = ACW_50_i / fs; % Divide by fs to convert samples to seconds
[~, ACW_0_i] = max(acf<=0);
acw_0 = ACW_0_i / fs;
lags = lags / fs;

if isplot
    plot(lags,acf,'k')
    xlim([0 max(lags)])
    hold on
    area(lags(1:ACW_50_i), acf(1:ACW_50_i),'FaceColor','r','FaceAlpha',0.3);
    area(lags(1:ACW_0_i), acf(1:ACW_0_i),'FaceColor','m','FaceAlpha',0.3);
    title(['ACW-0 = ', num2str(acw_0, '%.1f'), ' ACW-50 = ', num2str(acw_50, '%.1f')])
    xlabel('Lags (s)')
    ylabel('Autocorrelation')
end 
end
