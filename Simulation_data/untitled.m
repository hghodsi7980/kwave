input = zeros(31,100);
output = zeros(31,2);
k = 1;
for i = 0:30
    i
    file_name = sprintf("Mixed_1_%d.mat",i);
    load(file_name);
    hold on
    frequencies = frequencies(1:121);
    fft = fft_r(1000,1:121);
    plot(frequencies,fft/sum(fft));
    [peaks, indexes] = findpeaks(fft);
    freqs = frequencies(indexes);
    r1 (k,i+1) = freqs(1);
    r2(k,i+1) = freqs(2);
    % r1(k,i+1) = fft(1,6)/sum(fft);
    % r2(k,i+1) = fft(1,21)/sum(fft);
    % r3(k,i+1) = sum(fft);
    % por(k,i+1) = sum(sum(output_mask_clot))/(pi*250^2);
    Fibr(k,i+1) = 90-i*3;
    % por(k,i+1) = por(k,i+1)/2; 
    ax = gca;
    % Set font size for axes labels and title
    ax.FontSize = 24;  % Change the font size according to your preference
    ax.Title.FontSize = 24;  % Font size for the title
    xlabel Frequency(Hz)
    ylabel Amplitude(a.u.)
end
% k = k+1;
% for i = 0:30
%     i
%     file_name = sprintf("Fibrin_2_%d.mat",i);
%     load(file_name);
%     hold on
%     frequencies = frequencies(1:121);
%     fft = fft_r(1000,1:121);
%     plot(frequencies,fft/sum(fft));
%     r1(k,i+1) = fft(1,6)/sum(fft);
%     r2(k,i+1) = fft(1,21)/sum(fft);
%     r3(k,i+1) = sum(fft);
%     por(k,i+1) = sum(sum(output_mask_clot))/(pi*250^2);
%     por(k,i+1) = por(k,i+1)/2; 
%     ax = gca;
%     % Set font size for axes labels and title
%     ax.FontSize = 24;  % Change the font size according to your preference
%     ax.Title.FontSize = 24;  % Font size for the title
%     xlabel Frequency(Hz)
%     ylabel Amplitude(a.u.)
% end
% k = k+1;
% for i = 0:30
%     i
%     file_name = sprintf("Fibrin_3_%d.mat",i);
%     load(file_name);
%     hold on
%     frequencies = frequencies(1:121);
%     fft = fft_r(1000,1:121);
%     plot(frequencies,fft/sum(fft));
%     r1(k,i+1) = fft(1,6)/sum(fft);
%     r2(k,i+1) = fft(1,21)/sum(fft);
%     r3(k,i+1) = sum(fft);
%     por(k,i+1) = sum(sum(output_mask_clot))/(pi*250^2);
%     por(k,i+1) = por(k,i+1)/2; 
%     ax = gca;
%     % Set font size for axes labels and title
%     ax.FontSize = 24;  % Change the font size according to your preference
%     ax.Title.FontSize = 24;  % Font size for the title
%     xlabel Frequency(Hz)
%     ylabel Amplitude(a.u.)
% end
% k = k+1;
% for i = 0:30
%     i
%     file_name = sprintf("Fibrin_4_%d.mat",i);
%     load(file_name);
%     hold on
%     frequencies = frequencies(1:121);
%     fft = fft_r(1000,1:121);
%     plot(frequencies,fft/sum(fft));
%     r1(k,i+1) = fft(1,6)/sum(fft);
%     r2(k,i+1) = fft(1,21)/sum(fft);
%     r3(k,i+1) = sum(fft);
%     por(k,i+1) = sum(sum(output_mask_clot))/(pi*250^2);
%     por(k,i+1) = por(k,i+1)/2; 
%     ax = gca;
%     % Set font size for axes labels and title
%     ax.FontSize = 24;  % Change the font size according to your preference
%     ax.Title.FontSize = 24;  % Font size for the title
%     xlabel Frequency(Hz)
%     ylabel Amplitude(a.u.)
% end
% k = k+1;
% for i = 0:30
%     i
%     file_name = sprintf("Fibrin_5_%d.mat",i);
%     load(file_name);
%     hold on
%     frequencies = frequencies(1:121);
%     fft = fft_r(1000,1:121);
%     plot(frequencies,fft/sum(fft));
%     r1(k,i+1) = fft(1,6)/sum(fft);
%     r2(k,i+1) = fft(1,21)/sum(fft);
%     r3(k,i+1) = sum(fft);
%     por(k,i+1) = sum(sum(output_mask_clot))/(pi*250^2);
%     por(k,i+1) = por(k,i+1)/2; 
%     ax = gca;
%     % Set font size for axes labels and title
%     ax.FontSize = 24;  % Change the font size according to your preference
%     ax.Title.FontSize = 24;  % Font size for the title
%     xlabel Frequency(Hz)
%     ylabel Amplitude(a.u.)
% end

