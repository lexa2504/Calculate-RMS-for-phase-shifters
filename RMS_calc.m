clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Const
number_bits = 6; 
min_freq    = 6e9;
max_freq    = 16e9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read S-parameters

numfiles    = pow2(number_bits);
phase_step  = 360/numfiles;
filename    = "PS_test__"+(1:numfiles)+".s2p"; % Construct filenames
S           = sparameters(filename(1)); % Read file #1 for initial set-up
freq        = S.Frequencies; % Frequency values are the same for all files
numfreq     = numel(freq); % Number of frequency points
s21_data    = zeros(numfreq,numfiles); % Preallocate for speed
idx         = (freq >= min_freq) & (freq <= max_freq);

%% Read Touchstone files
for n   = 1:numfiles
    S   = sparameters(filename(n));
    s21 = rfparam(S,2,1);
    s21_data(:,n) = s21;
    s11 = rfparam(S,1,1);
    s11_data(:,n) = s11;
    s22 = rfparam(S,2,2);
    s22_data(:,n) = s22;
    s12 = rfparam(S,1,2);
    s12_data(:,n) = s12;
end
%% Work with data S21 db
s21_db = 20*log10(abs(s21_data));
s21_degree =(angle(s21_data));
s11_db = 20*log10(abs(s11_data));
s22_db = 20*log10(abs(s22_data));
s12_db = 20*log10(abs(s12_data));
%
%% Work with data S21 degree
s21_pass_data   = s21_data(idx,:);
s21_pass_degree = unwrap(s21_degree(idx,:))*180/pi;
freq_pass_ghz   = freq(idx)/1e9;    % Normalize to GHz  
%%
% %Graph S21_db, s11_db, s22_db, S12_db
% 
t = tiledlayout(2,2);
nexttile
plot(freq/1e9,s21_db)
xlabel('Frequency (GHz)')
ylabel('S_2_1 (dB)')
grid on
nexttile
plot(freq/1e9,s11_db)
xlabel('Frequency (GHz)')
ylabel('S_1_1 (dB)')
grid on
nexttile
plot(freq/1e9,s22_db)
xlabel('Frequency (GHz)')
ylabel('S_2_2 (dB)')
grid on
nexttile
plot(freq/1e9,s12_db)
xlabel('Frequency (GHz)')
ylabel('S_1_2 (dB)')
grid on


%%
%Graph S21_degree
figure
plot(freq_pass_ghz,s21_pass_degree)
xlabel('Frequency (GHz)')
ylabel('S_2_1 (degree)')
%title('S21 в градусах для 64 состояний')
axis on
grid on
%%
%%Normirovka and sortirovka
[value,position] = max(s21_pass_degree(1,:));
for i = 1:numfiles
    if i == position
        S_ph(:,i)=0;
        else 
            S_ph(:,i) = s21_pass_degree(:,position) - s21_pass_degree(:,i);
    end
end
numb = length(freq_pass_ghz);
for i = 1:numb
s21_data_norm(i,:) = sort(S_ph(i,:));
end
%%
% Plot norm phase
figure
plot(freq_pass_ghz,(s21_data_norm))
xlabel('Frequency (GHz)')
ylabel('S21 (degree)')
%title('S21 в градусах для 64 состояний')
grid on
%%
%ideal_phase_shifts
for j = 1:numfiles
    for k = 1:numb   %найти потом как 200 чтоб автоматом брало
    phase_step_ideal(k,j) = (j-1)*phase_step;
    end
end
%%
% RMS phase CALC
phase_error=s21_data_norm - phase_step_ideal;

for i=1:numb
    m_phase_error(i,:) =mean(phase_error(i,:));
end

for i=1:(numfiles)
    m_phase_error_all(:,i) = m_phase_error(:,1);
end

phase_error_abs = phase_error - m_phase_error;

for i=1:numfiles
      phase_error_abs_mul(:,i)= phase_error_abs(:,i).*phase_error_abs(:,i);
end

for i=1:k
      phase_error_abs_summ(i,:)= sum(phase_error_abs_mul(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMS_deviation = sqrt(phase_error_abs_summ/(numfiles-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Plot RMS phase
freq_pass_ghz   = freq(idx)/1e9;    % Normalize to GHz    
figure
plot(freq_pass_ghz,(RMS_deviation))
xlabel('Frequency (GHz)')
ylabel('RMS (degree)')
title('RMS deviation phase')
grid on
%%
%RMS magnitude CALC

k = length(freq);
for i=1:k
    mean_ampl_error(i,:) =mean(s21_db(i,:));
end

for i=1:(numfiles)
    mean_ampl_error_all(:,i) = mean_ampl_error(:,1);
end

Mag_error_db = s21_db - mean_ampl_error_all;

for i=1:numfiles
     Mag_error_db_mul(:,i)= Mag_error_db(:,i).*Mag_error_db(:,i);
end

for i=1:k
      Mag_error_db_summ(i,:)= sum(Mag_error_db_mul(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMS_mag_deviation = sqrt(Mag_error_db_summ/(numfiles-1));
%%%%%%%%%%%%
RMS_mag_deviation   = RMS_mag_deviation(idx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Plor RMS magnitude
freq_pass_ghz   = freq(idx)/1e9;    % Normalize to GHz    
figure
plot(freq_pass_ghz,RMS_mag_deviation)
xlabel('Frequency (GHz)')
ylabel('RMS (dB)')
title('RMS deviation magnitude')
grid on
%%