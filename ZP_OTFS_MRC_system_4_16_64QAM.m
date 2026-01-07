% Copyright (c) 2021, Tharaj Thaj, Emanuele Viterbo, and Yi Hong
close all; clear all; rng(1);

%% 1. Tham số chung
N = 64; M = 64;
M_mod_list = [4, 16, 64]; 
SNR_dB_cell = {10:2.5:20, 15:2.5:22.5, 25:2.5:35}; % Dải SNR riêng biệt cho từng QAM
N_fram = 1000; 

% Lưu trữ để vẽ đồ thị
results_BER = cell(1, length(M_mod_list));
plot_styles = {'-xb', '-or', '-dg'}; 

for idx_m = 1:length(M_mod_list)
    M_mod = M_mod_list(idx_m);
    M_bits = log2(M_mod);
    SNR_dB = SNR_dB_cell{idx_m};
    
    % Năng lượng và lưới dữ liệu
    eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
    length_ZP = M/16;
    M_data = M-length_ZP;
    data_grid = zeros(M,N);
    data_grid(1:M_data,1:N) = 1;
    N_syms_perfram = sum(sum(data_grid));
    N_bits_perfram = N_syms_perfram*M_bits;
    
    % Khởi tạo biến thống kê
    err_ber_MRC = zeros(1,length(SNR_dB));
    avg_ber_MRC = zeros(1,length(SNR_dB));
    no_of_detetor_iterations_MRC = zeros(length(SNR_dB),1);
    avg_no_of_iterations_MRC = zeros(1,length(SNR_dB));
    current_frame_number = zeros(1,length(SNR_dB));
    
    sigma_2 = (abs(eng_sqrt)^2)./(10.^(SNR_dB/10));
    Fn = dftmtx(N); Fn = Fn./norm(Fn);

    for iesn0 = 1:length(SNR_dB)
        for ifram = 1:N_fram
            current_frame_number(iesn0) = ifram;
            
            %% Transmitter
            trans_info_bit = randi([0,1], N_bits_perfram, 1);
            data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), M_mod, 'gray', 'InputType', 'bit');
            X = Generate_2D_data_grid(N, M, data, data_grid);
            s = reshape(X*Fn', N*M, 1);
            
            %% Channel (3GPP 500km/h)
            max_speed = 500; car_fre = 4e9; delta_f = 15e3; T = 1/delta_f;
            [chan_coef, delay_taps, Doppler_taps, taps] = Generate_delay_Doppler_channel_parameters(N, M, car_fre, delta_f, T, max_speed);
            L_set = unique(delay_taps);       
            gs = Gen_discrete_time_channel(N, M, taps, delay_taps, Doppler_taps, chan_coef);
            
            r = zeros(N*M, 1); l_max = max(delay_taps);
            for q = 1:N*M
                for l = (L_set+1)
                    if(q >= l), r(q) = r(q) + gs(l,q)*s(q-l+1); end
                end
            end
            r = r + sqrt(sigma_2(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
            
            %% Detector Parameters
            n_ite_MRC = 50; 
            % Tối ưu omega và khởi tạo cho từng bậc điều chế
            if M_mod == 4, omega = 1; init_estimate = 1;
            elseif M_mod == 16, omega = 0.8; init_estimate = 1;
            else, omega = 0.25; init_estimate = 0; end % 64-QAM
            
            [nu_ml_tilda] = Gen_delay_time_channel_vectors(N, M, l_max, gs);
            [H_tf] = Generate_time_frequency_channel_ZP(N, M, gs, L_set);
            
            [est_bits, det_iters, ~] = MRC_delay_time_detector(N,M,M_data,M_mod,sigma_2(iesn0),data_grid,r,H_tf,nu_ml_tilda,L_set,omega,1,init_estimate,n_ite_MRC);               

            %% Thống kê lỗi và Interactions (Sửa lỗi tên biến trans_info_bit)
            errors_frame = sum(xor(est_bits, trans_info_bit));                
            err_ber_MRC(iesn0) = err_ber_MRC(iesn0) + errors_frame;        
            avg_ber_MRC(iesn0) = err_ber_MRC(iesn0) / (N_bits_perfram * ifram);
            
            no_of_detetor_iterations_MRC(iesn0) = no_of_detetor_iterations_MRC(iesn0) + det_iters;
            avg_no_of_iterations_MRC(iesn0) = no_of_detetor_iterations_MRC(iesn0) / ifram;
            
            %% Hiển thị (DISP)
            if mod(ifram, 50) == 0 || ifram == N_fram
                clc; disp('####################################################################') 
                fprintf('ZP-OTFS Simulation - Current Modulation: %d-QAM\n', M_mod);
                display(SNR_dB, 'SNR (dB)');
                display(current_frame_number, 'Number of frames');
                display(avg_ber_MRC, 'Average BER - MRC Detector');        
                display(avg_no_of_iterations_MRC, 'Average iterations');       
                disp('####################################################################')
            end
        end
    end
    results_BER{idx_m} = avg_ber_MRC;
end

%% 2. Vẽ đồ thị tổng hợp
figure(1)
for idx_m = 1:length(M_mod_list)
    semilogy(SNR_dB_cell{idx_m}, results_BER{idx_m}, plot_styles{idx_m}, 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
end
grid on; xlabel('SNR(dB)'); ylabel('BER'); axis([10 35 1e-6 1e-1]);
legend('MRC-4QAM', 'MRC-16QAM', 'MRC-64QAM'); 
title('Performance Comparison of ZP-OTFS (MRC) for 4, 16, 64-QAM');