% Dựa trên nghiên cứu của Tharaj Thaj, Emanuele Viterbo, và Yi Hong
close all; clear all; rng(1);

%% 1. Tham số OTFS
N = 64; M = 64; 
M_mod = 64; % Bậc điều chế 64-QAM
M_bits = log2(M_mod);
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% 2. Cấu trúc lưới Delay-Doppler và Zero-Padding
length_ZP = M/16; 
M_data = M-length_ZP; 
data_grid = zeros(M,N);
data_grid(1:M_data,1:N) = 1;
N_syms_perfram = sum(sum(data_grid));
N_bits_perfram = N_syms_perfram*M_bits;

% Thông số tài nguyên
car_fre = 4*10^9; delta_f = 15*10^3; T = 1/delta_f; 

%% 3. Cấu hình SNR cho 64-QAM
SNR_dB = 25:2.5:35; % Dải SNR theo hình mẫu bạn mong muốn
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;

%% 4. Khởi tạo biến mô phỏng
N_fram = 1000; 
err_ber_MRC = zeros(1,length(SNR_dB));
avg_ber_MRC = zeros(1,length(SNR_dB));
no_of_detetor_iterations_MRC = zeros(length(SNR_dB),1);
avg_no_of_iterations_MRC = zeros(1,length(SNR_dB));
Fn = dftmtx(N); Fn = Fn./norm(Fn);
current_frame_number = zeros(1,length(SNR_dB));

for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
        current_frame_number(iesn0) = ifram;
        
        % Transmitter
        trans_info_bit = randi([0,1], N_bits_perfram, 1);
        data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), M_mod, 'gray', 'InputType', 'bit');
        X = Generate_2D_data_grid(N, M, data, data_grid);
        X_tilda = X*Fn'; s = reshape(X_tilda, N*M, 1);
        
        % Kênh truyền 3GPP (500 km/h)
        max_speed = 500;
        [chan_coef, delay_taps, Doppler_taps, taps] = Generate_delay_Doppler_channel_parameters(N, M, car_fre, delta_f, T, max_speed);
        L_set = unique(delay_taps);       
        gs = Gen_discrete_time_channel(N, M, taps, delay_taps, Doppler_taps, chan_coef);
        
        % Tín hiệu nhận
        r = zeros(N*M, 1); l_max = max(delay_taps);
        for q = 1:N*M
            for l = (L_set+1)
                if(q >= l), r(q) = r(q) + gs(l,q) * s(q-l+1); end
            end
        end
        r = r + sqrt(sigma_2(iesn0)/2) * (randn(size(s)) + 1i*randn(size(s)));
                
        %% 5. Bộ thu MRC (Phần sửa lỗi biến n_ite_MRC)
        [nu_ml_tilda] = Gen_delay_time_channel_vectors(N, M, l_max, gs);
        [H_tf] = Generate_time_frequency_channel_ZP(N, M, gs, L_set);
        
        % ĐỊNH NGHĨA BIẾN TẠI ĐÂY
        n_ite_MRC = 50;     % Số lần lặp tối đa
        omega=1;
        if(M_mod==64)
            omega=0.25;     % set omega to a smaller value (for example: 0.05) for modulation orders greater than 64-QAM
        end
        decision = 1;       % Hard decision
        init_estimate = 0;  % Bắt đầu từ 0 cho 64-QAM giúp ổn định hơn
        
        % Gọi hàm detector
        [est_info_bits_MRC, det_iters_MRC, ~] = MRC_delay_time_detector(N,M,M_data,M_mod,sigma_2(iesn0),data_grid,r,H_tf,nu_ml_tilda,L_set,omega,decision,init_estimate,n_ite_MRC);               

        %% 6. Thống kê lỗi và Hiển thị (DISP)
        errors_MRC = sum(xor(est_info_bits_MRC,trans_info_bit));                
        err_ber_MRC(iesn0) = err_ber_MRC(iesn0) + errors_MRC;        
        avg_ber_MRC(iesn0)=err_ber_MRC(iesn0)/length(trans_info_bit)/ifram;

        no_of_detetor_iterations_MRC(iesn0) = no_of_detetor_iterations_MRC(iesn0) + det_iters_MRC;
        avg_no_of_iterations_MRC(iesn0) = no_of_detetor_iterations_MRC(iesn0) / ifram;
        
        clc
        disp('####################################################################') 
        fprintf('ZP-OTFS-(N,M,QAM size)');disp([N,M,M_mod]);
        display(SNR_dB,'SNR (dB)');
        display(current_frame_number,'Number of frames');
        display(avg_ber_MRC,'Average BER - Delay-time domain MRC');        
        display(avg_no_of_iterations_MRC,'Average number of iterations for the MRC detector');       
        disp('####################################################################')
    end
end

% Vẽ đồ thị
figure(1); semilogy(SNR_dB, avg_ber_MRC, '-dg', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'y');
grid on; xlabel('SNR(dB)'); ylabel('BER'); title('Performance of ZP-OTFS with 64-QAM (MRC)');
legend('MRC-64QAM'); axis([25 35 1e-6 1e-1]);