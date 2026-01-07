% Copyright (c) 2021, Tharaj Thaj, Emanuele Viterbo, and Yi Hong
close all
clear all
rng(1)

%% OTFS parameters %%%%%%%%%%
N = 64;
M = 64;
% --- THAY ĐỔI 1: Bậc điều chế 16-QAM ---
M_mod = 16; 
M_bits = log2(M_mod);
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement
length_ZP = M/16;
M_data = M-length_ZP;
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
N_syms_perfram = sum(sum(data_grid));
N_bits_perfram = N_syms_perfram*M_bits;

% Time and frequency resources
car_fre = 4*10^9;
delta_f = 15*10^3; 
T = 1/delta_f; 

% --- THAY ĐỔI 2: Dải SNR cho 16-QAM ---
SNR_dB = 15:2.5:22.5; 
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;

%% Initializing simulation error count variables
N_fram = 1000;
est_info_bits_MRC=zeros(N_bits_perfram,1);
err_ber_MRC = zeros(1,length(SNR_dB));
% --- PHẦN BỔ SUNG: Khởi tạo biến đếm số lần lặp ---
no_of_detetor_iterations_MRC= zeros(length(SNR_dB),1);
avg_no_of_iterations_MRC=zeros(1,length(SNR_dB));
avg_ber_MRC = zeros(1,length(SNR_dB));

%% Normalized DFT matrix
Fn=dftmtx(N);  
Fn=Fn./norm(Fn);  
current_frame_number=zeros(1,length(SNR_dB));

for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
        current_frame_number(iesn0)=ifram;
        
        %% random input bits generation%%%%%
        trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        
        %% 2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        X = Generate_2D_data_grid(N,M,data,data_grid);
        
        %% OTFS modulation%%%%
        X_tilda=X*Fn';                     
        s = reshape(X_tilda,N*M,1);        
        
        %% OTFS channel generation (3GPP standard)
        max_speed=500;  
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        L_set=unique(delay_taps);       
        gs=Gen_discrete_time_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef);
        
        %% channel output%%%%%             
        r=zeros(N*M,1);        
        l_max=max(delay_taps);
        for q=1:N*M
            for l=(L_set+1)
                if(q>=l)
                    r(q)=r(q)+gs(l,q)*s(q-l+1);  
                end
            end
        end
        noise= sqrt(sigma_2(iesn0)/2)*(randn(size(s)) + 1i*randn(size(s)));
        r=r+noise;
                
        %% Generate channel vectors & TF channel
        [nu_ml_tilda]=Gen_delay_time_channel_vectors(N,M,l_max,gs);  
        [H_tf]=Generate_time_frequency_channel_ZP(N,M,gs,L_set);  
                
        %% MRC delay-time detection
        n_ite_MRC=50; 
        % --- THAY ĐỔI 3: Điều chỉnh Omega cho 16-QAM ---
        omega=1;
        if(M_mod>=64)
            omega=0.25;     % set omega to a smaller value (for example: 0.05) for modulation orders greater than 64-QAM
        end
        decision=1;         
        init_estimate=1;    
        
        [est_info_bits_MRC,det_iters_MRC,data_MRC] = MRC_delay_time_detector(N,M,M_data,M_mod,sigma_2(iesn0),data_grid,r,H_tf,nu_ml_tilda,L_set,omega,decision,init_estimate,n_ite_MRC);               

        %% errors count%%%%%
        errors_MRC = sum(xor(est_info_bits_MRC,trans_info_bit));                
        err_ber_MRC(iesn0) = err_ber_MRC(iesn0) + errors_MRC;        
        avg_ber_MRC(iesn0)=err_ber_MRC(iesn0)/length(trans_info_bit)/ifram;
        
        %% --- PHẦN BỔ SUNG: iterations count ---        
        no_of_detetor_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)+det_iters_MRC;
        avg_no_of_iterations_MRC(iesn0)=no_of_detetor_iterations_MRC(iesn0)/ifram;
        
        %% --- PHẦN BỔ SUNG: DISP error performance details ---        
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

figure(1)
semilogy(SNR_dB,avg_ber_MRC,'-or','LineWidth',2,'MarkerSize',8)
legend('MRC - 16QAM')
grid on; xlabel('SNR(dB)'); ylabel('BER');