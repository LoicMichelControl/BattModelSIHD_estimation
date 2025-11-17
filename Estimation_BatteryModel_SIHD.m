    
    %=================================================================
    %=================================================================
    
    % Semi-Implicit Homogeneous Differentiation (SIHD) estimation
    % applied to a simple battery model – Preliminary results
    %
    % This code illustrates the first application of the SIHD technique
    % to the estimation of a Li-ion battery equivalent circuit model
    % composed of four independent first-order ODEs.
    %
    % The SIHD estimator itself does NOT embed the battery equations. 
    % To obtain a prediction of the states, SIHD is supplied with a 
    % model-based duplicate (a "minimal" Luenberger observer), 
    % which acts as a “model surrogate” for the prediction step.
    % In this version, the terminal voltage Vt is not yet injected into
    % the SIHD estimation mechanism: only the duplicated model is used
    % to generate the prediction errors on the independent ODE branches.
    %
    % The purpose is to isolate the intrinsic behavior of SIHD  before 
    % introducing Vt-fed synchronization in future work.
    %
    % (c) 2025 — Nantes Université, Centrale Nantes, LS2N UMR 6004, France
    %
    % This is ongoing research and the code is provided for reproducible
    % preliminary proof-of-concept results.   
    %=================================================================
    %=================================================================
    
    clear all; clc; close all;
    
    %% ----------------- Simulation settings  -------
    TMax = 2000 ; % Maximum simulation time
    Deltah = 1e-3;  % Time-step
    
    NoiseAmp = 0.2;    % measurement noise amplitude on Vt [V]
    
    %% ----------------- Battery parameters (example values) -------------------
    
    % Model of the battery
    % dSoC = - I / Q_nom;
    % dV_RC1 = - V_RC1/(R1*C1) + I/C1;
    % dV_RC2 = - V_RC2/(R2*C2) + I/C2;
    % dT = (I^2*(R_s+R1+R2) - (T - T_amb)/R_th)/C_th;
    
    Q_nom = 5*3600;
    R_s = 0.01;
    R1 = 0.015; C1 = 2400;
    R2 = 0.02;  C2 = 10000;
    C_th = 500; R_th = 1.5; T_amb = 298;
    
    % OCV (nonlinear via pchip)
    SOC_data = [0 0.1 0.2 0.4 0.6 0.8 1.0];
    OCV_data = [3.0 3.2 3.4 3.7 3.9 4.1 4.2];
    Voc = @(SoC) interp1(SOC_data, OCV_data, SoC, 'pchip');
    
    % Input current (DC + multisine) - sign convention: I>0 discharge, I<0 charge
    I_DC = -1; f1 = 1/300; f2 = 1/60; f3 = 1/10;
    A1 = 0.3; A2 = 0.1; A3 = 0.05;
    
    I_input = @(t) -1*(1 - exp( -(t - Deltah)/0.1)) + A1*sin(2*pi*f1*(t - Deltah)) + A2*sin(2*pi*f2*(t - Deltah)) + A3*sin(2*pi*f3*(t - Deltah)) + 1*sin(2*pi*0.002*(t - Deltah));
    % Note: I_input(t) is evaluated inside the main loop.

    
    %% ----------------- Initialization of the battery model (x_model) ----------------
    
    SoC0 = 0;   % SoC
    V10 = -0.3; % V_RC1
    V20 = 0.1;  % V_RC2
    T0 = -0.4;  % Temperature
    
    %% ----------------- SIHD tuning parameters ----------------
    
    alpha_1 = 0.5;
    alpha_2 = 0.9;
    
    lambda_1 = 10;
    lambda_2 = 1e4;
    
    % Diagnostic vectors
    time = 0;
    uu = 0;
    rr = 0;
    time_vec = [];
    time_vec_sampled = [];
    
    % noise vector (pre-generate)
    eta_noise = NoiseAmp * randn(1, ceil(TMax/Deltah)+10);
    
    % logs
    maxIter = ceil(TMax/Deltah)+5;
    SoC_true_vec = zeros(1,maxIter);
    VRC1_true_vec = zeros(1,maxIter);
    VRC2_true_vec = zeros(1,maxIter);
    T_true_vec = zeros(1,maxIter);
    Vt_meas_vec = zeros(1,maxIter);
    I_vec = zeros(1,maxIter);
    
    % sampled logs
    cz1_d_vec = [];
    cz2_d_vec = [];
    bcz1_d_vec = [];
    bcz2_d_vec = [];
    ce_1_vec = [];
    ce_1_ymeas_vec = [];
    x_model = [SoC0; V10; V20; T0];
    x_model_ = [0; 0; 0; 0];
    
    %% ----------------- Observer initialization ------------------------
    
    % SIHD initial conditions (every parameter set to zero)
    cz1_d_a = 0;
    cz1_d_b = 0;
    cz1_d_c = 0;
    cz1_d_d = 0;
    
    cz2_d_a = 0;
    cz2_d_b = 0;
    cz2_d_c = 0;
    cz2_d_d = 0.2;
    
    cz2_p_a = 0;
    cz2_p_b = 0;
    cz2_p_c = 0;
    cz2_p_d = 0;
    
    % Internal projector flag init.
    save_cE_1_a = 0;
    save_cE_1_b = 0;
    save_cE_1_c = 0;
    save_cE_1_d = 0;
    
    % SIHD Parameters definition
    c_lambda_1_1_a = lambda_1;
    c_lambda_1_1_b = lambda_1;
    c_lambda_1_1_c = lambda_1;
    c_lambda_1_1_d = lambda_1;
    
    c_lambda_2_1_a = lambda_2;
    c_lambda_2_1_b = lambda_2;
    c_lambda_2_1_c = lambda_2;
    c_lambda_2_1_d = lambda_2;
    
    c_alpha_1_1_a = alpha_1;
    c_alpha_1_1_b = alpha_1;
    c_alpha_1_1_c = alpha_1;
    c_alpha_1_1_d = alpha_1;
    
    c_alpha_2_1_a = alpha_2;
    c_alpha_2_1_b = alpha_2;
    c_alpha_2_1_c = alpha_2;
    c_alpha_2_1_d = alpha_2;
    
    %% ----------------- Main loop -----------------
    while time <= TMax
        uu = uu + 1;
        time = Deltah * uu;
        time_vec(uu) = time;
    
        % % --- Input current definition (starts with a ramp and then a constant) ---
        if ( time < 1500)
            I = I_input(time) + eta_noise(uu);
        else
            I = -1e-3 + eta_noise(uu);
        end
    
        I_vec(uu) = I;

        % --- Battery model system update (Euler explicit) ---
        dx1  = 1e-4*x_model(1) - I / Q_nom;  % small leakage term to avoid zero dynamics
        dx2  = - x_model(2)/(R1*C1) + I/C1;
        dx3  = - x_model(3)/(R2*C2) + I/C2;
        dx4  = - ((x_model(4) - T_amb)/R_th)/C_th   + I^2*(R_s+R1+R2)/C_th ;
        x_model = x_model + Deltah*[dx1; dx2; dx3; dx4];
    
        % --- Duplication of the original Battery model (considering zero initial
        % conditions)
        dx1_ = 1e-4*x_model_(1) - I / Q_nom; % small leakage term to avoid zero dynamics
        dx2_ = - x_model_(2)/(R1*C1) + I/C1;
        dx3_ = - x_model_(3)/(R2*C2) + I/C2;
        dx4_ = - ((x_model_(4) - T_amb)/R_th)/C_th   + I^2*(R_s+R1+R2)/C_th ;
        x_model_ = x_model_ + Deltah*[dx1_; dx2_; dx3_; dx4_];
    
    
        xlog(:,uu) = x_model;
    
        % --- Calculate time-derivatives (differences) of each state
        if ( uu > 1 )
            Delta_(1,uu) = (xlog(1,uu) - xlog(1,uu-1))/ Deltah;
            Delta_(2,uu) = (xlog(2,uu) - xlog(2,uu-1))/ Deltah;
            Delta_(3,uu) = (xlog(3,uu) - xlog(3,uu-1))/ Deltah;
            Delta_(4,uu) = (xlog(4,uu) - xlog(4,uu-1))/ Deltah;
        end
    
        % --- Sampled observer prediction (every h_mult * Deltah) ---
    
        % Note: The prediction of each state knowing only the resulting Vt is a work in progress ;)
    
            rr = rr + 1;

            % Computation of each state error (difference between the duplicated model
            ce_1_a = x_model_(1,end) - cz1_d_a; 
            ce_1_b = x_model_(2,end) - cz1_d_b;
            ce_1_c = x_model_(3,end) - cz1_d_c;
            ce_1_d = x_model_(4,end) - cz1_d_d;
    
    
            c_lambda_p1_1_a = c_lambda_1_1_a;
            c_lambda_p1_1_b = c_lambda_1_1_b;
            c_lambda_p1_1_c = c_lambda_1_1_c;
            c_lambda_p1_1_d = c_lambda_1_1_d;
    
            c_lambda_p2_1_a = c_lambda_2_1_a;
            c_lambda_p2_1_b = c_lambda_2_1_b;
            c_lambda_p2_1_c = c_lambda_2_1_c;
            c_lambda_p2_1_d = c_lambda_2_1_d;

            % --- Projectors management: Force projectors flags to '1' if time >= 100
            if ( time >= 100)
        
                cE_1_a = 1;
                cE_1_b = 1;
                cE_1_c = 1;
                cE_1_d = 1;
    
            else
    
                cE_1_a = 0;
                cE_1_b = 0;
                cE_1_c = 0;
                cE_1_d = 0;

            end

            % and maintain the flag to '1'
            if ( cE_1_a == 1 || save_cE_1_a == 1)
                cE_1_a = 1;
                save_cE_1_a = 1;
            end
    
            if ( cE_1_b == 1 || save_cE_1_b == 1)
                cE_1_b = 1;
                save_cE_1_b = 1;
            end
    
            if ( cE_1_c == 1 || save_cE_1_c == 1)
                cE_1_c = 1;
                save_cE_1_c = 1;
            end
    
            if ( cE_1_d == 1 || save_cE_1_d == 1)
                cE_1_d = 1;
                save_cE_1_d = 1;
            end
    
            % --- First observer stage: correct SoC, V_RC1, V_RC2 and T  ---
    
            [ cE_1_a, cProj_1_a , borne_1_a ] = Proj_function( c_alpha_1_1_a, c_lambda_p1_1_a, 1, ce_1_a, 1, Deltah);
    
            [ cE_1_b, cProj_1_b , borne_1_b ] = Proj_function( c_alpha_1_1_b, c_lambda_p1_1_b, 1, ce_1_b, 1, Deltah);
    
            [ cE_1_c, cProj_1_c , borne_1_c ] = Proj_function( c_alpha_1_1_c, c_lambda_p1_1_c, 1, ce_1_c, 1, Deltah);
    
            [ cE_1_d, cProj_1_d , borne_1_d ] = Proj_function( c_alpha_1_1_d, c_lambda_p1_1_d, 1, ce_1_d, 1, Deltah);
    
    
            cz1_p_a = cz1_d_a + Deltah * ( cz2_p_a + c_lambda_1_1_a * ( abs( ce_1_a ) )^c_alpha_1_1_a * cProj_1_a );
    
            cz1_p_b = cz1_d_b + Deltah * ( cz2_p_b + c_lambda_1_1_b * ( abs( ce_1_b ) )^c_alpha_1_1_b * cProj_1_b );
    
            cz1_p_c = cz1_d_c + Deltah * ( cz2_p_c + c_lambda_1_1_c * ( abs( ce_1_c ) )^c_alpha_1_1_c * cProj_1_c );
    
            cz1_p_d = cz1_d_d + Deltah * ( cz2_p_d + c_lambda_1_1_d * ( abs( ce_1_d ) )^c_alpha_1_1_d * cProj_1_d );
    
            % --- Second observer stage: correct SoC, V_RC1, V_RC2 and T  ---
    
            [ cE_2_a, cProj_2_a , borne_2_a ] = Proj_function( c_alpha_2_1_a, c_lambda_p2_1_a, 1, ce_1_a, 1, Deltah);
    
            [ cE_2_b, cProj_2_b , borne_2_b ] = Proj_function( c_alpha_2_1_b, c_lambda_p2_1_b, 1, ce_1_b, 1, Deltah);
    
            [ cE_2_c, cProj_2_c , borne_2_c ] = Proj_function( c_alpha_2_1_c, c_lambda_p2_1_c, 1, ce_1_c, 1, Deltah);
            
            [ cE_2_d, cProj_2_d , borne_2_d ] = Proj_function( c_alpha_2_1_d, c_lambda_p2_1_d, 1, ce_1_d, 1, Deltah);
    

            cz2_p_a = cz2_d_a + ( cE_1_a ) * Deltah * ( c_lambda_2_1_a * ( abs( ce_1_a ) )^(2 * c_alpha_2_1_a - 1 ) * cProj_2_a );
    
            cz2_p_b = cz2_d_b + ( cE_1_b ) * Deltah * ( c_lambda_2_1_b * ( abs( ce_1_b ) )^(2 * c_alpha_2_1_b - 1 ) * cProj_2_b );

            cz2_p_c = cz2_d_c + ( cE_1_c ) * Deltah * ( c_lambda_2_1_c * ( abs( ce_1_c ) )^(2 * c_alpha_2_1_c - 1 ) * cProj_2_c );
    
            cz2_p_d = cz2_d_d + ( cE_1_d ) * Deltah * ( c_lambda_2_1_d * ( abs( ce_1_d ) )^(2 * c_alpha_2_1_d - 1 ) * cProj_2_d );
    
    
            % not used now for the estimation (see the note l.198):
            % Vt = Voc(x(1)) - R_s*I - x(2) - x(3);
            % Vt = Voc(cz2_p_a) - R_s*I - cz2_d_b - cz2_p_c;
    
            % --- Update observer states  ---
            cz2_d_a = cz2_p_a;
            cz2_d_b = cz2_p_b;
            cz2_d_c = cz2_p_c;
            cz2_d_d = cz2_p_d;
    
            cz1_d_a = cz1_p_a;
            cz1_d_b = cz1_p_b;
            cz1_d_c = cz1_p_c;
            cz1_d_d = cz1_p_d;
 
            % store sampled estimates
            cz1_d_a_vec(rr) = cz1_d_a;
            cz1_d_b_vec(rr) = cz1_d_b;
            cz1_d_c_vec(rr) = cz1_d_c;
            cz1_d_d_vec(rr) = cz1_d_d;
    
            cz2_d_a_vec(rr) = cz2_d_a;
            cz2_d_b_vec(rr) = cz2_d_b;
            cz2_d_c_vec(rr) = cz2_d_c;
            cz2_d_d_vec(rr) = cz2_d_d;
    
            cE_1_a_vec(uu) = cE_1_a;
            cE_1_b_vec(uu) = cE_1_b;
            cE_1_c_vec(uu) = cE_1_c;
            cE_1_d_vec(uu) = cE_1_d;
    
       
    end
    
    % ----------------- Post-processing / plotting ------------------------------
    
    FtSize = 30;
    
    
    figure('Name','States vs Estimates','NumberTitle','off','Position',[100 100 900 600]);
    subplot(4,1,1);
    plot(time_vec, 100*xlog(1,:),'b','LineWidth',2); hold on;
    plot(time_vec, 100*cz1_d_a_vec(1,:),'--r','LineWidth',2);
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    ylabel('$\mathrm{SoC}$ [$\%$]','fontsize', FtSize, 'Interpreter','latex');
    grid on;
    xlim([0, 2000])
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    subplot(4,1,2);
    plot(time_vec, xlog(2,:),'b','LineWidth',2); hold on
    plot(time_vec, cz1_d_b_vec(1,:),'--r','LineWidth',2);
    xlim([0, 2000])
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    ylabel('$V_{RC1}$ [V]','fontsize', FtSize, 'Interpreter','latex');
    grid on;
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    subplot(4,1,3);
    plot(time_vec, xlog(3,:),'b','LineWidth',2); hold on
    plot(time_vec, cz1_d_c_vec(1,:),'--r','LineWidth',2);
    xlim([0, 2000])
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    ylabel('$V_{RC2}$ [V]','fontsize', FtSize, 'Interpreter','latex');
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    grid on;
    subplot(4,1,4);
    plot(time_vec, xlog(4,:),'b','LineWidth',1.5); hold on;
    plot(time_vec, cz1_d_d_vec(1,:),'--r','LineWidth',1.5);
    xlim([0, 2000])
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    xlabel('Time [s]','FontSize', FtSize, 'Interpreter','latex');
    ylabel('$T$ [K]','fontsize', FtSize, 'Interpreter','latex');
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    
    
    figure('Name','DStates vs Estimates','NumberTitle','off','Position',[500 500 900 600]);
    subplot(4,1,1);
    plot(time_vec, Delta_(1,:),'b','LineWidth',2);hold on
    plot(time_vec, cz2_d_a_vec(1,:),'--r','LineWidth',2); hold on
    xlim([0, 2000])
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    ylabel('$\dot{\mathrm{SoC}}$ [$\%$/s]','fontsize', FtSize, 'Interpreter','latex');
    grid on;
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    subplot(4,1,2);
    plot(time_vec, Delta_(2,:),'b','LineWidth',2); hold on
    plot(time_vec, cz2_d_b_vec(1,:),'--r','LineWidth',2); hold on;
    xlim([0, 2000])
    ylim([-1e-3, 1e-3])
    ylabel('$\dot{V}_{RC1}$ [V/s]','fontsize', FtSize, 'Interpreter','latex');
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    grid on;
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    subplot(4,1,3);
    plot(time_vec, Delta_(3,:),'b','LineWidth',2); hold on
    plot(time_vec, cz2_d_c_vec(1,:),'--r','LineWidth',2); hold on;
    xlim([0, 2000])
    ylabel('$\dot{V}_{RC2}$ [V/s]','fontsize', FtSize, 'Interpreter','latex');
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    grid on;
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    subplot(4,1,4);
    plot(time_vec, Delta_(4,:),'b','LineWidth',2); hold on;
    plot(time_vec, cz2_d_d_vec(1,:),'--r','LineWidth',2); hold on
    xlim([0, 2000])
    legend('model','estimated (SIHD method)','fontsize', FtSize, 'Interpreter','latex')
    xlabel('Time [s]','FontSize', FtSize, 'Interpreter','latex');
    ylabel('$\dot{T}$ [K/s]','fontsize', FtSize, 'Interpreter','latex');
    set(gcf,'color',[1 1 1]);
    set(gca,'fontsize', FtSize);
    

    % Flags to plot
    % figure(4)
    % plot( cE_1_a_vec ,'b','LineWidth',3)
    % hold on
    % plot( cE_1_b_vec ,'--r','LineWidth',3)
    % plot( cE_1_c_vec ,'+g','LineWidth',3)
    % plot( cE_1_d_vec ,'oy','LineWidth',3)
