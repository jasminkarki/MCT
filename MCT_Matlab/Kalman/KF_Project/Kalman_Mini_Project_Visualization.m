%% Plot results
% Plot the predicted states
k = 1:i;
figure(1);
subplot(211);
plot(k, X_prior_store(1,:),'b', 'LineWidth', 1.2)
hold on
plot(k, X_prior_store(2,:),'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('height');grid on;
legend('h1', 'h2','Location','southeast');

subplot(212);
plot(k, X_prior_store(3,:),'b', 'LineWidth', 1.2)
hold on
plot(k, X_prior_store(4,:),'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('height');grid on;
legend('h3','h4','Location','southeast');
sgtitle("X prior plot")
print('X prior Graph','-dpng');

% Plot the updated states 
figure(2);
subplot(211);
plot(k, X_post_store(1,:),'b', 'LineWidth', 1.2)
hold on
plot(k, X_post_store(2,:),'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('height');grid on;
legend('h1', 'h2','Location','southeast');

subplot(212);
plot(k, X_post_store(3,:),'b', 'LineWidth', 1.2)
hold on
plot(k, X_post_store(4,:),'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('height');grid on;
legend('h3','h4','Location','southeast');
sgtitle("X posterior plot")
print('X posterior Graph','-dpng');

% Plot predicted and estimated states
figure(3);
subplot(211);
plot(k, X_prior_store(1,:),'r', 'LineWidth', 1.2)
hold on
plot(k, X_prior_store(2,:),'g', 'LineWidth', 1.2)
hold on
plot(k, X_post_store(1,:),'b', 'LineWidth', 1.2)
hold on
plot(k, X_post_store(2,:),'m', 'LineWidth', 1.2)
xlabel('step'); ylabel('height');grid on;
legend('h1 prior', 'h2 prior','h1 post','h2 post','Location','southeast');

subplot(212);
plot(k, X_prior_store(3,:),'r', 'LineWidth', 1.2)
hold on
plot(k, X_prior_store(4,:),'g', 'LineWidth', 1.2)
hold on
plot(k, X_post_store(3,:),'b', 'LineWidth', 1.2)
hold on
plot(k, X_post_store(4,:),'m', 'LineWidth', 1.2)
xlabel('step'); ylabel('height');grid on;
legend('h3 prior', 'h4 prior','h3 post','h4 post','Location','southeast');
sgtitle("X prior and X posterior plot")
print('X prior post Graph','-dpng');


% Plot the covariance matrix 
figure(4);
subplot(311);
plot(k, Trace_P_prior,'b', 'LineWidth', 1.2)
xlabel('step'); ylabel('Trace of P');grid on;
legend('P prior','Location','southeast');
subplot(312);
plot(k, Trace_P_post,'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('Trace of P');grid on;
legend('P post','Location','southeast');
subplot(313);
plot(k, Trace_P_prior,'b', 'LineWidth', 1.2)
hold on
plot(k, Trace_P_post,'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('Trace of P');grid on;
legend('P prior','P post','Location','southeast');
sgtitle("Covariance Matrix plot")
print('Covariance Matrix Graph','-dpng');

% Plot the Kalman Gain 
figure(5);
plot(1:100, K_store(1,1:2:200), '--*', 'LineWidth', 1.2)
xlabel('step'); ylabel('Kalman gain');grid on;
sgtitle("Kalman Gain plot")
print('Kalman Gain Graph','-dpng');

% Plot the Innovation and Residual 
figure(6);
subplot(211);
plot(1:50, Innov_store(1,1:50),'b', 'LineWidth', 1.2)
hold on
plot(1:50, Residual_store(1,1:50),'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('X(1)');grid on;
legend('InnovX(1)','ResidualX(1)');
subplot(212);
plot(1:50, Innov_store(2,1:50),'b', 'LineWidth', 1.2)
hold on
plot(1:50, Residual_store(2,1:50),'r', 'LineWidth', 1.2)
xlabel('step'); ylabel('X(2)');grid on;
legend('InnovX(2)','ResidualX(2)');
sgtitle("Innovation and Residual plot")
print('InnovResidue Graph','-dpng');
