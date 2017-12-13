rl=1:0.01:2;

figure;
subplot(1,2,1);
plot(r,in1s);
xlabel('r');
ylabel('g_{1,1}');
subplot(1,2,2);
plot(rl,U2(2)*30*ones(size(rl)));
axis([1,2,0,0.01]);
xlabel('r');
ylabel('homogeneous part');

figure;
subplot(1,2,1);
plot(r,in2s);
xlabel('r');
ylabel('g_{2,1}');
subplot(1,2,2);
plot(rl,-U2(2)*6./rl.^5);
xlabel('r');
ylabel('homogeneous part');

figure;
subplot(1,2,1);
plot(r,in3s);
xlabel('r');
ylabel('g_{3,1}');
subplot(1,2,2);
plot(rl,U2(2)*2./rl.^3+U3);
xlabel('r');
ylabel('homogeneous part');

figure;
subplot(1,2,1);
plot(r,in4s);
xlabel('r');
ylabel('g_{4,1}');
subplot(1,2,2);
plot(rl,-U2(2)./rl+U3*rl+U1(2));
xlabel('r');
ylabel('homogeneous part');

figure;
plot(r,sl(:,2));