
% N.O. November 08
% GPCA matrix rank minimization using theory of moments

% modified on 17 Jan 01

% reweight by using adaptive step size, choose delta according to minimum
% singular varlue of the last solution

% solve for a valid noise via polynomial root finding

%l_inf norm bounded noise

clear all;
%close all;


addpath('../CodeGPCAPDASpectral/helper_functions')

% fix random number generator so we can repeat the experiment
seed = 1; %5;
rng(seed)

%profile on
% ARX Systems:

% I.   y_t = 0.2y_(t-1)+0.24y_(t-2)+2u_(t-1)+ w_t 
% poles at -0.4, 0.6
% II.  y_t = -1.4y_(t-1)-0.53y_(t-2)+u_(t-1)+w_t
% poles at -0.7+0.2i, -0.7-0.2i
% III. y_t = 1.7y_(t-1)-0.72y_(t-2)+0.5u_(t-1)+w_t
% poles at 0.9, 0.8

% Hybrid system: 
% y(t) = p1(t)*sys1+p2(t)*sys2+p3(t)*sys3+w(t)

% -2<p1<2, -1<p2<1, -3<p3<3 for all t

sys1_a1 = 1; sys1_a2 = 0; sys1_b1 = 0;
sys2_a1 = 0; sys2_a2 = 1; sys2_b1 = 0;
sys3_a1 = 0; sys3_a2 = 0; sys3_b1 = 1;

N = 96; % time horizon
%N = 60;
rand_switch = 0;
sw_seq = 1;
num_sys = 3;

epsilon_eta = 0.25;

if sw_seq == 1
p1 = [0, 0, 0.2*ones(1,N/3), -1.4*ones(1,N/3), 1.7*ones(1,N/3) ];
p2 = [0, 0, 0.24*ones(1,N/3), -0.53*ones(1,N/3), -0.72*ones(1,N/3) ];
p3 = [0, 0, 2*ones(1,N/3), 1*ones(1,N/3), 0.5*ones(1,N/3) ];
end

if sw_seq == 2
p1 = [0, 0, 0.2*ones(1,N/4), -1.4*ones(1,N/4), 0.2*ones(1,N/4), -1.4*ones(1,N/4) ];
p2 = [0, 0, 0.24*ones(1,N/4), -0.53*ones(1,N/4), 0.24*ones(1,N/4), -0.53*ones(1,N/4) ];
p3 = [0, 0, 2*ones(1,N/4), 1*ones(1,N/4), 2*ones(1,N/4), 1*ones(1,N/4) ];
end

if sw_seq == 3
p1 = [0, 0, 0.2*ones(1,N)];
p2 = [0, 0, 0.24*ones(1,N)];
p3 = [0, 0, 2*ones(1,N)];
end

if sw_seq == 4
p1 = [0, 0, 1.7*ones(1,N/4), -1.4*ones(1,N/4), 1.7*ones(1,N/4), -1.4*ones(1,N/4) ];
p2 = [0, 0, -0.72*ones(1,N/4), -0.53*ones(1,N/4), -0.72*ones(1,N/4), -0.53*ones(1,N/4) ];
p3 = [0, 0, 0.5*ones(1,N/4), 1*ones(1,N/4), 0.5*ones(1,N/4), 1*ones(1,N/4) ];
end

if sw_seq == 5
p1 = [0, 0, 1.7*ones(1,N/4), 0.2*ones(1,N/4), 1.7*ones(1,N/4), 0.2*ones(1,N/4) ];
p2 = [0, 0, -0.72*ones(1,N/4), 0.24*ones(1,N/4), -0.72*ones(1,N/4), 0.24*ones(1,N/4) ];
p3 = [0, 0, 0.5*ones(1,N/4), 2*ones(1,N/4), 0.5*ones(1,N/4), 2*ones(1,N/4) ];
end

if sw_seq == 6
p1 = [0, 0, -1.4*ones(1,N/3), 0.2*ones(1,N/3), -1.4*ones(1,N/3) ];
p2 = [0, 0, -0.53*ones(1,N/3), 0.24*ones(1,N/3), -0.53*ones(1,N/3) ];
p3 = [0, 0, 1*ones(1,N/3), 2*ones(1,N/3), 1*ones(1,N/3) ];
end

if sw_seq == 7
p1 = [0, 0, 0.2*ones(1,N/3), -1.4*ones(1,N/3), 0.2*ones(1,N/3) ];
p2 = [0, 0, 0.24*ones(1,N/3), -0.53*ones(1,N/3), 0.24*ones(1,N/3) ];
p3 = [0, 0, 2*ones(1,N/3), 1*ones(1,N/3), 2*ones(1,N/3) ];
end

if rand_switch
order = [1, 2, randperm(N)+2];
p1 = p1(order);
p2 = p2(order);
p3 = p3(order);
end


w = [0, 0, (rand(1,N)-0.5)*2*epsilon_eta];

%w = [0, 0, 0.2*ones(1,N)];

u = [0, 0, randn(1,N)];

m = 3; %order of the system


% generate data
y(1) = 0;
y(2) = 0;
for t = 3:N+2
    s1 = sys1_a1*y(t-1)+sys1_a2*y(t-2)+sys1_b1*u(t-1);
    s2 = sys2_a1*y(t-1)+sys2_a2*y(t-2)+sys2_b1*u(t-1);
    s3 = sys3_a1*y(t-1)+sys3_a2*y(t-2)+sys3_b1*u(t-1);
    y(t) = p1(t)*s1+p2(t)*s2+p3(t)*s3+w(t); %process noise
end

Hy = toeplitz(y(4:-1:2),y(4:N+1)); 
Hu = u(3:N);
H = [Hy;Hu];

% for i = 3 :N+2
%     data(i-2,:) = [y(i) y(i-1) y(i-2) u(i-1)];
% end

[vermap,powers] = veronese(H,num_sys);
rank(vermap)
[u_n, d_n, v_n] =svd(vermap);
d_n = diag(d_n)';

reg_size = m+1;
[row, col] = size(vermap);

mom_coef = cell(1,num_sys);
for i = 1:num_sys
    mom_coef{i} = zeros(size(vermap));
end

for i = 1:row
    for j = 1:num_sys
        if powers(i,1)==j %need to change y(t)^j with [y(t)-w_t]^j
            pow_now = powers(i,2:end);
            pow_now = pow_now'*ones(1,col);
            coefs_from_other_yu = prod(H(2:end,:).^pow_now,1);
            yt_wt_coefs = pascal(j+1,1);
            yt_wt_coefs = yt_wt_coefs(j+1,1:(j+1));
            for k = 1:j
                yt_power = j-k;
                coefs_all = coefs_from_other_yu.*(H(1,:).^yt_power);
                mom_coef{k}(i,:) = mom_coef{k}(i,:)+ coefs_all*yt_wt_coefs(k+1);
            end
        end
    end
end

%%
%Debug: check if moment coefficient matrices are formed correctly 
clean_y = y-w;
H_new = H;
H_new(1,:) = clean_y(4:N+1); 
[vermap_c, powers] = veronese(H_new,num_sys);
rank(vermap_c)
[u_c, d_c, v_c] =svd(vermap_c);
d_c = diag(d_c)';
% %%
% 
% try_err = ones(row,1)*w(4:4+col-1);
% m1s = mom_coef{1}.*try_err;
% m2s = mom_coef{2}.*(try_err.^2);
% 
% clean_mat = vermap+m1s+m2s;
% rank(clean_mat)
% [u d v] =svd(clean_mat);
% u
% th = diag(d)'
% sum(th)
%%

seq_size = num_sys; %if this is 2n, this means using nXn moment matrices
half_size = floor(seq_size/2);
delta = 10^-4;
Wy = eye(row);
Wz = eye(col);
nnzs = [];
min_sv = zeros(1,25);
null_err_norm = zeros(1,25);
matrix_dist = zeros(1,25);

% cvx_solver sdpt3
cvx_solver sedumi

cvx_quiet(true);
% Turn off the weird warning issued by hankel() function.
warning('off','MATLAB:hankel:AntiDiagonalConflict');
for iter = 1:25 %reweighting iterations
  cvx_quiet(true);
  cvx_begin sdp
    variable Y(row,row) symmetric
    variable Z(col,col) symmetric
    variable Mij(seq_size,col)
    X = vermap;
    for l = 1:num_sys
        X = X + mom_coef{l}.*(ones(row,1)*Mij(l,:));
    end

    minimize( trace(Wy*Y)+trace(Wz*Z))
    subject to
        [Y X; X' Z] >= 0;
        if ~mod(seq_size,2) %even conditions A(k)>=0 \eps^2A(k-1)>=C
            for k = 1:col
                hankel([1;Mij(1:half_size,k)],Mij(half_size:seq_size,k))>=0;
                if half_size == 1
                    epsilon_eta^2-Mij(2,k)>=0;
                else
                    epsilon_eta^2*hankel([1;Mij(1:half_size-1,k)],Mij(half_size-1:seq_size-2,k))-hankel(Mij(2:half_size+1,k),Mij(half_size+1:seq_size,k))>=0;
                end
            end
        else %odd conditions \eps A(k)>=B(k)>= -\eps A(k)
            for k = 1:col
                epsilon_eta*hankel([1;Mij(1:half_size,k)],Mij(half_size:seq_size-1,k)) - hankel(Mij(1:half_size+1,k),Mij(half_size+1:seq_size,k)) >= 0;
                hankel(Mij(1:half_size+1,k),Mij(half_size+1:seq_size,k)) + epsilon_eta*hankel([1;Mij(1:half_size,k)],Mij(half_size:seq_size-1,k)) >= 0;
            end
        end 
    cvx_end
  cvx_quiet(false);
  
  
  if sum(isnan(X(:)))>0; X = X_old; Mij = Mij_old; 
      break;
  else
      X_old = X; Mij_old = Mij;
  end
  [u, d, v] =svd(X);
  d=diag(d);
  min_sv(iter) = d(end);
  sign_c = sign(u_c(1,end));
  sign_new = sign(u(1,end));
  if sign_c~=sign_new; u(:,end)=-u(:,end); end
  null_err_norm(iter) = norm(u(:,end)'-u_c(:,end)');
  matrix_dist(iter) = norm(vermap_c-X);
  if min_sv(iter)<=10^-8; break; end; %small enough
  %Mij(1,:)
  %w
  %pause
  Wy = inv(Y+max(delta,d(end))*eye(row));
  Wz = inv(Z+max(delta,d(end))*eye(col));
  
end
warning('on','MATLAB:hankel:AntiDiagonalConflict');
%cvx_quiet(false);
figure; plot(min_sv,'*'); title('The change in min Singular Value in Log-det iterations');
%figure; plot(null_err_norm,'*'); title('The error between the true null space and estimated null space');
sign_c = sign(u_c(1,end));
sign_n = sign(u_n(1,end));
if sign_c~=sign_n; u_n(:,end)=-u_n(:,end); end
%hold on;plot(norm(u_n(:,end)'-u_c(:,end)')*ones(1,iter),'r*');
rank(X);
[u, d, v] =svd(X);
d = diag(d)';
%%
% prob = 1-(Mij(1,:).^2)./Mij(2,:); %prob of zero
% point_mass = Mij(2,:)./Mij(1,:);
%%
if 0
    for i=1:col
        figure; hold on;
        stem(Mij(1,i),1,'linewidth',2);
        stem(w(3+i),1,'r','linewidth',2);
        axis([-epsilon_eta epsilon_eta 0 1]);
    end
end

%%
%replace mean

null_est = u(:,end);

poly_coefs = zeros(num_sys+1,col);
poly_coefs(num_sys+1,:) = null_est'*vermap;
for i = 1:num_sys
    poly_coefs(i,:) =  null_est'*mom_coef{num_sys+1-i};
end

%%
noise_est = zeros(1, col);
for i =1:col
    dummy = roots(poly_coefs(:,i));
    [xx, dum_ind] = min(abs(dummy));
    noise_est(i) = real(dummy(dum_ind));
end
%%

%     X = vermap;
%     for l = 1:num_sys
%         X = X + mom_coef{l}.*(ones(row,1)*Mij(l,:));
%     end

clean_y_est = y(4:col+3)-noise_est;
H_new_est = H;
H_new_est(1,:) = clean_y_est; 
[vermap_est,powers] = veronese(H_new_est,num_sys);
rank(vermap_est)
[u_est, d_est, v_est] =svd(vermap_est);
d_est = diag(d_est)';

c = u_est(:,end);
x = H_new_est;
vermap = vermap_est';

[b_moment, group_moment] = construct_models(c,powers,x,vermap,num_sys);

figure; plot(group_moment,'*'); title('Moment Clustering')

b_m = [];
for i = 1:num_sys
    b_m = [b_m, b_moment(:,i)/b_moment(1,i)];
end
max_error_moment = max(min(abs(b_m'*H_new_est)));
disp(['Max Error of Moment Method is ' num2str(max_error_moment)]);

c = u_n(:,end);
x = H;
[vermap,powers] = veronese(H,num_sys);
vermap = vermap';

[b_gpca, group_gpca] = construct_models(c,powers,x,vermap,num_sys);
figure; plot(group_gpca,'*'); title('GPCA Clustering');

b_g = [];
for i=1:num_sys
    b_g = [b_g, b_gpca(:,i)/b_gpca(1,i)];
end
max_error_gpca = max(min(abs(b_g'*H)));
disp(['Max Error of GPCA is ' num2str(max_error_gpca)]);