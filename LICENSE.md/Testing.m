# Geometric-Clustering-based-on-Information-Bottleneck
clear
format compact
close all

%% 
% Definition of mu's and Sigma
% Mean vectors and covariance matrix
m1=[-1 1]';  m2=[1 1]';  m3=[1 -1]';  m4=[-1 -1]';
S=[.3 0; 0 .3];
% Number of data points
n_points_per_class=625;

% (i) Data point generation
X=[mvnrnd(m1',S,n_points_per_class); mvnrnd(m2',S,n_points_per_class);...
    mvnrnd(m3',S,n_points_per_class);mvnrnd(m4',S,n_points_per_class)]';
label=[ones(1,n_points_per_class) 2*ones(1,n_points_per_class)...
    3*ones(1,n_points_per_class) 4*ones(1,n_points_per_class)];
[l,N]=size(X);
%% Initialize
beta = 10;
%pIX = eye(N)/sum(N(:));
pI = 1/N*ones(N,1);
unnorm_p_X_I = exp((-1/(2)) * distance_matrix(X));
Z_X_I = repmat(sum(unnorm_p_X_I), [N 1]);
p_X_I = unnorm_p_X_I ./ Z_X_I;
pIX = p_X_I' .* repmat(pI, [1 N]);
XC0 = repmat(mean(X,2), [1 4]) + rand(2,4);

my_distance = zeros(N,4);
for i = 1:N
    for j = 1:4
        my_distance(i,j) = mydistance(X(:,i),XC0(:,j));
    end
end

unnorm_p0_X_C = exp(-1*my_distance);
Z0_C_beta = repmat(sum(unnorm_p0_X_C), [N 1]);
p0X_C = unnorm_p0_X_C ./ Z0_C_beta;
p0C = ones(4,1)*0.25;
%% Implement information bottleneck to clustering

maxiter = 200;
tol = 0.00000001;
[ pC_I, pX_C, PC, XC ] = Geo_all_iteration( pIX, beta, p0X_C,p0C,X,maxiter, tol );
% [pC_I, pX_C, PC,XC,L] = Geo_IB_per_iteration(pIX, beta, p0X_C,p0C,X);
%  for i = 1:200
%      pX_C_pre = pX_C;
%      pC_pre = PC;
%      %XC_pre = XC;
%      [pC_I, pX_C, PC, XC,L] = Geo_IB_per_iteration(pIX, beta, pX_C_pre,pC_pre,X);
%  end
 
 %% Verify its significance
 
 hat_label = zeros(N,1);
 for i = 1:N
     hat_label(i) = find(pC_I(:,i) == max(pC_I(:,i)));
 end
 
 %% Visualize
 label_1 = repmat({'1'}, [n_points_per_class 1]);
  label_2 = repmat({'2'}, [n_points_per_class 1]);
   label_3 = repmat({'3'}, [n_points_per_class 1]);
    label_4 = repmat({'4'}, [n_points_per_class 1]);
    marker_color = [1 0 0; .7 .1 0; .5 .5 .5; .1 .8 .2];
 figure
 for i = 1:N
     if hat_label(i) == 1
         plot(X(1,i),X(2,i),'.','MarkerSize',10','Color',marker_color(1,:));
         hold on
     end
     if hat_label(i) == 2
         plot(X(1,i),X(2,i),'.','MarkerSize',10','Color',marker_color(2,:));
         hold on
     end
     if hat_label(i) == 3
         plot(X(1,i),X(2,i),'.','MarkerSize',10','Color',marker_color(3,:));
         hold on
     end
     if hat_label(i) == 4
         plot(X(1,i),X(2,i),'.','MarkerSize',10','Color',marker_color(4,:));
         hold on
     end
 end
 





