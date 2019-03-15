


%%% Laura Ekroot, Member,IEEE, and Thomas M. Cover, Fellow,IEEE.
%%% The Entropy of Markov Trajectories.
%%% IEEE TRANSACTIONS ON INFORMATIONTHEORY, VOL. 39, NO. 4, JULY 1993


% You can run my MATLAB code with the example given in the original article
% 
P = [ 0    0.9000         0    0.1000 ;
    0.2500    0.2500    0.2500    0.2500 ;
    0.5000    0.5000         0         0 ; 
         0         0    1.0000         0];

%%%%%%%% STATIONARY DISTRIBUTION
m=zeros(1,len);

[V D] = eig( P.' );
st = V(:,1).';
m=abs(st)./sum(abs(st));


%%% THE ENTROPY RATE
H=0;

sum2=0;
for k=1:len
    sum1=zeros(1,len);
    for l=1:len
        sum1(l)=P(k,l)*log2(P(k,l));
    end
        sum2 = sum2 + m(k)*nansum(sum1);
end

H = - sum2;


%% THE MATRIX OF THE FIRST STEP ENTROPY
H1=zeros(len,len);

sum2=zeros(1,len);
for k=1:len
    sum1=zeros(1,len);
    c=find(P(k,:) > 0);
    for l=1:length(c)
        sum1(l) = -P(k,c(l))*log2(P(k,c(l)));
    end
    sum2(1,k)= nansum(sum1);
end


for k=1:len
    H1(k,:)=sum2(k);
end

%%%%% CREATE MATRIX HDelta
Hdelta=zeros(len,len);


%% We first have to estimate the entropy rate H(X) - equation 2 
vec=zeros(1,len);
vec=H./m;

for k=1:len
    Hdelta(k,k)=vec(k);
end

%%% ESTIMATE K - equation 28
K=zeros(len,len);
I=eye(len);
A=zeros(len,len);

for k=1:len
    A(k,:)=m;
end

K=((I - P + A)^(-1))*(H1 - Hdelta);

%%% Estimate K' equation 29
K1=zeros(len,len);

for k=1:len
    K1(:,k)=K(k,k);
end

%%%%%%%%%%% equation 27 - estimate the matrix H of trajectory entropies
Hmarktraj=zeros(len,len);
Hmarktraj = K - K1 + Hdelta ;

%%%% PLOT THE MATRICS
Hdelta
K
K1
Hmarktraj





