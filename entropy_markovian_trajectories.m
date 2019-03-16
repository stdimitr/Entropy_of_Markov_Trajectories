function Hmarktraj = entropy_markovian_trajectories(symbseq)


%estimate the entropy reduction rate from a symbolic time serie

%INPUT: symbseq     : symbolic time series that describes the Markovian Chain
%OUTPUT: Hmarktraj  : the estimated H matrix of Markovian Trajectories

%DIMITRIADIS STAVROS v1.0 3/5/2012 / v1.1 14/6/2015

%%% This code is an implementation of entropy of Markov Trajectories
%%% original presented in:

%%% Laura Ekroot, Member,IEEE, and Thomas M. Cover, Fellow,IEEE.
%%% The Entropy of Markov Trajectories.
%%% IEEE TRANSACTIONS ON INFORMATIONTHEORY, VOL. 39, NO. 4, JULY 1993


%%% If you use this code, please cite the following article:
%%% Dimitriadis SI, Salis C.Mining Time-Resolved Functional Brain Graphs
%%%  to an EEG-Based Chronnectomic Brain Aged Index (CBAI).
%%% Front. Hum. Neurosci., 07 September 2017 | 
%%% https://doi.org/10.3389/fnhum.2017.00423

r1=unique(symbseq);



len=length(r1);
P=zeros(len,len);

for i=1:length(symbseq)-1
    symb1=find(r1==symbseq(i));
    symb2=find(r1==symbseq(i+1));
    P(symb1,symb2) = P(symb1,symb2) + 1; 
end


%normalize each row
for i=1:len
    sum1=sum(P(i,:));
    P(i,:)=P(i,:)/sum1;
end

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










