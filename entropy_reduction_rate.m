function entredrate=entropy_reduction_rate(symbseq)


%estimate the entropy reduction rate from a symbolic time serie

%INPUT: symbolic time serie
%OUTPUT: the estimated entropy reduction rate

%DIMITRIADIS STAVROS 3/2012

r1=unique(symbseq);



len=length(r1);
conprob=zeros(len,len);

for i=1:length(symbseq)-1
    symb1=find(r1==symbseq(i));
    symb2=find(r1==symbseq(i+1));
    conprob(symb1,symb2) = conprob(symb1,symb2) + 1; 
end



%normalize each row
for i=1:len
    sum1=sum(conprob(i,:));
    conprob(i,:)=conprob(i,:)/sum1;
end


%estimate entropy
entropy=0;
prob=zeros(1,len);

for i=1:len
    p=find(symbseq==r1(i));
    prob(i)=length(p)/length(symbseq);
    entropy=entropy + prob(i)*log(prob(i));
end

entropy=-entropy;

    

condentropy=zeros(1,len);
%estimate conditional entropy

sum1=0;
for i=1:len
    [r c]=find(conprob(i,:) > 0);
    for j=1:length(c)
        sum1 = sum1 + conprob(i,c(j))*log(conprob(i,c(j)));
    end
    condentropy(i)= - sum1;
    sum1=0;
end


entredrate=0;
sum1=0;
for i=1:len
    sum1 = sum1 + prob(i)*condentropy(i);
end
    entredrate = (entropy - sum1)/entropy;
    
    
    
    
    