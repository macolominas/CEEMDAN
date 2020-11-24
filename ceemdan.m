%Function for CEEMDAN

%WARNING: for this code works it is necessary to include in the same
%directoy the file emd.m developed by Rilling and Flandrin.
%This file is available at %http://perso.ens-lyon.fr/patrick.flandrin/emd.html
%We use the default stopping criterion.
%We use the last modification: 3.2007

%   Syntax

%modes=ceemdan(x,Nstd,NR,MaxIter,SNRFlag)
%[modes its]=ceemdan(x,Nstd,NR,MaxIter,SNRFlag)

%   Description

%OUTPUT
%modes: contain the obtained modes in a matrix with the rows being the modes        
%its: contain the sifting iterations needed for each mode for each realization (one row for each realization)

%INPUT
%x: signal to decompose
%Nstd: noise standard deviation
%NR: number of realizations
%MaxIter: maximum number of sifting iterations allowed.
%SNRFlag: if equals 1, then the SNR increases for every stage, as in [1].
%           If equals 2, then the SNR is the same for all stages, as in [2]. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The current is an improved version, introduced in:

%[1] Colominas MA, Schlotthauer G, Torres ME. "Improve complete ensemble EMD: A suitable tool for biomedical signal processing" 
%       Biomedical Signal Processing and Control vol. 14 pp. 19-29 (2014)

%The CEEMDAN algorithm was first introduced at ICASSP 2011, Prague, Czech Republic

%The authors will be thankful if the users of this code reference the work
%where the algorithm was first presented:

%[2] Torres ME, Colominas MA, Schlotthauer G, Flandrin P. "A Complete Ensemble Empirical Mode Decomposition with Adaptive Noise"
%       Proc. 36th Int. Conf. on Acoustics, Speech and Signa Processing ICASSP 2011 (May 22-27, Prague, Czech Republic)

%Author: Marcelo A. Colominas
%contact: macolominas@bioingenieria.edu.ar
%Last version: 25 feb 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [modes,its]=ceemdan(x,Nstd,NR,MaxIter,SNRFlag)
x=x(:)';
desvio_x=std(x);
x=x/desvio_x;

modes=zeros(size(x));
temp=zeros(size(x));
aux=zeros(size(x));
iter=zeros(NR,round(log2(length(x))+5));

for i=1:NR
    white_noise{i}=randn(size(x));%creates the noise realizations
end;

for i=1:NR
    modes_white_noise{i}=emd(white_noise{i});%calculates the modes of white gaussian noise
end;

for i=1:NR %calculates the first mode
    xi=x+Nstd*modes_white_noise{i}(1,:)/std(modes_white_noise{i}(1,:));
    [temp, o, it]=emd(xi,'MAXMODES',1,'MAXITERATIONS',MaxIter);
    temp=temp(1,:);
    aux=aux+(xi-temp)/NR;
    iter(i,1)=it;
end;

modes= x-aux; %saves the first mode
medias = aux;
k=1;
aux=zeros(size(x));
es_imf = min(size(emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter)));

while es_imf>1 %calculates the rest of the modes
    for i=1:NR
        tamanio=size(modes_white_noise{i});
        if tamanio(1)>=k+1
            noise=modes_white_noise{i}(k+1,:);
            if SNRFlag == 2
                noise=noise/std(noise); %adjust the std of the noise
            end;
            noise=Nstd*noise;
            try
                [temp,o,it]=emd(medias(end,:)+std(medias(end,:))*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
            catch    
                it=0;
                temp=emd(medias(end,:)+std(medias(end,:))*noise,'MAXMODES',1,'MAXITERATIONS',MaxIter);
            end;
            temp=temp(end,:);
        else
            try
                [temp, o, it]=emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter);
            catch
                temp=emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter);
                it=0;
            end;
            temp=temp(end,:);
        end;
        aux=aux+temp/NR;
    iter(i,k+1)=it;    
    end;
    modes=[modes;medias(end,:)-aux];
    medias = [medias;aux];
    aux=zeros(size(x));
    k=k+1;
    es_imf = min(size(emd(medias(end,:),'MAXMODES',1,'MAXITERATIONS',MaxIter)));
end;
modes = [modes;medias(end,:)];
modes=modes*desvio_x;
its=iter;


   