
```
function small_lms_test()
clear all; close all;
% make_pcm();
%near is micphone captured signal
fid=fopen('mixed.pcm', 'rb'); % Load far end
near_f=fread(fid,inf,'int16');
fclose(fid);

%far is speaker played music
fid=fopen('noise.pcm', 'rb'); % Load fnear end
far_f=fread(fid,inf,'int16');
fclose(fid);

len=length(near_f);
% fs = 16000;
mufb = 0.6;
threshold = 0.00001;
M = 16;
N = 64;
% L = M*N;

% step = 0.1875;
% w=zeros(L,1);
WFb=zeros(N+1,M);
YFb=zeros(N+1,M);
XFm=zeros(N+1,M);
% YFm=zeros(N+1,M);
NN=len;
Nb=floor(NN/N)-M;
start=1;
xo=zeros(N,1);
do=xo;
pn0=10*ones(N+1,1);
alp = 0.15; % Power estimation factor 
zm=zeros(N,1);

erfb=zeros(len,1);
ddecrfb=zeros(len,1);
aaa = zeros(len+Nb,1);
for kk=1:Nb
    pos = N * (kk-1) + start;

    % FD block method
    % ---------------------- Organize data

    %far is speaker played music
    xk = far_f(pos:pos+N-1);
    %near is micphone captured signal
    dk = near_f(pos:pos+N-1);

    %----------------------- far end signal process
    xx = [xo;xk];
    xo = xk;
    tmp = fft(xx);
    XX = tmp(1:N+1);

    dd = [do;dk]; % Overlap
    do = dk;
    tmp = fft(dd); % Frequency domain
    DD = tmp(1:N+1);

    % ------------------------far end Power estimation
    pn0 = (1 - alp) * pn0 + alp * real(XX.* conj(XX));
    pn = pn0;
  pn = (1 - alp) * pn + alp * M * pn0;

    % ---------------------- Filtering
    XFm(:,1) = XX;
    for mm=0:(M-1)
        m=mm+1;
        
        YFb(:,m) = XFm(:,m) .* WFb(:,m);
    end
    yfk = sum(YFb,2);
    tmp = [yfk ; flipud(conj(yfk(2:N)))];
    ykt = real(ifft(tmp));
    ykfb = ykt(end-N+1:end);

    tmp = [yfk ; zm];
    yktdjh = real(ifft(tmp));
    ykfbdjh = yktdjh(end-N+1:end);
    
    % ---------------------- Error estimation
    ekfb = dk - ykfb;
%     dddecrfb = dk - dk;
    ddecrfb(pos:pos+N-1) = ykfbdjh;
    erfb(pos:pos+N-1) = ekfb;
   
    tmp = fft([zm;ekfb]); % FD version for cancelling part (overlap-save)
    Ek = tmp(1:N+1);

    % ------------------------ Adaptation
    %Ek2 = Ek ./(M*pn + 0.001); % Normalized error
    Ek2 = Ek ./(pn + 0.001); % Normalized error
%     Ek2 = Ek; % not use pn 
    absEf = max(abs(Ek2), threshold);
    absEf = ones(N+1,1)*threshold./absEf;
    Ek2 = Ek2.*absEf;

    mEk = mufb.*Ek2;
    PP = conj(XFm).*(ones(M,1) * mEk')';
    tmp = [PP ; flipud(conj(PP(2:N,:)))];
    IFPP = real(ifft(tmp));
    PH = IFPP(1:N,:);
    tmp = fft([PH;zeros(N,M)]);
    FPH = tmp(1:N+1,:);
    WFb = WFb + FPH;
    aa(:,1) = WFb(:,1);
    aaa(pos:pos+N) = aa;
end

% plot(real(aaa));
% plot_task(near_f,far_f,erfb);
plot_task(near_f,'near_f',far_f,'far_f',erfb,'erfb');
% plot_task(near_f,'near_f',far_f,'far_f',ddecrfb,'dddecrfb');
writespeech('en345.pcm',erfb);
writespeech('yfb.pcm',ddecrfb);
return


function make_pcm()
fclear = readspeech('clear.pcm', 100000);
FL = length(fclear); 
noise=0.1*randn(1,FL);
fclearT=fclear-mean(fclear);                          
fclearT=fclearT/max(abs(fclearT));  
fclearT=fclearT';
mixed = fclearT+noise*0.5;
% plot_task(fclear,noise,mixed);
plot_task(fclear,'fclear',noise,'noise',mixed,'mixed');
writespeech('mixed.pcm',mixed*32767);
writespeech('noise.pcm',noise*32767);
return



function plot_task(s,s_title,exc,exc_title,s_rec,s_rec_title)
figure;
subplot(3,1,1);
stem(s/max(s),'Marker','none');
title(s_title);
xlabel('n (samples)');
ylabel('Amplitude');

subplot(3,1,2);
stem(exc/max(exc),'Marker','none');
title(exc_title);
xlabel('n (samples)');
ylabel('Amplitude');

subplot(3,1,3);
stem(s_rec/max(s_rec),'Marker','none');
title(s_rec_title);
xlabel('n (samples)');
ylabel('Amplitude');

return

function s = readspeech(filename, L)
    fid = fopen(filename, 'r');
    s = fread(fid, L, 'int16');
    fclose(fid);
return

function writespeech(filename,s)
    fid = fopen(filename,'w');
    fwrite(fid, s, 'int16');
    fclose(fid);
return
```
