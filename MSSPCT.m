function  [Tc,w_CT,ff] =MSSPCT(x,fs,t,hlength,WL,num)


%% 
[xrow,xcol] = size(x);
Siglength=xrow;
% hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);
ht=ht';
%% Gaussian window
h =exp(-pi/0.32^2*ht.^2);%g
dh=-2*pi/(0.32^2)*ht.*exp(-pi/0.32^2*ht.^2);%g'
th=ht.*h;%tg
tth=ht.*ht.*h;%t.^2g
tdh=ht.*dh;%tg'
%%
Frelength=round(Siglength/2);
wopt2.fs=fs;
xi1h=-3.9265;
t1h=-0.1159;
wopt2.wp.t2h=-t1h;
wopt2.wp.t1h=t1h;
wopt2.wp.xi1h=xi1h;
wopt2.wp.xi2h=-xi1h;
t = (0:Siglength-1)/fs;%111
% M=551;%real_test
% dt = 0.002;
% m = 0:M-1;
% t = dt*m+0.7;
[N,Ratio] =PRE(x,fs,t,WL);

%%
    ff=0:fs/Siglength:(fs/2)-(fs/2)/Siglength;
    tfsupp=zeros(N,Siglength);
    [Spectemp,f]=Polychirplet_l(x,fs,hlength,Ratio,h);%Polychirplet2(Sig,SampFreq,R1,N1,h);Polychirplet(x,fs,Ratio,Siglength,hlength,h)
    test=ecurve(Spectemp,ff,wopt2,'method',1);%脊线提取
    tfsupp(1,:)=test(1,:);
    [p, z] = polylsqr(t,test(1,:),N);
    Ratio=z(2:end);
 %%
%     i=2;
% while N>1
for i=2:N
    [Spectemp,f]=Polychirplet_l(x,fs,hlength,Ratio,h);%Polychirplet2(Sig,SampFreq,R1,N1,h);Polychirplet(x,fs,Ratio,Siglength,hlength,h)
    test2=ecurve(Spectemp,ff,wopt2,'method',1);%脊线提取
    tfsupp(i,:)=test2(1,:);
    [p, z] = polylsqr(t,test2(1,:),N);
    Ratio=z(2:end);
%     er=max((tfsupp(i,:)-tfsupp(i-1,:))/tfsupp(i-1,:));
%     i=i+1;
%     if er<1e-4
%         break;
%     end        
end 
    Spec = Spectemp;
    Spec_d = Polychirplet_l(x,fs,hlength,Ratio,dh);
    Spec_th = Polychirplet_l(x,fs,hlength,Ratio,th);
    Spec_tth = Polychirplet_l(x,fs,hlength,Ratio,tth);
    Spec_tdh = Polychirplet_l(x,fs,hlength,Ratio,tdh);

%%
www = ( Spec_tth.*Spec_d-Spec_th.*Spec_tdh-Spec_th.*Spec ) ./ ( Spec_th.^2-Spec.*Spec_tth );
w_CT=repmat((1:Frelength).',1,Siglength)+round(imag(www));
% www = -Spec_th.*Spec ./ ( Spec_th.^2-Spec.*Spec_tth );
% w_CT=repmat((1:Frelength).',1,Siglength)+round(imag(www));
%%
for ite=1:num
[Tc]=SST(Spec,w_CT);
Spec=Tc;
end
% Tc=Tc/(xrow/2);
end

function [Ts_f]=SST(tfr_f,omega_f)
[tfrm,tfrn]=size(tfr_f);
Ts_f= zeros (tfrm,tfrn) ; 
for b=1:tfrn%time
    % Reassignment step
    for eta=1:tfrm%frequency
        %if abs(tfr_f(eta,b))>0.001*mx%you can set much lower value than this.
            k = omega_f(eta,b);
            if k>=1 && k<=tfrm
%             if k==eta
                Ts_f(k,b) = Ts_f(k,b) + tfr_f(eta,b);
            end
        %end
    end
end
end

function  [N,Ratio] =PRE(x,fs,t,WL)
    Sig=x;
    wopt2.fs=fs;
    xi1h=-3.9265;
    t1h=-0.1159;
    wopt2.wp.t2h=-t1h;
    wopt2.wp.t1h=t1h;
    wopt2.wp.xi1h=xi1h;
    wopt2.wp.xi2h=-xi1h;
    WinLen=WL;
    WinLen=WinLen+1-rem(WinLen,2);
    tt = linspace(-0.5,0.5,WinLen)';
    h = exp(-pi/0.32^2*tt.^2);
    N11=length(Sig);
    ff1=0:fs/N11:(fs/2)-(fs/2)/N11;
    Ratio = [0,0];
    [pct1_matrix_chirp, fre] =Polychirplet_l(Sig',fs,WinLen,Ratio,h);%PCT
    tfsupp=ecurve(pct1_matrix_chirp,ff1',wopt2,'method',1);%脊线提取
    [p, z] = polylsqr(t,tfsupp(1,:),21);
    Ratio=z(2:end);
    N=length(Ratio);
end
