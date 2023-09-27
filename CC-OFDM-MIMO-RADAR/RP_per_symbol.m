%% 19 Januari 2023
clear all;clc;close all;
load Data.mat 'teta_c' 'd' 'c' 'teta' 'Aptayl' 'lamda' 'Nt';
Nc = 2048;      % Number of sub carriers
Nb=63;          % Number of beam Nb = 63 beams
Nsym=12;        % Numer of symbol
Nsb = 8;        % Number of sub bands
Ncb = Nc/Nsb;   % Number of sub carriers per beam
Nsb_beam = 3;   
Nsym_sb = 4;
Nsamp = 1e4;%8*Nc;                  % Number of samples in 1 symbol period
samp=1;
g=golay(Ncb);
teta_c = teta_c(2:Nb+1); 
bim = 1:Nb;
% p=0:Nt-1;
p = -(Nt-1)/2:(Nt-1)/2;
f_subc = 1:256;
T0=0.1E-3;         % Symbol duration in seconds
delf=1/T0;       % sub carrier spacing in Hz
fif=20E6;           
fup=3E9-fif; 
% fup=3e9-f0;
fn=fif+(0:Nc-1).'*delf;
nftb=Ncb/2;
fnb = reshape(fn,Ncb,[]);
fnb_c = fnb(Ncb/2,:); 
Att_max = 1./((15e3:15e3:180e3).^4);
dmin = 2*(0.5*lamda*Nt)^2/lamda;
%% initialization
AFb_tapered_iso = zeros(Ncb,length(teta),Nb);
AFb_tapered_cos = zeros(Ncb,length(teta),Nb);
AFb_tapered_uni_iso = zeros(Ncb,length(teta),Nb);
AFb_tapered_uni_cos = zeros(Ncb,length(teta),Nb);
fnb_int = zeros(Ncb,Nb,Nsym);
sb63_sym = zeros(Nb,Nsym);
f_band = zeros(Ncb,Nb,Nsym);
alfa = zeros(Ncb,length(teta),Nb);
beta = zeros(Ncb,Nb,Nsym);
beta_n = zeros(Ncb,Nb,Nsym);
beta_corr = zeros(Ncb,Nb,Nsym);
delt_n = zeros(Ncb,Nb,Nsym);
f_band_c = zeros(Nb,Nsym);
alfa_c = zeros(Nb,Nsym);
beta_c = zeros(Nb,Nsym);
delt_bc = zeros(Nb,Nsym);
AFb_tay_cod_corr = zeros(Ncb,length(teta),Nt);
%%
teta_rnd = round(teta*180/pi*10)/10;
teta_c_rnd = round(teta_c*180/pi*10)/10;
% for itc = 1:Nb
%     aa(itc) = find(teta_rnd==teta_c_rnd(itc));
% end

sb63 = sb_beam_schedule(Nc,Ncb,Nb,Nsb,Nsb_beam);
gb = [1:Ncb/2;Ncb/2+1:Ncb];
gb_sb1 = reshape(gb,8,[]).';
gb1 = circshift(gb,4,2);
gb_sb2 = reshape(gb1,8,[]).';
gb2 = circshift(gb,8,2);
gb_sb3 = reshape(gb2,8,[]).';
sb_beam_I = [1:8:Nb;5:8:Nb;2:8:Nb;6:8:Nb;3:8:Nb;7:8:Nb;4:8:Nb;8:8:Nb+1;]; % ommit beam 64 in sub band 8
sb_beam_II = [7:8:Nb;3:8:Nb;8:8:Nb+1;4:8:Nb;1:8:Nb;5:8:Nb;2:8:Nb;6:8:Nb]; % ommit beam 64 in sub band 3
sb_beam_III = [5:8:Nb;1:8:Nb;6:8:Nb;2:8:Nb;7:8:Nb;3:8:Nb;8:8:Nb+1;4:8:Nb];% ommit beam 64 in sub band 7
sb_beam (:,:,1)=sb_beam_I;
sb_beam (:,:,2)=sb_beam_II;
sb_beam (:,:,3)=sb_beam_III;

gb_sym=zeros(Nb+1,Nsym);
imod = zeros(Nsym);
for isym = 1:Nsym
    if isym<5
        sb_beam_i = sb_beam_I;
        gb_sb_i = gb_sb1;
    elseif and(isym>4,isym<9)
        sb_beam_i = sb_beam_II;
        gb_sb_i = gb_sb2;
    elseif isym>8
        sb_beam_i = sb_beam_III;
        gb_sb_i = gb_sb3;
    end
    for isb = 1: Nsb
        imod (isym)= mod(isym,4);
        if imod (isym) == 0 
            imod (isym) = 4;
        end
        gb_sym(sb_beam_i(isb,:),isym)=gb_sb_i((imod(isym)-1)*Nsb+isb,:);
    end
end

for ib = 1:Nb
    for isym = 1:Nsym
        sb63_sym (ib,:)= reshape(ones(Nsym_sb,1)*sb63(ib,:),1,[]);              % 12 symbols with 3 sub band per beam
    end
end
for i = 1:Nb
    fieldname = sprintf('beam_%i', i);
    Tx.(fieldname).sb_sym = sb63_sym(i,:);
    Tx.(fieldname).gb_sym = gb_sym(i,:);
    Tx.(fieldname).teta_c = teta_c(i);
%     Tx.(fieldname).teta3dBu = teta3dBu(i);
end

%% ARRAY PARAMETERS
for isym = 1:Nsym
    f_band(:,:,isym) = fnb(:,sb63_sym(:,isym));
    f_band_c(:,isym)= fnb_c(:,sb63_sym(:,isym));
    alfa_c(:,isym)= d/c.*(f_band_c(:,isym)+fup)*2*pi.*sin(teta_c(1:Nb).');  %% teta= negative 
    beta_c(:,isym)= -alfa_c(:,isym);
    delt_bc(:,isym) = beta_c(:,isym)./(2*pi.*f_band_c(:,isym));
    symbol = sprintf('symbol_%i', isym);
    Tx.(symbol).alfa_c = alfa_c(:,isym);
    Tx.(symbol).beta_c = beta_c(:,isym);
    Tx.(symbol).delt_bc = delt_bc(:,isym);

%% WAVEFORM PARAMETERS
    for inb = [1 2 9 23 24 25 31 32 33 39 40 41 55 62 63]%1:17%Nb
        beam = sprintf('beam_%i', inb);
%         delt_n(:,inb,isym) = f_band_c(inb,isym)./f_band(:,inb,isym).*delt_bc(inb,isym);    
        delt_n(:,inb,isym) = (f_band_c(inb,isym)./f_band(:,inb,isym)).*delt_bc(inb,isym);    
        beta_n(:,inb,isym) = 2*pi*f_band(:,inb,isym).*delt_n(:,inb,isym);
        beta(:,inb,isym)=2*pi*f_band(:,inb,isym).*delt_bc(inb,isym);  
        alfa(:,:,inb)=d/c.*(f_band(:,inb,isym)+fup)*2*pi*sin(teta);
        beta_corr (:,inb,isym)= 2*pi.*f_band(:,inb,isym).*(delt_n(:,inb,isym)-Tx.(symbol).delt_bc(inb));
        Tx.(symbol).(beam).f_band = f_band(:,inb,isym);
        Tx.(symbol).(beam).f_band_c = f_band_c(inb,isym);
        Tx.(symbol).(beam).delt_bc = delt_bc(inb,isym);
        Tx.(symbol).(beam).delt_n = delt_n(:,inb,isym);
        Tx.(symbol).(beam).beta_c = beta_c(inb,isym);
        Tx.(symbol).(beam).beta_n = beta_n(:,inb,isym);
        Tx.(symbol).(beam).beta = beta(:,inb,isym);
        Tx.(symbol).(beam).beta_corr = beta_n(:,inb,isym)-beta(:,inb,isym);
        Tx.(symbol).(beam).code = g(gb_sym(inb,isym),:);
        Tx.(symbol).(beam).alfa = alfa(:,:,inb);
    end
end

for isym = 1:Nsym
    tic;
    symbol = sprintf('symbol_%i', isym);
    for inb = [1 2 9 24 31 32 33 40 55 62 63]%1:17%Nb %[24 31 32 33 40]
        beam = sprintf('beam_%i', inb);
        for iel=1:Nt
            AFb_tay_cod_corr(:,:,iel)=Aptayl(iel).*(exp(1i*p(iel)*Tx.(symbol).(beam).beta_corr).*ones(1,length(teta))).*(exp(1i*p(iel)*(Tx.(symbol).(beam).alfa+Tx.(symbol).(beam).beta)));% codedT   
        end
        Tx.(symbol).(beam).AFb = sum(AFb_tay_cod_corr,3);
    end
    toc;
end
for ib = 1:Nb
    iteta_c (ib) = find(teta_rnd==teta_c_rnd(ib));
end

%% Initialization
k = zeros(Nsym,Nsamp);              % initialize discrete time samples
st_sc = zeros(Ncb,Nsamp,Nsym);      % initialize discrete frequency samples
st = zeros(Nsamp,Nsym,Nb);          % initialize time domain signal
st2 = zeros(Ncb,Nsym,Nb);           % initialize multicarrier signal 
st_seri = zeros(Nb,Nsym*Nsamp);     % initialize time domain signal in series
st_ref = zeros(Nb,2*Nsym*Nsamp);    % initialize the reference time domain signal 

%% Target Parameter
tta0 = [teta_c_rnd(32) teta_c_rnd(32)];    % Target direction
R = [100000 100060];                  % Target Range in meters.
v = [0 0];                      % Target velocity in m/s
rcs = [1 1];
Ntgt = length(tta0);
Afb_max = max(max(Tx.symbol_1.beam_32.AFb));
sr_seri = zeros(Ntgt,2*Nsym*Nsamp); % initialize the receive time domain signal
sr_seri2 = zeros(Ntgt,2*Nsym*Nsamp);
%% Generating time domain signal reference
for itgt =1:Ntgt
    tau0(itgt) = 2*R(itgt)/c;                       % round trip delay of the target
    ntau(itgt) = round(tau0(itgt)/T0*Nsamp);        % number of bits correspond to tau0
    tta0_rnd(itgt) = round(tta0(itgt)*10)/10;
    itta0(itgt) = find(teta_rnd==tta0_rnd(itgt));
    b0(itgt) =  find(teta_c_rnd == tta0(itgt));         % b0 define the beam index based on the assumed target direction
    beam0 = sprintf('beam_%i', itgt);                 % specify the beam_b reference
    for inb = [1 2 9 24 31 32 33 40 55 62 63]           % this supposed to be from beam 1 to Nb (setting to this specific beams just to speed up the simulation time
        beam = sprintf('beam_%i', inb);                 % specify the beam_b reference
        for isym = 1:Nsym
            k(isym,:) = (isym-1)*Nsamp:isym*Nsamp-1;    % specify the k time samples
            symbol = sprintf('symbol_%i', isym);        % specify the symbol index
            % transmit time domain signal from beam b0
%             st_sc(:,:,isym) =  1/Nsamp.*Tx.(symbol).(beam).AFb(:,itta0(itgt)).*Tx.(symbol).(beam).code.'.*(exp(1i*2*pi.*Tx.(symbol).(beam).f_band*(k(isym,:).*T0/Nsamp)));
            st_sc(:,:,isym) =  1/Nsamp.*Tx.(symbol).(beam).AFb(:,iteta_c(inb)).*Tx.(symbol).(beam).code.'.*(exp(1i*2*pi.*Tx.(symbol).(beam).f_band*(k(isym,:).*T0/Nsamp)));
            % transmit multicarrier time domain signal from beam b0      
            st (:,isym,inb)= sum(st_sc(:,:,isym),1).';
            % checker 1: fft of transmit multicarrier time domain signal from beam b0      
            st2 (:,isym,inb) = Nsamp.*(exp(-1i*2*pi.*Tx.(symbol).(beam).f_band*(k(isym,:).*T0/Nsamp)))*st (:,isym,inb);
        end
        % paralel to serial of the multiple symbol series
        st_seri(inb,:)= reshape(st(:,:,inb),1,[]);
        % zero padding at the tail series twice of the length for receiver receiving time period
        st_ref(inb,:) = [st_seri(inb,:) zeros(1,Nsym*Nsamp)];
        % when there is a target at a distance R, the round trip delay
        % equivalent to number of bits of ntau
        if inb == b0(itgt)
            % shifting the signal ntau bits
            %% IMPORTANT
            sr_seri(itgt,:) = 1/R(itgt)^4.*circshift(st_ref(inb,:),ntau(itgt));
%             sr_seri(itgt,:) = circshift(st_ref(inb,:),ntau(itgt));
            % undo shifting sr_seri2(itgt,:) must be the same as st_ref(inb,:) 
            sr_seri2(itgt,:) = circshift(sr_seri(itgt,:),-ntau(itgt));
            % serial to paralel with column index = symbol index 
            sr(:,:,itgt) = reshape(sr_seri2(itgt,:),Nsamp,[]);      
        end
    end
end
A_kor = max(xcorr(st_ref(32,:)));
%% Plotting The beampattern
figure
kara = {'y','--r',':b'};
fscp = [1 128 256];
for itgt = 1:Ntgt
    beam0 = sprintf('beam_%i', b0(itgt));         % specify the beam_b0
    for ifff = 1:length(fscp)
        AFb_plot(ifff,:,itgt)= cos(teta).*Tx.symbol_1.(beam0).AFb(fscp(ifff),:);
        figure(itgt);
    %     plot(teta*180/pi,abs(cos(teta).*Tx.symbol_1.(beam0).AFb(fscp(ifff),:)),kara{ifff},'LineWidth',2);hold on;
    %     plot(teta*180/pi,20*log10(abs(cos(teta).*Tx.symbol_1.(beam0).AFb(fscp(ifff),:))),kara{ifff},'LineWidth',2);hold on;
        plot(teta*180/pi,20*log10(abs(AFb_plot(ifff,:,itgt))),kara{ifff},'LineWidth',2);hold on;
        xlabel('\theta (deg) ');
        ylabel('Magnitude (dB)');
    %     title(['Array Factor of coded CC-OFDM beam b = ',num2str(b0-32),' (\theta_b',' = ',num2str(round(teta_c(b0)*180/pi*10)/10),char(176),')']);
        title(['Power pattern of beam b = ',num2str(b0(itgt)-32),' (\theta_b',' = ',num2str(round(teta_c(b0(itgt))*180/pi*10)/10),char(176),')']);
%         subtitle('with phase correction')
        legend(['f_1_,_',num2str(b0(itgt)-32),'    ';'f_1_2_8_,_',num2str(b0(itgt)-32);'f_2_5_6_,_',num2str(b0(itgt)-32)]);
        axis([180/pi*teta_c(b0(itgt))-15 180/pi*teta_c(b0(itgt))+15 -20 38])
        grid on;
    end
end
% hold on;plot(-0.63,33.2,'bo');
% hold on;plot(0.63,33.2,'ro');
% legend(['${f}_{1  }$';'${f}_{128}$';'${f}_{256}$';'$3dB^{-  }$';'$3dB^{+  }$'],'interpreter','latex');
% for multiple target this is where the echo signals of the targets are mixed together 
sr_seri_Ntgt = sum(sr_seri,1);
% this is the discrete time samples at the receiver
k_seri = (0:2*Nsym*Nsamp-1)*T0/Nsamp;
format long
aa = zeros(Nb,2*2*Nsym*Nsamp-1);
cor = zeros(Nb,2*Nsym*Nsamp);
sr_seri_Ntgt_ret = zeros(Nb,2*Nsym*Nsamp);
sr_par = zeros(Nsamp,2*Nsym,Nb);
% est_b = zeros(Nb,Nsym);
est_range = zeros(1,Nb);
est_ntau = zeros(1,Nb);
kor = zeros(Nsym,2*2*Nsym*Nsamp-1,Nb);
for inb = [1 2 9 24 31 32 33 40 55 62 63]%1:17%Nb
    beam = sprintf('beam_%i', inb);
    for isym = 1:Nsym
        kor1(isym,:,inb) = xcorr(sr_seri_Ntgt,st(:,isym,inb));
        if isym == 1
%             kor(isym,:,inb) = xcorr(sr_seri_Ntgt,st(:,isym,inb));
            kor_acc(isym,:,inb) = kor1(isym,:,inb);
            max_kor (inb)= max(kor_acc(isym,:,inb));
        else
%             kor1(isym,:,inb) = 
            kor(isym,:,inb) = circshift(kor1(isym,:,inb),-(isym-1)*Nsamp);           
            kor_acc(isym,:,inb) = 1/max_kor(inb).*(kor_acc(isym-1,:,inb).*kor(isym,:,inb));
        end
    end
    [a1(inb),p1(inb)] = max(kor_acc(Nsym,:,inb));
    kor_fin (inb,:) = kor_acc(Nsym,:,inb);
end
Rmax = c/2/delf*Nsym;
thr = 1/Rmax^4*Afb_max^2;
thr2 = 1./Rmax.^4.'.*Nsym;
est_b = a1>thr2*8;            % 2 targets at co-channel beams and same beams
% est_b = abs(a1./max(a1))>0.75;  % 2 targets at neighboring beams
kor_dis = est_b.'.*kor_fin;

figure(3);
subplot(211);
plot(k_seri*1e3,real(st_ref(b0(itgt),:)));
title(['Transmitted symbol of beam ',num2str(b0-32)]);
xlabel('time (ms)');ylabel('real magnitude');
subplot(212);
plot(k_seri*1e3,180/pi*angle(st_ref(b0(itgt),:)));
% subtitle('Transmitted symbol of beam 0');
xlabel('time (ms)');ylabel('Phase (deg)');

figure(4)
subplot(211);
plot(k_seri*1e3,real(st_ref(b0(itgt),:)));
title(['Transmitted symbol of beam ',num2str(b0(itgt)-32)]);
xlabel('time (ms)');ylabel('real magnitude');
subplot(212);
plot(k_seri*1e3,real(sr_seri_Ntgt));
title(['Echo signal of beam ',num2str(b0(itgt)-32)]);
xlabel('time (ms)');ylabel('real magnitude');

figure(5)
subplot(211);
plot(k_seri*1e3,real(sr_seri_Ntgt));
subtitle('Delayed signal of beam 0');
xlabel('time (ms)');ylabel('Real Magnitude');
subplot(212);
plot(k_seri*1e3,180/pi*angle(sr_seri_Ntgt));
% subtitle('Delayed signal of beam 0');
xlabel('time (ms)');ylabel('Phase (deg)');

figure(6)
subplot(311);
plot(k_seri*1e3,abs(sr_seri2(1,:)));xlabel('time (ms)');ylabel('Magnitude');
title('test-point 1')
subplot(312);
plot(k_seri*1e3,abs(sr_seri2(2,:)));xlabel('time (ms)');ylabel('Magnitude');
subplot(313);
plot(k_seri*1e3,abs(sum(sr_seri2,1)));xlabel('time (ms)');ylabel('Magnitude');

figure(7);
subplot(311);
plot(k_seri*1e3,abs(sr_seri(1,:)));xlabel('time (ms)');ylabel('Magnitude');
subtitle('signal echo from target 1');
% axis([min(tau0)*1e3-0.005 max(tau0)*1e3+0.005 -inf 1.1*abs(a1(b0(1)))]);
subplot(312);
plot(k_seri*1e3,abs(sr_seri(2,:)));xlabel('time (ms)');ylabel('Magnitude');
subtitle('signal echo from target 2');
% axis([min(tau0)*1e3-0.005 max(tau0)*1e3+0.005 -inf 1.1*abs(a1(b0(1)))]);
subplot(313);
plot(k_seri*1e3,abs(sr_seri_Ntgt));xlabel('time (ms)');ylabel('Magnitude');
subtitle('signal echo from all the targets');
% axis([min(tau0)*1e3-0.005 max(tau0)*1e3+0.005 -inf 1.1*abs(a1(b0(1)))]);

kkor = (-2*Nsym*Nsamp+1:2*Nsym*Nsamp-1)*T0/Nsamp*1e3;

% figure(8);
% subplot(311)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1)-1)));
% subtitle(['correlator output of beam ',num2str(b0(1)-1)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(312)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1))));
% subtitle(['correlator output of beam ',num2str(b0(1))]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(313)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1)+1)));
% subtitle(['correlator output of beam ',num2str(b0(1)+1)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% 
% figure(9);
% subplot(311)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1)-8)));
% subtitle(['correlator output of beam ',num2str(b0(1)-8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(312)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1))));
% subtitle(['correlator output of beam ',num2str(b0(1))]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(313)
% plot(c/2*1e-3.*kkor,abs(kor_acc(Nsym,:,b0(1)+8)));
% subtitle(['correlator output of beam ',num2str(b0(1)+8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% 
% 
% figure(10);
% subplot(311)
% plot(kkor,abs(kor_acc(Nsym,:,b0(2)-1)));
% subtitle(['correlator output of beam ',num2str(b0(2)-1)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(312)
% plot(kkor,abs(kor_acc(Nsym,:,b0(2))));
% subtitle(['correlator output of beam ',num2str(b0(2))]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(313)
% plot(kkor,abs(kor_acc(Nsym,:,b0(2)+1)));
% subtitle(['correlator output of beam ',num2str(b0(2)+1)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);

% figure(9);
% subplot(311)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1)-8)));
% subtitle(['correlator output of beam ',num2str(b0(1)-8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(312)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1))));
% subtitle(['correlator output of beam ',num2str(b0(1))]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(313)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1)+8)));
% subtitle(['correlator output of beam ',num2str(b0(1)+8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% 
% figure(10);
% subplot(311)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1)-8)));
% subtitle(['correlator output of beam ',num2str(b0(1)-8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(312)
% plot(kkor,abs(kor_acc(Nsym,:,b0(1))));
% subtitle(['correlator output of beam ',num2str(b0(1))]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(313)
% plot(c/2*1e-3.*kkor,abs(kor_acc(Nsym,:,b0(1)+8)));
% subtitle(['correlator output of beam ',num2str(b0(1)+8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);

% figure(11);
% subplot(311)
% plot(kkor,real(a1(b0(2))).*abs(kor_acc(Nsym,:,b0(2)-8)));
% subtitle(['correlator output of beam ',num2str(b0(2)-8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(312)
% plot(kkor,real(a1(b0(2))).*abs(kor_acc(Nsym,:,b0(2))));
% subtitle(['correlator output of beam ',num2str(b0(2))]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(313)
% plot(c/2*1e-3.*kkor,real(a1(b0(2))).*abs(kor_acc(Nsym,:,b0(2)+8)));
% subtitle(['correlator output of beam ',num2str(b0(2)+8)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% 
% figure(12);
% subplot(311)
% plot(kkor,real(a1(b0(1)-1)).*abs(kor_acc(Nsym,:,b0(1)-1)));
% subtitle(['correlator output of beam ',num2str(b0(1)-1)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(312)
% plot(kkor,real(a1(b0(1))).*abs(kor_acc(Nsym,:,b0(1))));
% subtitle(['correlator output of beam ',num2str(b0(1))]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);
% subplot(313)
% plot(c/2*1e-3.*kkor,real(a1(b0(1)+1)).*abs(kor_acc(Nsym,:,b0(1)+1)));
% subtitle(['correlator output of beam ',num2str(b0(1)+1)]);
% axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3,-inf abs(max(max_kor))+abs(max(max_kor))*0.1]);

figure(12)
subplot(311)
plot(kkor,10*log10(abs(kor_dis(b0(1)-1,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)-1-32)]);
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
subplot(312)
plot(kkor,10*log10(abs(kor_dis(b0(1),:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)-32)]);
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
subplot(313)
plot(kkor,10*log10(abs(kor_dis(b0(1)+1,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)+1-32)]);
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');

figure(13);
subplot(311)
plot(kkor,10*log10(abs(kor_dis(b0(1)-8,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)-8-32)]);
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
subplot(312)
plot(kkor,10*log10(abs(kor_dis(b0(1),:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)-32)]);
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
subplot(313)
plot(kkor,10*log10(abs(kor_dis(b0(1)+8,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)+8-32)]);
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');

figure(14);
subplot(311)
plot(kkor,10*log10(abs(kor_dis(b0(2)-1,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)-1-32)]);
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
subplot(312)
plot(kkor,10*log10(abs(kor_dis(b0(2),:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(313)
plot(kkor,10*log10(abs(kor_dis(b0(2)+1,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)+1-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);

figure(15);
subplot(311)
plot(kkor,10*log10(abs(kor_dis(b0(2)-8,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)-8-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(312)
semilogy(kkor,10*log10(abs(kor_dis(b0(2),:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(313)
plot(c/2*1e-3.*kkor,10*log10(abs(kor_dis(b0(2)+8,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)+8-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);

figure(16)
subplot(211);stem(a1);ylabel('Magnitude');xlabel('beam index - b');
subplot(212);stem(est_b);ylabel('Magnitude');xlabel('beam index - b');

figure(17);
subplot(411)
plot(kkor,10*log10(abs(kor_dis(b0(1)-1,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)-1-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(412)
plot(kkor,10*log10(abs(kor_dis(b0(1),:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(413)
plot(kkor,10*log10(abs(kor_dis(b0(1)+1,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)+1-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(414)
plot(kkor,10*log10(abs(kor_dis(b0(1)+2,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(2)+2-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);

figure(18);
plot(kkor,10*log10(abs(kor_dis(b0(2),:))),'LineWidth',1);grid on;
%title(['correlator output of beam ',num2str(b0(2))]);
axis([(min(tau0)-4e-6)*1e3 (max(tau0)+4e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau (ms)');
hAx(1)=gca;
hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
hold(hAx(2),'on')
plot(hAx(2),kkor*1e-6*c/2,10*log10(abs(kor_dis(b0(2),:))));
axis([(min(tau0)-4e-6)*c/2*1e-3 (max(tau0)+4e-6)*c/2*1e-3 -400 -150]);
xlabel('Range (km)');

figure(19);
plot(kkor,10*log10(abs(kor_dis(b0(1),:))),'LineWidth',1);grid on;
%title(['correlator output of beam ',num2str(b0(2))]);
axis([(min(tau0)-4e-6)*1e3 (max(tau0)+4e-6)*1e3 -400 -150]);
xlabel('time-lag, \tau (ms)');
hAx(1)=gca;
hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
hold(hAx(2),'on')
plot(hAx(2),kkor*1e-6*c/2,10*log10(abs(kor_dis(b0(1),:))));
axis([(min(tau0)-4e-6)*c/2*1e-3 (max(tau0)+4e-6)*c/2*1e-3 -400 -150]);
xlabel('Range (km)');

figure(20);
subplot(411)
plot(kkor,10*log10(abs(kor_dis(b0(1)-8,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)-8-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(412)
plot(kkor,10*log10(abs(kor_dis(b0(1),:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(413)
plot(kkor,10*log10(abs(kor_dis(b0(1)+8,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)+8-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);
subplot(414)
plot(kkor,10*log10(abs(kor_dis(b0(1)+16,:))),'LineWidth',1);grid on;
subtitle(['correlator output of beam ',num2str(b0(1)+16-32)]);
xlabel('time-lag, \tau ( ms )');ylabel('Magnitude (dB)');
axis([(min(tau0)-12e-6)*1e3 (max(tau0)+12e-6)*1e3 -400 -150]);