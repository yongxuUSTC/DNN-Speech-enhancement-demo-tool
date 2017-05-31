

clc;
 clear all;
 if ~isdeployed
addpath('help');
 end
 
 system('perl help\\step1_wav2raw.pl');
 system('perl help\\step2_get_scp_raw.pl');
system('perl help\\step3_get_lsp.pl');
system('perl help\\step4_get_scp_lsp.pl');
 
%%%%compute trainset global mean and variance
expansion_frames=7;
lsp_dim=257;
mfcc_dim=93;
all_dim=lsp_dim+mfcc_dim;
output_frames_num=lsp_dim+mfcc_dim;
% load('config\timit_aurora4_115NT_7SNRs_each190_80utts_noisy_lsp_be_random_linux_global_mv.mat');
% load('config\timit_115NT_7SNRs_each190utts_noisy_lsp_be_random_linux_global_mv.mat');
load('config\timit_115NT_7SNRs_each190utts_noisy_lspANDmfcc93b40_be_random_linux_global_mv.mat');
% load('config\timit_115NT_7SNRs_each190utts_noisy_lspANDmfcc_be_random_linux_global_mv.mat');
% ratios=ave_array_ratios;
m_global_mean=mean(global_mean(1:lsp_dim));
m_global_var=mean(global_var(1:lsp_dim));

for epoch=50:1:50
% wts_name=sprintf('config\\se_weights%d',epoch);
% weights=load(wts_name);
% weights=load('config\\se_weights50(TIMIT-16k-115NT-80H,ReLU-F6NAT-hid2500-bMFCC-DropoutV0.1H0.1-2out0.8S0.2N-energyN,epoch50,err46.15).mat');
weights=load('config\\se_weights50(TIMIT-16k-115NT-80H,ReLU-F6NAT-hid2500-bMFCC93-Mel40-DropoutV0.1H0.1-2out0.8S0.2N-energyN,epoch50,err46.51).mat');
% weights=load('config\\se_weights50(TIMIT-16k-115NT-80H,ReLU-F6NAT-DropoutV0.1H0.1-2out0.9S0.1N-energy300N,epoch50,err48.24).mat');
% weights=load('config\\se_weights50(TIMIT-16k-115NT-80H,ReLU-F6NAT-DropoutV0.1H0.1-2out0.8S0.2N-energy300N,epoch50,err48.51).mat');

noisy_lsp_list='scp\yqzu_lsp_16k_noisy.scp';
% noisy_lsp_list='scp\test.scp';
flsp=fopen(noisy_lsp_list);
tline=fgetl(flsp);
line_num=0;
while(tline~=-1)
    line_num=line_num+1;
% %     filename=regexp(tline,'(.+)\\(F01-.+)\.lsp','tokens');
%     filename=regexp(tline,'(.+)\\(.+16k.+)\.lsp','tokens');
    filename=regexp(tline,'(.+v_lsp)\\(.+)\.lsp','tokens');
    path=filename{1}{1};
    fname=filename{1}{2};
%     
%     noise_type=regexp(path,'.+(N[\d]+_SNR[\d\-]+)','tokens');
%     noise_type=noise_type{1}{1};
%    
    

%  cmd=[];
%  fid = fopen('test.list', 'wt');
%  fprintf(fid, '%s\n', tline);
%  fclose(fid);
% cmd=sprintf('MergeFeatures.exe 5 test.list input_lsp.txt');
% system(cmd);

    [htkdata,nframes,sampPeriod,sampSize,paramKind]=readhtk_new(tline,'le');
    htkdata=htkdata';
    htkdata_noisy=htkdata;
    
%     cmd=sprintf('copy %s raw_noisy.lsp',tline);
%     system(cmd);

    old_tline=tline;
%     raw_noisy_path=strcat('\\172.16.45.18\yongxu_F\TIMIT_16kHz_SEDNN_data_100h\NoisyData_en_cmu\',noise_type);
%     raw_noisy_tline=strrep(tline,path,raw_noisy_path);
%     tline=old_tline;
%     raw_noisy_tline=strrep(raw_noisy_tline,'.lsp','.raw'); 
raw_noisy_tline=strrep(tline,'.lsp','.raw'); 
    tline=old_tline;
%                    cmd=sprintf('copy %s raw_noisy.raw',raw_noisy_tline);
%     system(cmd);
        cmd=sprintf('help\\RAW2WAV.exe 1 16000 %s tmp\\noisy.wav',raw_noisy_tline);
    system(cmd);

%         raw_clean_tline=strrep(tline,path,'\\172.16.45.18\yongxu_F\TIMIT_16kHz_SEDNN_data_100h\en_cmu_test50');
%     tline=old_tline;
%     raw_clean_tline=strrep(raw_clean_tline,'.lsp','.raw'); 
%         cmd=sprintf('Raw2Wav.exe 1 16000 %s wav_clean.wav',raw_clean_tline);
%     system(cmd);
%     wav_clean_tline='wav_clean.wav';
    
    system('help\\HCopy.exe -C help\\config_mel40.mfcc tmp\\noisy.wav tmp\\noisy.mfcc');
        [htkdata_mfcc,nframes_mfcc,sampPeriod_mfcc,sampSize_mfcc,paramKind_mfcc]=readhtk_new('tmp\\noisy.mfcc','be');
    htkdata_mfcc=htkdata_mfcc';
    
    htkdata_in=[htkdata_noisy htkdata_mfcc];
    frame_expand(htkdata_in,expansion_frames); %%%写到‘input_lsp.txt’文件中了

%     out_tline=strrep(tline,'.lsp','_DNNenh_mfcc.wav'); 
%     old_tline=tline;
    
   out_noisy_tline=strrep(tline,'.lsp','.wav'); 
    old_tline=tline;
%         cmd=sprintf('d:\\bin\\RAW2WAV.exe 1 16000 %s %s',raw_noisy_tline,out_noisy_tline);
%     system(cmd);
    
%     system('Enhancement_logmmse_new_plus_1.exe raw_noisy.wav
%     raw_noisy_lmmse.wav 2');


 sources=load('tmp\\input_lsp.txt');
    [hang,lie] = size(sources);
    %%%estimated noise spectrum
    noise=zeros(1,output_frames_num);
    for t=1:1:6
    noise=noise+htkdata_in(t,:);
    end
%     for t=hang-5:1:hang
%     noise=noise+htkdata(t,:);
%     end
    noise=noise./6;
    for t=1:hang
        for d=(expansion_frames*output_frames_num+1):1:(expansion_frames*output_frames_num+output_frames_num)
        sources(t,d)=noise(1,d-expansion_frames*output_frames_num);
        end
    end
    
    
%     for i=1:expansion_frames*output_frames_num
%     sources(:,i) = ((sources(:,i)-global_mean(i))-2)/global_var(i);
%     end;
%     for i=(expansion_frames*output_frames_num+1):1:(expansion_frames*output_frames_num+257)
%     sources(:,i) = (sources(:,i)-global_mean(i-expansion_frames*output_frames_num)-2)/global_var(i-expansion_frames*output_frames_num);
%     end;

totnum=size(sources,1);
% fprintf(1, 'Size of the test dataset= %5d \n', totnum);

% [batchsize,lie] = size(sources);
batchsize=1;
[hang,lie] = size(sources);
numbatches=floor(totnum/batchsize);

numdims  =  size(sources,2);
testbatchsources = zeros(batchsize, numdims, numbatches);

for b=1:numbatches
  testbatchsources(:,:,b) = sources((1+(b-1)*batchsize:b*batchsize), :);
end;

[N testnumdims testnumbatches]=size(testbatchsources);

data=zeros(testnumdims);
dataout=zeros(testnumbatches,lsp_dim*2+mfcc_dim);

%%%the estimation of mean bias and discriminative normalization
data=zeros(nframes,expansion_frames*output_frames_num+output_frames_num);
gmbias=zeros(nframes,1);
for batch = 1:testnumbatches
  data(batch,:) = [testbatchsources(:,:,batch)];


%    mdata=data(((expansion_frames-1)/2)*output_frames_num+1:1:((expansion_frames-1)/2+1)*output_frames_num);
% %    mdata=data;
%    if (mean(mdata)<mean(global_mean))
%        gmbias=0;
%    else 
%     gmbias  =  mean(mdata)-mean(global_mean);
%    end
   
%    mdata=data(batch,((expansion_frames-1)/2)*output_frames_num+1:1:((expansion_frames-1)/2+1)*output_frames_num);
%    m_mdata=mean(mdata);
% %    mdata=data(batch,:);
% if((m_mdata<(m_global_mean+m_global_var/2))&(m_mdata>(m_global_mean-m_global_var/2)))
%     gmbias(batch)=0;
% elseif (m_mdata>=(m_global_mean+m_global_var/2))
%     gmbias(batch)  =  m_mdata-(m_global_mean+m_global_var/2);
% else
%    gmbias(batch)  = m_mdata-(m_global_mean-m_global_var/2);
% end
% %    aa(batch)=gmbias;

noisy_dim=output_frames_num-mfcc_dim;
  mdata=data(batch,((expansion_frames-1)/2)*output_frames_num+1:1:((expansion_frames-1)/2)*output_frames_num+noisy_dim);
   m_mdata=mean(mdata);
if((m_mdata<=(m_global_mean+m_global_var/2))&(m_mdata>=(m_global_mean-m_global_var/2)))
    gmbias(batch)=0;
elseif((m_mdata<=(m_global_mean+m_global_var))&(m_mdata>(m_global_mean+m_global_var/2)))
    gmbias(batch)=m_mdata-(m_global_mean+m_global_var/2); 
elseif (m_mdata>(m_global_mean+m_global_var))
    gmbias(batch)  =  m_mdata-(m_global_mean+m_global_var);
 elseif((m_mdata<(m_global_mean-m_global_var/2))&(m_mdata>=(m_global_mean-m_global_var)))
    gmbias(batch)=m_mdata-(m_global_mean-m_global_var/2);  
else
   gmbias(batch)  = m_mdata-(m_global_mean-m_global_var);
end
%    aa(batch)=gmbias;

    for i=1:expansion_frames*output_frames_num
        if((i>lsp_dim & i<=all_dim) || (i>(2*lsp_dim+1*mfcc_dim) & i<=2*all_dim) ||  (i>(3*lsp_dim+2*mfcc_dim) & i<=3*all_dim) || (i>(4*lsp_dim+3*mfcc_dim) & i<=4*all_dim) || (i>(5*lsp_dim+4*mfcc_dim) & i<=5*all_dim) || (i>(6*lsp_dim+5*mfcc_dim) & i<=6*all_dim) || (i>(7*lsp_dim+6*mfcc_dim) & i<=7*all_dim))
        data(batch,i) = ((data(batch,i)-global_mean(i)))/global_var(i);
        else
        data(batch,i) = ((data(batch,i)-global_mean(i))-gmbias(batch))/global_var(i);    
        end
    end;
    for i=(expansion_frames*output_frames_num+1):1:(expansion_frames*output_frames_num+lsp_dim)
    data(batch,i) = (data(batch,i)-global_mean(i-expansion_frames*output_frames_num)-gmbias(batch))/global_var(i-expansion_frames*output_frames_num);
    end
       for i=(expansion_frames*output_frames_num+lsp_dim+1):1:(expansion_frames*output_frames_num+lsp_dim+mfcc_dim)
    data(batch,i) = (data(batch,i)-global_mean(i-expansion_frames*output_frames_num))/global_var(i-expansion_frames*output_frames_num);
    end
end

data=[data ones(nframes,1)];
%%%%%DNN decoding
  %%%% dropout + ReLU
  [size1,size2]=size(weights.w1);
  new_w1=[weights.w1(1:1:size1-1,:).*0.9;weights.w1(size1,:)];
%   w1probs = 1./(1 + exp(-data*weights.w1));
  w1probs=data*new_w1;
  w1probs(w1probs<=0)=0;
  w1probs = [w1probs  ones(nframes,1)];

    [size1,size2]=size(weights.w2);
  new_w2=[weights.w2(1:1:size1-1,:).*0.9;weights.w2(size1,:)];
%  w2probs = 1./(1 + exp(-w1probs*weights.w2)); 
   w2probs=w1probs*new_w2;
   w2probs(w2probs<=0)=0;
%    clear w1probs;
  w2probs = [w2probs ones(nframes,1)];
  
    [size1,size2]=size(weights.w3);
  new_w3=[weights.w3(1:1:size1-1,:).*0.9;weights.w3(size1,:)];  
%   w3probs = 1./(1 + exp(-w2probs*weights.w3));
  w3probs=w2probs*new_w3;
  w3probs(w3probs<=0)=0;
%    clear w2probs;
  w3probs = [w3probs ones(nframes,1)];
  
    [size1,size2]=size(weights.w4);
  new_w4=[weights.w4(1:1:size1-1,:).*0.9;weights.w4(size1,:)];
  dataout = w3probs*new_w4;
  
  %%%%inverse norm
for t=1:nframes
for i=1:1:(all_dim+lsp_dim)
%     dataout(t,i) = (dataout(t,i).*global_var(i)+global_mean(i)+gmbias(t));
dataout(t,i) = (dataout(t,i).*global_var(i)+global_mean(i));
end
end

dataout_p1=dataout(:,1:1:lsp_dim);
dataout_p2=dataout(:,all_dim+1:1:(all_dim+lsp_dim));


%%%%%%%%IBMpp
[frames,dims]=size(dataout_p1);
DNNenh_f=zeros(frames,dims);
irm=zeros(frames,dims);
% alpha=zeros(frames,1);
for t=1:frames 
    for d=1:dims
irm(t,d)=sqrt(exp(dataout_p1(t,d))/(exp(dataout_p1(t,d))+exp(dataout_p2(t,d))));  
% irm(t,d)=sqrt(exp(dataout_p1(t,d))/(exp(htkdata_noisy(t,d))));  
    end
%     alpha(t)=mean(irm(t,:));
end
for t=1:frames 
%     if(mean(irm(t,:))<0.1) %%%%若这一帧，整体语音存在概率很低，说明非语音帧，直接置0，扔掉
%     DNNenh_f(t,:)=dataout_p1(t,:)*0.8;
%     continue;
%     end
    for d=1:dims
%     %irm=sqrt(exp(dataout_p1(t,d))/(exp(dataout_p1(t,d))+exp(dataout_p2(t,d))));    
%     if(irm(t,d)>alpha(t)*1.3)
%         DNNenh_f(t,d)=htkdata_noisy(t,d); 
%     elseif(irm(t,d)<0.2)
%         DNNenh_f(t,d)=dataout_p1(t,d)*0.8;
%     else
%         DNNenh_f(t,d)=dataout_p1(t,d); 
%     end

% % %%%%思路二：irm当做weina gain
    if(irm(t,d)>0.75)
% %         DNNenh_f(t,d)=htkdata_noisy(t,d);
  %%%%inverse norm
    DNNenh_f(t,d) =( (dataout(t,d)+gmbias(t))+htkdata_noisy(t,d))/2;

    elseif(irm(t,d)<0.1)
%         DNNenh_f(t,d)=dataout_p1(t,d)*0.8;      
%         DNNenh_f(t,d)=dataout(t,d).*global_var(d)+global_mean(d)+gmbias(t)/2;
          DNNenh_f(t,d)=dataout(t,d);
    else
% %         DNNenh_f(t,d)=log((irm(t,d))^2*exp(htkdata_noisy(t,d)));
         DNNenh_f(t,d)= (dataout(t,d)+gmbias(t)-1);
    end

%     else
%          DNNenh_f(t,d)= (dataout(t,d)+gmbias(t)*irm(t,d));
%     end

% DNNenh_f(t,d)= dataout_p1(t,d);
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % %%%%%%IBMpp2
% % DNNenh_f2=zeros(frames,dims);
% % irm2=zeros(frames,dims);
% % alpha2=zeros(frames,1);
% % for t=1:frames 
% %     for d=1:dims
% %    irm2(t,d)=sqrt(exp(DNNenh_f(t,d))/(exp(DNNenh_f(t,d))+exp(dataout_p2(t,d))));  
% %     end
% %     alpha2(t)=mean(irm2(t,:));
% % end
% % for t=1:frames 
% %     if(alpha2(t)<0.2) %%%%若这一帧，整体语音存在概率很低，说明非语音帧，直接置0，扔掉
% %     DNNenh_f2(t,:)=dataout_p1(t,:)*0.8;
% %     else
% %     DNNenh_f2(t,:)=DNNenh_f(t,:);    
% %     end
% % end
% % DNNenh_f=DNNenh_f2;
% 
% % %%%%%最大值修正，纯粹为了ppt好看
% % [nT,nF]=size(DNNenh_f);
% % max_DNNenh_f=max(DNNenh_f(:));
% % for t=1:nT
% %     for f=1:nF
% %     if(DNNenh_f(t,f) > max_DNNenh_f-1.5)
% %     DNNenh_f(t,f)=DNNenh_f(t,f)-1.5;
% %     end
% %     end
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writeHTK_new('tmp\\out.htk', DNNenh_f, testnumbatches, 160000, lsp_dim*4, 9, 'le');
% writeHTK_new('out.htk', dataout_p1, testnumbatches, 160000, 257*4, 9, 'le');
%cmd=sprintf('LogSpec2Wav.exe %s %s out.htk info.txt outfile_proc.wav -F RAW -fs 8',raw_clean_tline,raw_noisy_tline);
cmd=sprintf('help\\LogSpec2raw_16bit_withoutXF.exe %s %s tmp\\out.htk tmp\\info.txt tmp\\DNN_enh.raw -F RAW -fs 16',raw_noisy_tline,raw_noisy_tline);
system(cmd);

system('help\\Raw2Wav.exe 1 16000 tmp\\DNN_enh.raw tmp\\DNN_enh.wav');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmd=sprintf('copy tmp\\DNN_enh.wav %s\\%s_DNNenh.wav',path,fname);
system(cmd);




% %%%%%%% plot correpondance figure 
% diffv=zeros(frames,dims);
%     cmd=sprintf('Wav2LogSpec.exe -F RAW -fs 16 %s clean.lsp',raw_clean_tline);
%     system(cmd);       
%     [clean_lsp,nframes,sampPeriod,sampSize,paramKind]=readhtk_new('clean.lsp','le');
%     clean_lsp=clean_lsp';   
% for t=1:frames 
%     for d=1:dims     
%    diffv(t,d)= clean_lsp(t,d)-DNNenh_f(t,d) ;
%     end   
% end
% I=rot90(diffv); 
% imagesc(I);
% %%%%%



% break;
tline=fgetl(flsp);  %% goto next line
end

%  ave_pesq=mean(pesq_mat);
%  ave_pesq=roundn(ave_pesq,-4);
%  ave_segSNR=mean(segSNR_mat);
%   ave_segSNR=roundn(ave_segSNR,-4);
%  ave_lsd=mean(LSD_mat);
%   ave_lsd=roundn(ave_lsd,-4);
%    ave_STOI=mean(STOI_mat);
%   ave_STOI=roundn(ave_STOI,-4);
%  save statistics_N5_N6_6SNRs ave_pesq ave_segSNR ave_lsd ave_STOI pesq_mat segSNR_mat LSD_mat STOI_mat;
% % save statistics_N5_N6_6SNRs ave_pesq    pesq_mat  ;

end

% clear all;
% load('statistics_N5_N6_6SNRs.mat');
