clear all;


%% Define simulation setup
%processing capabilities
B = 20*10^6; %bandwidth
f =  1*10^6;
C = 20;
%%%D = 10*10^6;
S_u= 1*10^3;
debug_1 = zeros(1,1000);
debug_2 = zeros(1,1000);
%Number of Monte Carlo setups
nbrOfSetups = 4;%3 30;
aaaaaa = 0;
outputs = zeros(1,1000);
const_outputs = zeros(1,1000);
best_outputs = zeros(1,1000);
indices = zeros(1,nbrOfSetups);
%Number of channel realizations per setup
nbrOfRealizations = 1; %1000;
%Number of APs in the cell-free network
L = 10;%25;%20; %25; %20;% 10; %original: 10

%Number of UEs
K = 6;

%Number of antennas per AP
N = 1;

%Length of the coherence block
tau_c = 200;

%Number of pilots per coherence block
tau_p = 20;

%Uplink transmit power per UE (mW)
p = 100;

%power ratio
ada_1 = 1;
%beta1 = 100;
beta1 = 1; %keep it 10 for all channels
beta2 = 1;  %% beta2 = 100 for chan 1 keep beta2 10 for chans2 and 1 for chan 3 
%beta2 = 1000;
% beta3 = 10^3;
% beta4 = 10^7;
% beta5 = 10^7;
% beta3 = 10^3;
% beta4 = 10^7;
% % % % % % % % % fileidfun = fopen('BPL3_1000m_VER15_20AP_6UE_SHAD_Hmat_scenario3_10APs_Multi_2_users_FUNFINALACTUREG_complex_part2_7_28_ExtendingTesting_ver62.txt');
%%%%%%%%%%%%%%%%%fileidfun = fopen('STOP_bbPLS_vREDO_bbvMOO_PLS_STOP_VER819_NOMOvACTUAL_LARGE_BPL3_1000m_VER15_25AP_6UE_SHAD_Hmat_scenario3_10APs_Multi_2_users_FUNFINALACTUREG_complex_part2_7_28_ExtendingTesting_ver62.txt');
fileidfun = fopen('EFBACTUAL_RREDO_100meters_PRETTY_PLS_VER551_MULTIPLE_PATHS_LARGE_VER1_10AP_6_UE_50000SHAD_part2.txt');
id_count = 0;
INDEX_COUNTER = 0;
all_states = [];
all_gain = [];
%log2(1+((p*ada_1*alphasss^2*abs(H_AP).^2)/(1+p*ada_1*((1-alphasss)^2)*abs(H_AP)^2)))
for iiii = 1:6000%4000%300000%1000
    d = str2num(fgetl(fileidfun));
    
    for slen = 1:2:(L*2)-1%19
        all_gain = [all_gain abs(d(slen)+i*d(slen+1))];        
    end
    
    all_states = [all_states all_gain.'];
    all_gain = [];
    id_count = id_count + 1;
    if mod(id_count,6) ~=0
        continue
    end
% % % %     if iiii == 1 %initialize both user and AP positions
% % % %         [gainOverNoisedB,R,pilotIndexCF,pilotIndexSC,APpositions,UEpositions] = generateSetup_threeslope_rev(L,K,N,tau_p,1,p);
% % % %     else
% % % %          [gainOverNoisedB,R,pilotIndexCF,pilotIndexSC] = generateSetup_threeslope_rev_justuserpos_change22(L,K,N,tau_p,1,p,APpositions,UEpositions); 
% % % %     end
% % % %         betaVal = db2pow(gainOverNoisedB);
% % % %    [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndexCF,p);
% % % %    %Hhat_AP = reshape(Hhat_AP(:,nbrOfRealizations,:),[N*L K]);
% % % %    H_AP = reshape(H_AP(:,nbrOfRealizations,:),[N*L K]);
    for n = 1:nbrOfSetups

%        [gainOverNoisedB,R,pilotIndexCF,pilotIndexSC] = generateSetup_threeslope(L,K,N,tau_p,1,p);
%         betaVal = db2pow(gainOverNoisedB);


        %Full transmit power case

        %Generate channel realizations, channel estimates, and estimation
        %error correlation matrices for all UEs to the APs 
%         [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndexCF,p);
        %old
        %fun = @(alphasss) max(((alphasss*D*C)/f) + ((alphasss*S_u)/(B*log2(1+((p*ada_1*alphasss.^2*abs(H_AP).^2)/(1+p*ada_1*((1-alphasss).^2)*abs(H_AP).^2))))));
        %initial condition
       
%         fun = @(alphasss) ((p*ada_1*alphasss.^2*abs(H_AP).^2))/(1+p*ada_1*((1-alphasss).^2)*abs(H_AP).^2);
%         fun = @(alphasss) sum(((log2(1+((p*ada_1*alphasss.^2*abs(H_AP).^2)/(1+p*ada_1*((1-alphasss).^2)*abs(H_AP).^2))))));
       %max
        
         %added latency of sending H_AP - may need to revise rate for this?
         
         %Hhat_AP + sending fragement from AP to central AP + sending
         %alphas!!!
         %(beta3*(sum(alphasss.*D)/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'))))))) 
%            fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (beta4*(((L-1)*(length(H_AP)+length(alphasss)))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))))  +   (beta5*(((L-1)*(length(alphasss)))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'))))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));    
          %%%%%%%%%%%fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));    
% % % % % % % % % % % % % % % % % % %          fun = @(alphasss) max((beta1*(sum(((alphasss.*D*C)/f),2))).' + (beta2*(sum((alphasss*S_u),2).'./sum((B*log2(1+(p* abs(sum(ada_1*alphasss.*H_AP,2).').^2./(1+p* mean(abs(sum(ada_1*(1-alphasss).*H_AP,2).').^2))))),1))));  
% % % % % % % % % % % % % % % % % % % % % % % % % % % %          fun = @(alphasss) max((beta1*(sum(((alphasss.*D*C)/f),2))).' + (beta2*(sum((alphasss*S_u),2).'./sum((B*log2(1+(p* abs(sum(ada_1*alphasss.*H_AP,2).').^2./(1+p* mean(abs(sum(ada_1*(1-alphasss).*H_AP,2).').^2))))),1))));  

% % % % % % % % % % % % % % % % % % %              fun = @(alphasss) max((beta1*(sum(((alphasss.*S_u*C)/f),2))).' + (beta2*(sum((alphasss*S_u),2).'./sum((B*log2(1+(p* abs(sum(ada_1*alphasss.*all_states,2).').^2./(1+p* mean(abs(sum(ada_1*(1-alphasss).*all_states,2).').^2))))),1))));  
             
             
             
             fun = @(alphasss) max((beta1*(sum(((alphasss.*S_u*C)/f),2)))) + max(max(beta2*((alphasss*S_u)./(B*log2(1+(sum(p*ada_1*alphasss.^2.*all_states.^2,2)./(1+sum(p*ada_1*((1-alphasss).^2).*all_states.^2,2)))))),[],1));  
%          fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (beta3*((alphasss.*(L-1)*D)/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))))  + (beta4*(((L-1)*(length(Hhat_AP)+length(alphasss)))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.'))))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))));
         
         








         %Replace H_AP with Hhat_AP
%                            fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (((L-1)*length(Hhat_AP))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))));
         
         %Replace Length(H_AP) with product of dimensions of H_AP or * 2?
%                   fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (((L-1)*length(H_AP))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));
         
          %added latency of AP to AP sending decoded fragments
%          fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + ((alphasss.*(L-1)*D)/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));
         
         
% % % % %         fun = @(alphasss) max(((alphasss.*D*C)/f) + ((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'))))));
% % % % %         fun = @(alphasss) max(((alphasss.*D*C)/f) + ((alphasss*S_u)./(B*log2(1+SINR)))); 
        if n == 1
        alphasss_0= rand(L,K);
           for jajj = 1:K %normalize each column so alphas sum up to 1
             alphasss_0(:,jajj) = alphasss_0(:,jajj)/sum(alphasss_0(:,jajj));
           end
        %%%outputs(n) = outputs(n) + fun(alphasss_0);

        else
        alphasss_0 = alphasss;
        end

%         lb = zeros(1,L);%0;
%         ub = ones(1,L);%;1;
%         Aeq = ones(1,L);
%         beq = 1;
        A = [];
        b = [];
        
% % % % % % % % % % %        tde = max(max(beta2*((alphasss_0*S_u)./(B*log2(1+(sum(p*ada_1*alphasss_0.^2.*all_states.^2,1)./(1+sum(p*ada_1*((1-alphasss_0).^2).*all_states.^2,1)))))),[],1));

        lb = zeros(1,L*K);%0;
        ub = ones(1,L*K);%;1;
        Aeq = [[ones(1,L) zeros(1,K*L-L)];[zeros(1,L) ones(1,L) zeros(1,K*L-2*L)]; [zeros(1,2*L) ones(1,L) zeros(1,K*L-3*L)]; [zeros(1,3*L) ones(1,L) zeros(1,K*L-4*L)];  [zeros(1,4*L) ones(1,L) zeros(1,K*L-5*L)];[zeros(1,K*L-L) ones(1,L)]];
%  3         Aeq = [[ones(1,L) zeros(1,K*L-L)];[zeros(1,L) ones(1,L) zeros(1,L)]; [zeros(1,K*L-L) ones(1,L)]];
%  2       Aeq = [[ones(1,L) zeros(1,K*L-L)];[zeros(1,K*L-L) ones(1,L)]];    
        beq = ones(K,1);

%         if n > 1
%             outputs(n) = outputs(n) + fun(alphasss);
%         end
% % %        alphasss = fmincon(fun,alphasss_0,A,b,Aeq,beq,lb,ub);
% % %        %'MaxFunEvals',4500,'TolCon', 1e-10'  6500


      options = optimoptions(@fmincon,'MaxFunEvals',4500,'StepTolerance',1e-10,'TolCon',1e-10);%'MaxIter',2000,'TolCon',1e-10); %2000 TolCon 1e-10 TolX
       alphasss = fmincon(fun,alphasss_0,A,b,Aeq,beq,lb,ub,[],options);
 

        d= 1;
        
       


    end
     outputs(iiii/K) =  fun(alphasss);
     
     %uniform alphas
     %const_alphasss = 1/L*ones(1,L);
     const_alphasss = 1/L*ones(L,K);
     const_outputs(iiii/K) = max((beta1*(sum(((const_alphasss.*S_u*C)/f),2)))) + max(max(beta2*((const_alphasss*S_u)./(B*log2(1+(sum(p*ada_1*const_alphasss.^2.*all_states.^2,2)./(1+sum(p*ada_1*((1-const_alphasss).^2).*all_states.^2,2)))))),[],1));
% % % % % % % % % %      const_outputs(iiii) = max((beta1*(sum(((const_alphasss.*S_u*C)/f),2))).' + (beta2*(sum((const_alphasss*S_u),2).'./sum((B*log2(1+(p* abs(sum(ada_1*const_alphasss.*all_states,2).').^2./(1+p* mean(abs(sum(ada_1*(1-const_alphasss).*all_states,2).').^2))))),1))));
% % % % % % % % % % %      const_outputs(iiii) = max((beta1*(sum(((const_alphasss.*D*C)/f),2))).' + (beta2*(sum((const_alphasss*S_u),2).'./sum((B*log2(1+(p* abs(sum(ada_1*const_alphasss.*H_AP,2).').^2./(1+p* mean(abs(sum(ada_1*(1-const_alphasss).*H_AP,2).').^2))))),1))));
     %const_outputs(iiii)= max((beta1*((const_alphasss.*D*C)/f)) + (beta2*((const_alphasss*S_u)./(B*log2(1+((p*ada_1*const_alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-const_alphasss).^2).*abs(H_AP).^2.')))))));
     
     %strategy of choosing AP with largest channel gain (best channel)
     %best_alphasss = zeros(1,L);
     best_alphasss = zeros(L,K);
     [ti,tp] = max(all_states);
     %[ti,tp] = max(abs(H_AP));
     best_alphass(tp(1),1) = 1;
     best_alphass(tp(2),2) = 1;
     best_alphass(tp(3),3) = 1;
     best_alphass(tp(4),4) = 1;
     best_alphass(tp(5),5) = 1;
     best_alphass(tp(6),6) = 1;
     best_outputs(iiii/K) = max((beta1*(sum(((best_alphasss.*S_u*C)/f),2)))) + max(max(beta2*((best_alphasss*S_u)./(B*log2(1+(sum(p*ada_1*best_alphasss.^2.*all_states.^2,2)./(1+sum(p*ada_1*((1-best_alphasss).^2).*all_states.^2,2)))))),[],1));
     %%%%best_outputs(iiii) = max((beta1*(sum(((best_alphasss.*S_u*C)/f),2))).' + (beta2*(sum((best_alphasss*S_u),2).'./sum((B*log2(1+(p* abs(sum(ada_1*best_alphasss.*all_states,2).').^2./(1+p* mean(abs(sum(ada_1*(1-best_alphasss).*all_states,2).').^2))))),1))));
% % % % % % % % % % %      best_outputs(iiii) = max((beta1*(sum(((best_alphasss.*D*C)/f),2))).' + (beta2*(sum((best_alphasss*S_u),2).'./sum((B*log2(1+(p* abs(sum(ada_1*best_alphasss.*H_AP,2).').^2./(1+p* mean(abs(sum(ada_1*(1-best_alphasss).*H_AP,2).').^2))))),1))));
     %best_outputs(iiii) = max((beta1*((best_alphasss.*D*C)/f)) + (beta2*((best_alphasss*S_u)./(B*log2(1+((p*ada_1*best_alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-best_alphasss).^2).*abs(H_AP).^2.')))))));
     
     all_states = [];
     id_count = 0;
     
     if iiii == 3996
         dde= 1;
     end
     
% % % %      debug_1(iiii) = max(beta1*((alphasss.*D*C)/f));
% % % %      debug_2(iiii) = max((beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));
% % % %      %%%%New Stuff put Hhat or no
% % % %         SINRssss = ((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'));
% % % %         [IIIIIII,BBBBBBB1] = sort(SINRssss,'ascend');
% % % %         [IIIIIII2,BBBBBBB2] = sort(alphasss,'ascend');
% % % %         
% % % %         if any(BBBBBBB1 ~= BBBBBBB2)
% % % %             aaaaaa = aaaaaa + 1;
% % % %         end
% % % %         %2 APs
% % % % %            fprintf(fileID,'%.5f %.5f\n',SINRssss(1),...
% % % % %             SINRssss(2));
% % % % %         
% % % % %          fprintf(fileID2,'%.5f %.5f \n',alphasss(1),...
% % % % %             alphasss(2));
% % % %         
% % % %         fprintf(fileID,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',SINRssss(1),...
% % % %             SINRssss(2),SINRssss(3),SINRssss(4),SINRssss(5),SINRssss(6),SINRssss(7),...
% % % %             SINRssss(8),SINRssss(9),SINRssss(10));
% % % %         
% % % %          fprintf(fileID2,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',alphasss(1),...
% % % %             alphasss(2),alphasss(3),alphasss(4),alphasss(5),alphasss(6),alphasss(7),...
% % % %             alphasss(8),alphasss(9),alphasss(10));
    
    
end


%for consistency with other methods
outputs(1) = [];
const_outputs(1) = [];
best_outputs(1) = [];


%CDF calculation 
%statsss = 10:10:300; %% chan 1
%%%%statsss = 10:50:1000;  %%% chan 2/3
%%%%statsss = 5:5:1000;  %%% chan 2/3
statsss = 0:0.001:0.2;
cdffff = zeros(1,length(statsss));
cdffff_const = zeros(1,length(statsss));
cdffff_best = zeros(1,length(statsss));
for j = 1:length(statsss)
    tata = length(find(outputs<statsss(j)));
    tat2 = length(find(const_outputs<statsss(j)));
    tat3 = length(find(best_outputs<statsss(j)));
    cdffff(j) = (tata/999);   %1-
    cdffff_const(j) = (tat2/999);
    cdffff_best(j) = (tat3/999);
end

plot(statsss, cdffff,'-r');
hold on
plot(statsss, cdffff_const,'-b');
plot(statsss, cdffff_best,'-g');
hold off
%title('Latency (\mu_1 = 10, \mu_2 = 1) vs. Timesteps');
xlabel('Max Latency (s)');
ylabel('CDF');
legend('CDF - Interior-Point method', 'CDF - Uniform Method', 'CDF - Best Channel Method');