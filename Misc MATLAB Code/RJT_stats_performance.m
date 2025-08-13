%%
%Stat Analysis of discrete performance and stiffness data for 10-5 study
Data=readtable('Moments and performance variables New.xlsx');
%Data = standardizeMissing(Data,-99);
col=numel(Data.Properties.VariableNames);
%%

for i=1:col
    mean_values(i)=mean(Data{:,i});
    sd_values(i)=std(Data{:,i});

end
%Check for Normality
%norm=[];
%for i=1:col
    %a=normplot(Data{:,i});
    %b=lillietest(Data{:,i});
    %norm=[norm,b];
    %ax=gcf
    %exportgraphics(ax,'LinePlot.jpg')   
%end
%Rearrange data for 2x2 Anova
%SLDJ vs SLRJT

%for i= 1:7
  % [stats,tbl,p,c,m, SSQs, DFs, MSQs, Fs, Ps]=set_anova(Data,mean_values,sd_values);
%end
%DLDJ VS DLRJT joint stiff
%for i=1:3
    %[stats,tbl,p,c,m, SSQs, DFs, MSQs, Fs, Ps]=set_anova(Data,mean_values,sd_values);
%end
%JH
%sldjinj=Data.JumpHeight_SLDJ_Inj;
%sldjnon=Data.JumpHeight_SLDJ_Non;
%sljrtinj=Data.JumpHeight_SLRJT_Inj;
%sljrtnon=Data.JumpHeight_SLRJT_Non;
%sldjasy=Data.JumpHeight_SLDJ_Asy;
%slrjtasy=Data.JumpHeight_SLRJT_Asy;
%stat_analysis(sldjinj,sldjnon,sljrtinj,sljrtnon,sldjasy,slrjtasy);
%CT
%sldjinj=Data.CT_SLDJ_Inj;
%sldjnon=Data.CT_SLDJ_Non;
%sljrtinj=Data.CT_SLRJT_Inj;
%sljrtnon=Data.CT_SLRJT_Non;
%sldjasy=Data.CT_SLDJ_Asy;
%slrjtasy=Data.CT_SLRJT_Asy;
%stat_analysis(sldjinj,sldjnon,sljrtinj,sljrtnon,sldjasy,slrjtasy);
%RSI
%sldjinj=Data.RSI_SLDJ_Inj;
%sldjnon=Data.RSI_SLDJ_Non;
%sljrtinj=Data.RSI_SLRJT_Inj;
%sljrtnon=Data.RSI_SLRJT_Non;
%sldjasy=Data.RSI_SLDJ_Asy;
%slrjtasy=Data.RSI_SLRJT_Asy;
%stat_analysis(sldjinj,sldjnon,sljrtinj,sljrtnon,sldjasy,slrjtasy);


%SLDJ Stiffness
%sldj_ank_inj= Data.Ankle_Stiff_SLDJ_Inj;
%sldj_ank_non= Data.Ankle_Stiff_SLDJ_Non;
%sldj_knee_inj= Data.Knee_Stiff_SLDJ_Inj;
%sldj_knee_non= Data.Knee_Stiff_SLDJ_Non;
%sldj_hip_inj= Data.Hip_Stiff_SLDJ_Inj;
%sldj_hip_non= Data.Hip_Stiff_SLDJ_Non;
%stat_analysis(sldj_ank_inj,sldj_ank_non,sldj_knee_inj,sldj_knee_non,sldj_hip_inj,sldj_hip_non);

%SLRJT Stiffness
%slrjt_ank_inj= Data.Ankle_Stiff_SLRJT_Inj;
%slrjt_ank_non= Data.Ankle_Stiff_SLRJT_Non;
%slrjt_knee_inj= Data.Knee_Stiff_SLRJT_Inj;
%slrjt_knee_non= Data.Knee_Stiff_SLRJT_Non;
%slrjt_hip_inj= Data.Hip_Stiff_SLRJT_Inj;
%slrjt_hip_non= Data.Hip_Stiff_SLRJT_Non;
%stat_analysis(slrjt_ank_inj,slrjt_ank_non,slrjt_knee_inj,slrjt_knee_non,slrjt_hip_inj,slrjt_hip_non);

%%
%DLDJ Stiffness
%dldj_ank_inj= Data.Ankle_Stiff_DLDJ_Inj;
%dldj_ank_non= Data.Ankle_Stiff_DLDJ_Non;
%dldj_knee_inj= Data.Knee_Stiff_DLDJ_Inj;
%dldj_knee_non= Data.Knee_Stiff_DLDJ_Non;
%dldj_hip_inj= Data.Hip_Stiff_DLDJ_Inj;
%dldj_hip_non= Data.Hip_Stiff_DLDJ_Non;
%stat_analysis(dldj_ank_inj,dldj_ank_non,dldj_knee_inj,dldj_knee_non,dldj_hip_inj,dldj_hip_non);
%%
%DLRJT Stiffness
%dlrjt_ank_inj= Data.Ankle_Stiff_RJT_Inj;
%dlrjt_ank_non= Data.Ankle_Stiff_RJT_Non;
%dlrjt_knee_inj= Data.Knee_Stiff_RJT_Inj;
%dlrjt_knee_non= Data.Knee_Stiff_RJT_Non;
%dlrjt_hip_inj= Data.Hip_Stiff_RJT_Inj;
%dlrjt_hip_non= Data.Hip_Stiff_RJT_Non;
%stat_analysis(dlrjt_ank_inj,dlrjt_ank_non,dlrjt_knee_inj,dlrjt_knee_non,dlrjt_hip_inj,dlrjt_hip_non);

%%
%DLDJ v DLRJT
%dldj_jh= Data.JumpHeight_DLDJ;
%dlrjt_jh=Data.JumpHeight_RJT;
%dldj_ct= Data.CT_DLDJ; 
%dlrjt_ct=Data.CT_RJT;
%dldj_rsi=Data.RSI_DLDJ;
%dlrjt_rsi=Data.RSI_RJT;
%stat_analysis(dldj_jh,dlrjt_jh,dldj_ct,dlrjt_ct,dldj_rsi,dlrjt_rsi);

% COM Stiffness
%dldj_stiff= Data.COM_Stiff_DLDJ;
%dlrjt_stiff=Data.COM_Stiff_RJT;
%sldj_stiff_acl= Data.COM_SLDJ_Inj; 
%slrjt_stiff_acl=Data.COM_SLRJT_Inj;
%sldj_stiff_nonacl=Data.COM_SLDJ_Non;
%slrjt_stiff_nonacl=Data.COM_SLRJT_Non;
%stat_analysis(dldj_stiff,dlrjt_stiff,sldj_stiff_acl,sldj_stiff_nonacl,slrjt_stiff_acl,slrjt_stiff_nonacl);


% DLDJ v DLRJT joint stiffness aclr
%dldj_ank_inj= Data.Ankle_Stiff_DLDJ_Inj;
%dlrjt_ank_inj= Data.Ankle_Stiff_RJT_Inj;
%dldj_knee_inj= Data.Knee_Stiff_DLDJ_Inj;
%dlrjt_knee_inj= Data.Knee_Stiff_RJT_Inj;
%dldj_hip_inj= Data.Hip_Stiff_DLDJ_Inj;
%dlrjt_hip_inj= Data.Hip_Stiff_RJT_Inj;
%stat_analysis(dldj_ank_inj,dlrjt_ank_inj,dldj_knee_inj,dlrjt_knee_inj,dldj_hip_inj,dlrjt_hip_inj);

% DLDJ v DLRJT joint stiffness nonaclr
%dldj_ank_non= Data.Ankle_Stiff_DLDJ_Non;
%dlrjt_ank_non= Data.Ankle_Stiff_RJT_Non;
%dldj_knee_non= Data.Knee_Stiff_DLDJ_Non;
%dlrjt_knee_non= Data.Knee_Stiff_RJT_Non;
%dldj_hip_non= Data.Hip_Stiff_DLDJ_Non;
%dlrjt_hip_non= Data.Hip_Stiff_RJT_Non;
%stat_analysis(dldj_ank_non,dlrjt_ank_non,dldj_knee_non,dlrjt_knee_non,dldj_hip_non,dlrjt_hip_non);

% SLDJ v SLRJT joint stiffness aclr
%sldj_ank_inj= Data.Ankle_Stiff_SLDJ_Inj;
%slrjt_ank_inj= Data.Ankle_Stiff_SLRJT_Inj;
%sldj_knee_inj= Data.Knee_Stiff_SLDJ_Inj;
%slrjt_knee_inj= Data.Knee_Stiff_SLRJT_Inj;
%sldj_hip_inj= Data.Hip_Stiff_SLDJ_Inj;
%slrjt_hip_inj= Data.Hip_Stiff_SLRJT_Inj;
%stat_analysis(sldj_ank_inj,slrjt_ank_inj,sldj_knee_inj,slrjt_knee_inj,sldj_hip_inj,slrjt_hip_inj);

% SLDJ v SLRJT joint stiffness nonaclr
%sldj_ank_non= Data.Ankle_Stiff_SLDJ_Non;
%slrjt_ank_non= Data.Ankle_Stiff_SLRJT_Non;
%sldj_knee_non= Data.Knee_Stiff_SLDJ_Non;
%slrjt_knee_non= Data.Knee_Stiff_SLRJT_Non;
%sldj_hip_non= Data.Hip_Stiff_SLDJ_Non;
%slrjt_hip_non= Data.Hip_Stiff_SLRJT_Non;
%stat_analysis(sldj_ank_non,slrjt_ank_non,sldj_knee_non,slrjt_knee_non,sldj_hip_non,slrjt_hip_non);



%COM Stiffness
dldj_stiff= Data.COM_Stiff_DLDJ;
dlrjt_stiff=Data.COM_Stiff_RJT;
sldj_stiff_acl= Data.COM_SLDJ_Inj; 
slrjt_stiff_acl=Data.COM_SLRJT_Inj;
sldj_stiff_nonacl=Data.COM_SLDJ_Non;
slrjt_stiff_nonacl=Data.COM_SLRJT_Non;
stat_analysis(dldj_stiff,dlrjt_stiff,sldj_stiff_acl,slrjt_stiff_acl,sldj_stiff_nonacl,slrjt_stiff_nonacl);

