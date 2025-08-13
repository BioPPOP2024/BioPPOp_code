%read data
Data=readtable('Moments and performance variables New.xlsx');
%factors
injleg=ones(1,21)'
nonleg=[ones(1,21)*2]'
leg=[injleg;nonleg;injleg;nonleg]
dj=[ones(1,42)*3]'
rjt=[ones(1,42)*4]'
jump=[dj;rjt]

%Single-leg jump stats

%Jump Height
SLDJinj_JH= Data.JumpHeight_SLDJ_Inj;
SLDJnon_JH= Data.JumpHeight_SLDJ_Non;
SLRJTinj_JH= Data.JumpHeight_SLRJT_Inj;
SLRJTnon_JH= Data.JumpHeight_SLRJT_Non;
SLDJ=[Data.JumpHeight_SLDJ_Inj;Data.JumpHeight_SLDJ_Non];
SLRJT= [Data.JumpHeight_SLRJT_Inj; Data.JumpHeight_SLRJT_Non];
SL_JH=[SLDJ;SLRJT];
[p,tbl,stats]=anovan(SL_JH,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause

%Contact Time
SLDJinj_CT= Data.CT_SLDJ_Inj;
SLDJnon_CT= Data.CT_SLDJ_Non;
SLRJTinj_CT= Data.CT_SLRJT_Inj;
SLRJTnon_CT= Data.CT_SLRJT_Non;
SL_CT = [SLDJinj_CT;SLDJnon_CT;SLRJTinj_CT;SLRJTnon_CT];
[p,tbl,stats]=anovan(SL_CT,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause

%RSI
SLDJinj_RSI= Data.RSI_SLDJ_Inj;
SLDJnon_RSI= Data.RSI_SLDJ_Non;
SLRJTinj_RSI= Data.RSI_SLRJT_Inj;
SLRJTnon_RSI= Data.RSI_SLRJT_Non;
SL_RSI = [SLDJinj_RSI;SLDJnon_RSI;SLRJTinj_RSI;SLRJTnon_RSI];
[p,tbl,stats]=anovan(SL_RSI,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause
%Double-leg jump stats
DLRJT_JH= Data.JumpHeight_RJT;
DLDJ_JH= Data.JumpHeight_DLDJ;
[h1,p1,stats1]=ttest(DLDJ_JH,DLRJT_JH,'Alpha',0.01);
pause
DLDJ_CT= Data.CT_DLDJ;
DLRJT_CT= Data.CT_RJT;
[h2,p2,stats2]=ttest(DLDJ_CT,DLRJT_CT,'Alpha',0.01);
pause
DLDJ_RSI= Data.RSI_DLDJ;
DLRJT_RSI= Data.RSI_RJT;
[h3,p3,stats3]=ttest(DLDJ_RSI,DLRJT_RSI,'Alpha',0.01);

%lsi-jh
SLDJ_JH_LSI=Data.JumpHeight_SLDJ_Asy;
SLRJT_JH_LSI= Data.JumpHeight_SLRJT_Asy;
[h4,p4,stats4]=ttest(SLDJ_JH_LSI,SLRJT_JH_LSI);
%lsi-ct
SLDJ_CT_LSI=Data.CT_SLDJ_Asy;
SLRJT_CT_LSI= Data.CT_SLRJT_Asy;
[h5,p5,stats5]=ttest(SLDJ_CT_LSI,SLRJT_CT_LSI);

%lsi-rsi
SLDJ_RSI_LSI=Data.RSI_SLDJ_Asy;
SLRJT_RSI_LSI= Data.RSI_SLRJT_Asy;
[h6,p6,stats6]=ttest(SLDJ_RSI_LSI,SLRJT_RSI_LSI);
pause

%Vertical Stiffness
%Single-leg Jumps
SLDJinj_COM= Data.COM_SLDJ_Inj;
SLDJnon_COM= Data.COM_SLDJ_Non;
SLRJTinj_COM= Data.COM_SLRJT_Inj;
SLRJTnon_COM= Data.COM_SLRJT_Non;

SL_COM=[SLDJinj_COM;SLDJnon_COM;SLRJTinj_COM;SLRJTnon_COM];
[p,tbl,stats]=anovan(SL_COM,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause
%Double-leg Jumps
DLDJ_COM= Data.COM_Stiff_DLDJ;
DLRJT_COM= Data.COM_Stiff_RJT;
[h7,p7,stats7]=ttest(DLDJ_COM,DLRJT_COM);

%Ankle Stiff
%Single-leg Jumps
SLDJinj_ANK= Data.Ankle_Stiff_SLDJ_Inj;
SLDJnon_ANK= Data.Ankle_Stiff_SLDJ_Non;
SLRJTinj_ANK= Data.Ankle_Stiff_SLRJT_Inj;
SLRJTnon_ANK= Data.Ankle_Stiff_SLRJT_Non;
SL_ANK=[SLDJinj_ANK;SLDJnon_ANK;SLRJTinj_ANK;SLRJTnon_ANK];
[p,tbl,stats]=anovan(SL_ANK,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause
%Double-leg Jumps
DLDJinj_ANK= Data.Ankle_Stiff_DLDJ_Inj;
DLDJnon_ANK= Data.Ankle_Stiff_DLDJ_Non;
DLRJTinj_ANK= Data.Ankle_Stiff_RJT_Inj;
DLRJTnon_ANK= Data.Ankle_Stiff_RJT_Non;
DL_ANK=[DLDJinj_ANK;DLDJnon_ANK;DLRJTinj_ANK;DLRJTnon_ANK];
[p,tbl,stats]=anovan(DL_ANK,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause

%Knee Stiff
%Single-leg Jumps
SLDJinj_KNE= Data.Knee_Stiff_SLDJ_Inj;
SLDJnon_KNE= Data.Knee_Stiff_SLDJ_Non;
SLRJTinj_KNE= Data.Knee_Stiff_SLRJT_Inj;
SLRJTnon_KNE= Data.Knee_Stiff_SLRJT_Non;
SL_KNE=[SLDJinj_KNE;SLDJnon_KNE;SLRJTinj_KNE;SLRJTnon_KNE];
[p,tbl,stats]=anovan(SL_KNE,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause


%Double-leg Jumps
DLDJinj_KNE= Data.Knee_Stiff_DLDJ_Inj;
DLDJnon_KNE= Data.Knee_Stiff_DLDJ_Non;
DLRJTinj_KNE= Data.Knee_Stiff_RJT_Inj;
DLRJTnon_KNE= Data.Knee_Stiff_RJT_Non;
DL_KNE=[DLDJinj_KNE;DLDJnon_KNE;DLRJTinj_KNE;DLRJTnon_KNE];
[p,tbl,stats]=anovan(DL_KNE,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause
%Hip Stiff
%Single-leg Jumps
SLDJinj_HIP= Data.Hip_Stiff_SLDJ_Inj;
SLDJnon_HIP= Data.Hip_Stiff_SLDJ_Non;
SLRJTinj_HIP= Data.Hip_Stiff_SLRJT_Inj;
SLRJTnon_HIP= Data.Hip_Stiff_SLRJT_Non;
SL_HIP=[SLDJinj_HIP;SLDJnon_HIP;SLRJTinj_HIP;SLRJTnon_HIP];
[p,tbl,stats]=anovan(SL_HIP,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause
%Double-leg Jumps
DLDJinj_HIP= Data.Hip_Stiff_DLDJ_Inj;
DLDJnon_HIP= Data.Hip_Stiff_DLDJ_Non;
DLRJTinj_HIP= Data.Hip_Stiff_RJT_Inj;
DLRJTnon_HIP= Data.Hip_Stiff_RJT_Non;
DL_HIP=[DLDJinj_HIP;DLDJnon_HIP;DLRJTinj_HIP;DLRJTnon_HIP];
[p,tbl,stats]=anovan(DL_HIP,{leg,jump},'model','interaction');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2])
pause