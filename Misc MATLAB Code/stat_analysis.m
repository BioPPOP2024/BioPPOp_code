% %SLDJ & SLRJT performance variables
function stat_analysis(Inj_DJ,Non_DJ,Inj_RJT,Non_RJT,Asya,Asyb)
answer = inputdlg('ENTER NAME OF COMPARISON:',...
             'Sample', [1 50]);
Name=cell2mat(answer);
%t-test of performance metric SLDJ
[h,p,ci,stats] = ttest(Inj_DJ,Non_DJ);
d = computeCohen_d(Inj_DJ,Non_DJ,'paired');
xlswrite('Stats',h,Name,'A1')
xlswrite('Stats',p,Name,'A2')
xlswrite('Stats',ci,Name,'A3')
xlswrite('Stats',d,Name,'A6')
averinj=mean(Inj_DJ);
avernon=mean(Non_DJ);
sdinj=std(Inj_DJ);
sdnon=std(Non_DJ);
xlswrite('Stats',averinj,Name,'C1')
xlswrite('Stats',avernon,Name,'C2')
xlswrite('Stats',sdinj,Name,'D1')
xlswrite('Stats',sdnon,Name,'D2')

answer2 = inputdlg('ENTER NAME OF COMPARISON:',...
             'Sample', [1 50]);
Name2=cell2mat(answer2);
%t-test of performance metric SLRJT
[h,p,ci,stats] = ttest(Inj_RJT,Non_RJT);
d = computeCohen_d(Inj_RJT,Non_RJT,'paired');
xlswrite('Stats',h,Name2,'A1')
xlswrite('Stats',p,Name2,'A2')
xlswrite('Stats',ci,Name2,'A3')
xlswrite('Stats',d,Name2,'A6')
averinj=mean(Inj_RJT);
avernon=mean(Non_RJT);
sdinj=std(Inj_RJT);
sdnon=std(Non_RJT);
xlswrite('Stats',averinj,Name2,'C1')
xlswrite('Stats',avernon,Name2,'C2')
xlswrite('Stats',sdinj,Name2,'D1')
xlswrite('Stats',sdnon,Name2,'D2')
%Assymetry sldj & slrjt
answer3 = inputdlg('ENTER NAME OF COMPARISON:',...
             'Sample', [1 50]);
Name3=cell2mat(answer3);
[h,p,ci,stats] = ttest(Asya,Asyb);
d = computeCohen_d(Asya,Asyb,'paired');
xlswrite('Stats',h,Name3,'A1')
xlswrite('Stats',p,Name3,'A2')
xlswrite('Stats',ci,Name3,'A3')
xlswrite('Stats',d,Name3,'A6')
aversldj=mean(Asya)
averslrjt=mean(Asyb)
sdsldj=std(Asya)
sdslrjt=std(Asyb)
xlswrite('Stats',aversldj,Name3,'C1')
xlswrite('Stats',averslrjt,Name3,'C2')
xlswrite('Stats',sdsldj,Name3,'D1')
xlswrite('Stats',sdslrjt,Name3,'D2')


end