function SL=SL_report(SL,out_name,subjN)
% Setup the QA directory
QA_dir=fullfile(SL.dir.outpath,SL.dir.QA,filesep);
if ~exist(QA_dir,'dir'), mkdir(QA_dir); end
 
% File specific information 
sdisp('I/O',1);
display(['Reading from: ' SL.dir.stpath]);
display(['Saving to   : ' SL.dir.outpath]);
switch SL.dir.overwrite
    case 0, display('Overwrite   : OFF');
    case 1, display('Overwrite   : ON');
end
switch SL.region.use_mask
    case 0, display('Masking     : OFF');
    case 1, display('Masking     : ON');
end

for ii=1:length(SL.design.save_str)
    if ~exist(fullfile(QA_dir,SL.design.save_str{ii}),'dir')
       mkdir(fullfile(QA_dir,SL.design.save_str{ii}));
    end
    sdisp(['Model ' n2sp(ii,2) ': Type ' SL.design.calc{ii} ' : ' SL.design.save_str{ii}],1);
    model_type=SL.design.calc{ii};
    if strcmp(model_type, 'Anova1'),    model_type='Anova'; end
    if strcmp(model_type, 'Identity1'), model_type='Anova'; end
    if strcmp(model_type, 'Spear'),     model_type='cont'; end
    if strcmp(model_type, 'Euclid'),    model_type='cont'; end
    if strcmp(model_type, 'Kendall'),   model_type='cont'; end
    RSA_pdm(SL.design.matrix{ii},SL.design.Box,...
        [SL.dir.subjects{subjN} '_design'],...
        fullfile(QA_dir,SL.design.save_str{ii},[SL.dir.subjects{subjN} '_design']));
    switch model_type
        case 'Anova'
            for jj=1:length(SL.design.anova{ii}.f)
                RSA_pdm(SL.design.anova{ii}.f{jj},SL.design.Box,...
                    [SL.dir.subjects{subjN} '_f' SL.design.anova{ii}.names{jj}],...
                    fullfile(QA_dir,SL.design.save_str{ii},[SL.dir.subjects{subjN} '_f', ...
                    SL.design.anova{ii}.names{jj}]));
            end
    end
end
     
if SL.region.noSL==0
    display('--Save Files S(1)--');
    display('DNE: Does not exist');
    display('FND: Found');
    display('--Save Files S(1)--');
    for ii=1:length(out_name)
        f=fullfile(SL.dir.outpath,SL.dir.subjects{1},[out_name{ii} '.img']);
        switch exist(f,'file')
            case 0, display(['DNE: ' f]);
            case 2, display(['FND: ' f]);
        end
    end
end               
                        



