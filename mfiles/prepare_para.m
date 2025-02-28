function para = prepare_para(para)

fprintf('Loading parameters...');tic

matObj = matfile([para.dir.load_kSpace_dir,para.dir.load_kSpace_name]);
varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name],'kSpace_info')
    para.Recon.nSMS = max(kSpace_info.phase_mod) + 1;
    para.kSpace_info = kSpace_info;
end

para.time = datestr(clock,'yymmdd_hhMMSS');

disp('RawData:')
disp([para.dir.load_kSpace_dir, para.dir.load_kSpace_name])

para.dir.save_recon_img_name= para.dir.load_kSpace_name;

if isempty(dir([pwd,'/RawData/']))
    mkdir([pwd,'/RawData/']);
end

para.dir.save_recon_img_dir = [para.dir.load_kSpace_dir(1:end-8),'ReconData/'];
if isempty(dir(para.dir.save_recon_img_dir))
    mkdir(para.dir.save_recon_img_dir)
end

if para.Recon.nSMS == 1
    para.Recon.type = '2D';
else
    para.Recon.type = 'separate SMS';
end

para.Recon.epsilon = eps('single');
para.Recon.break = 0;
para.step_size = 2;

para.Recon.non_steady_state_rays = 300;

para.CPUtime.load_para_time = toc;toc;fprintf('\n');

end
