
% based on Zaixu's code "Step_8th_Visualize_Workbench_Atlas"
%% Initialize SPM and load in files
%

% Clear workspace
clear

% Add paths to Connectome Workbench and SPM (if needed)
addpath('C:\workbench\bin_windows64'); % Add the directory containing Connectome Workbench binaries
addpath('C:\workbench\spm12'); % Add the directory containing SPM (if needed)

% Initialize SPM (if needed)
spm('Defaults','fmri');
spm_jobman('initcfg');

% Read in PFN_L and PFN_R
PFN_L=load('Rel_max_maps/Max_Neg_PRS_2_L_clust_5.mat');
PFN_R=load('Rel_max_maps/Max_Neg_PRS_2_R_clust_5.mat');


VisualizeFolder = ['Rel_max_maps'];
Map =['Max_Neg_PRS_2_clust_5']

sbj_AtlasLabel_lh = PFN_L.x.vertex_ID_lh;
sbj_AtlasLabel_rh = PFN_R.x.vertex_ID_rh;


%% PFN ID labeling
ColorInfo_Atlas = [VisualizeFolder '\name_Atlas.txt'];
system(['del ' ColorInfo_Atlas]);

SystemName = {'Net 1: DM', 'Net 2: SM', 'Net 3: FP', 'Net 4: SM', 'Net 5: DA', ...
              'Net 6: VS', 'Net 7: VA', 'Net 8: DM', 'Net 9: VA', 'Net 10: VS', 'Net 11: SM', ...
              'Net 12: DM', 'Net 13: SM', 'Net 14: DA', 'Net 15: FP', 'Net 16: AU', 'Net 17: FP'};
ColorPlate = {'242 139 168', '173 216 230', '244 197 115', '73 143 191', ...
              '65 171 93', '137 63 153', '217 117 242', '226 57 93', ...
              '206 28 249', '102 5 122', '33 113 181', '170 12 61', ...
              '7 69 132', '0 109 44', '216 144 72', '78 49 168', '204 109 14'};
for i = 1:17
  system(['echo ' SystemName{i} ' >> ' ColorInfo_Atlas]);
  system(['echo ' num2str(i) ' ' ColorPlate{i} ' 1 >> ' ColorInfo_Atlas]); 
end

% left hemi
V_lh = gifti;
V_lh.cdata = sbj_AtlasLabel_lh;
V_lh_File = [VisualizeFolder '/' Map '_lh.func.gii'];
save(V_lh, V_lh_File);
pause(1);
V_lh_Label_File = [VisualizeFolder '/' Map '_lh_AtlasLabel.label.gii'];
cmd = ['wb_command -metric-label-import ' V_lh_File ' ' ColorInfo_Atlas ' ' V_lh_Label_File];
system(cmd);
% right hemi
V_rh = gifti;
V_rh.cdata = sbj_AtlasLabel_rh;
V_rh_File = [VisualizeFolder '/' Map '_rh.func.gii'];
save(V_rh, V_rh_File);
pause(1);
V_rh_Label_File = [VisualizeFolder '/' Map '_rh_AtlasLabel.label.gii'];
cmd = ['wb_command -metric-label-import ' V_rh_File ' ' ColorInfo_Atlas ' ' V_rh_Label_File];
system(cmd);
% convert into cifti file
cmd = ['wb_command -cifti-create-label ' VisualizeFolder '/' Map '_AtlasLabel' ...
       '.dlabel.nii -left-label ' V_lh_Label_File ' -right-label ' V_rh_Label_File];
system(cmd);
pause(1);
%system(['rm -rf ' V_lh_File ' ' V_rh_File ' ' V_lh_Label_File ' ' V_rh_Label_File]);
