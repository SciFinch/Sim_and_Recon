function [datasetInfoStruct] = GetDatasetInfo(dataStr)
%GETDATASETINFO [datasetInfoStruct] = GetDatasetInfo(dataStr)
%	This function stores details on the path- and filenames of all data
%	for each dataset compatible with the ME-MCIR code.
%	On providing an identification string, dataStr, this function
%	will return a struct of (mostly) strings, which are stratified
%	to pathnames and filenames to motion models, respiratory data,
%	and so on.
%
%	Currently supported [warning: ca:
%	 - 'stub' - generic testing dataset with no meaningful data
%	 - 'simple' - testing dataset with simple geometric phantom
%	 - 'thesis1' - volunteer data used in thesis (volunteer 1)
%	 - 'thesis2' - volunteer data used in thesis (volunteer 2)
%	 - 'thesis3' - volunteer data used in thesis (volunteer 3)
%	 - 'thesis4' - volunteer data used in thesis (volunteer 4)
%	 - 'CK1' - collab patient data from Christoph Kolbitsch (patient 1)
%	 - 'CK2' - collab patient data from Christoph Kolbitsch (patient 2)
%	NOTE: thesisN data exist with lesions in different positions. This
%	      code could be extended to include this (at present this is not
%	 	  the case, defaults have been selected instead).

% Warning: This is a long file!

% Universal paths and filenames
datasetInfoStruct.resultWatch.pn = 'C:\Users\db12\Dropbox\ResultWatch\';
datasetInfoStruct.reconSave.pn = [];
% At some point I need to formally check motion orientations for each case (in case of
% (x,y,z) -> (hf,ap,rl) mixups, and to be _completely_ certain once and for all).
% Until then, the following will report a warning if orientations /not/ checked,
% for motion fields (MF) and motion models (MM);
%   NOTE: I think the thesisN datasets have some orientation checking already,
%	      check the motion model files
datasetInfoStruct.flags.orientationChecked.MF.fwd = false;
datasetInfoStruct.flags.orientationChecked.MF.bck = false;
datasetInfoStruct.flags.orientationChecked.MF.inv = false;
datasetInfoStruct.flags.orientationChecked.MM.fwd = false;
datasetInfoStruct.flags.orientationChecked.MM.bck = false;
datasetInfoStruct.flags.orientationChecked.MM.inv = false;

thesisN_basePath = 'C:\Users\db12\Repository_WIN\Raw_Simulation_Data\';
simpleData_basePath = 'C:\Users\db12\Dropbox\Work_Docs\SimpleData\';
% Dataset-specific paths and filenames:
switch dataStr
    
    %% STUB
    % Info: This is mostly empty. Only to be used for testing or as template for new datasets.
    case 'stub'
        
        % Dataset ID:
        datasetInfoStruct.dataStr = 'stub';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = [];
        datasetInfoStruct.maps.attenuation.human = [];
        datasetInfoStruct.maps.attenuation.hardware = [];
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = [];
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 1;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        
        % Motion fields
        datasetInfoStruct.nDynamics = 1;
        datasetInfoStruct.refDynNum = 1;
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = [];
            datasetInfoStruct.motionFields.HF.fwd{it} = [];
            datasetInfoStruct.motionFields.AP.fwd{it} = [];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [];
            datasetInfoStruct.motionFields.HF.bck{it} = [];
            datasetInfoStruct.motionFields.AP.bck{it} = [];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = [];
            datasetInfoStruct.motionFields.HF.inv{it} = [];
            datasetInfoStruct.motionFields.AP.inv{it} = [];
        end
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = [];
        % - forward
        datasetInfoStruct.motionModels.all.fwd = [];
        datasetInfoStruct.motionModels.RL.fwd = [];
        datasetInfoStruct.motionModels.HF.fwd = [];
        datasetInfoStruct.motionModels.AP.fwd = [];
        % - backward
        datasetInfoStruct.motionModels.all.bck = [];
        datasetInfoStruct.motionModels.RL.bck = [];
        datasetInfoStruct.motionModels.HF.bck = [];
        datasetInfoStruct.motionModels.AP.bck = [];
        % - inverse
        datasetInfoStruct.motionModels.all.inv = [];
        datasetInfoStruct.motionModels.RL.inv = [];
        datasetInfoStruct.motionModels.HF.inv = [];
        datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = [];
        
        %% SIMPLE
        % Info: Simple geometric shaped phantoms with simple motion trajectories for testing code at runtime.
    case 'simple'
        datasetInfoStruct.dataStr = 'simple';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = [simpleData_basePath 'EmissionMap_x100.gipl'];
        datasetInfoStruct.maps.attenuation.human = [simpleData_basePath 'MuMap_x10000.gipl'];
        datasetInfoStruct.maps.attenuation.hardware = [simpleData_basePath 'MuMap_hardware_empty.gipl'];
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = [simpleData_basePath 'SimulatedPETData.mat'];
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 1;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        
        % Motion fields
        datasetInfoStruct.nDynamics = 10;
        datasetInfoStruct.refDynNum = 1;
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = [simpleData_basePath ];
            datasetInfoStruct.motionFields.HF.fwd{it} = [simpleData_basePath ];
            datasetInfoStruct.motionFields.AP.fwd{it} = [simpleData_basePath ];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [];
            datasetInfoStruct.motionFields.HF.bck{it} = [];
            datasetInfoStruct.motionFields.AP.bck{it} = [];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = [simpleData_basePath ];
            datasetInfoStruct.motionFields.HF.inv{it} = [simpleData_basePath ];
            datasetInfoStruct.motionFields.AP.inv{it} = [simpleData_basePath ];
        end
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = [simpleData_basePath 'SimpleMotionModel.mat'];
        % - forward
%         datasetInfoStruct.motionModels.all.fwd = [];
%         datasetInfoStruct.motionModels.RL.fwd = [];
%         datasetInfoStruct.motionModels.HF.fwd = [];
%         datasetInfoStruct.motionModels.AP.fwd = [];
%         % - backward
%         datasetInfoStruct.motionModels.all.bck = [];
%         datasetInfoStruct.motionModels.RL.bck = [];
%         datasetInfoStruct.motionModels.HF.bck = [];
%         datasetInfoStruct.motionModels.AP.bck = [];
%         % - inverse
%         datasetInfoStruct.motionModels.all.inv = [];
%         datasetInfoStruct.motionModels.RL.inv = [];
%         datasetInfoStruct.motionModels.HF.inv = [];
%         datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = [simpleData_basePath 'respsig_all.mat'];
        
        %% THESIS1
        % Info: First volunteer dataset, used in thesis. Good middle-range motion displacements.
    case 'thesis1'
        datasetInfoStruct.dataStr = 'thesis1';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = [thesisN_basePath 'DB\staticMaps\'];
        datasetInfoStruct.maps.attenuation.human = [thesisN_basePath '\DB\staticMaps\'];
        datasetInfoStruct.maps.attenuation.hardware = [];
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = [];
        
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 6;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        % Alternatively, pre-generated PET simulations:
        datasetInfoStruct.sims.preExisting = ['\\nas14\motion_repository\PET_sims\Sinograms\DB\21001\DB_' num2str(datasetInfoStruct.realPET.nGates) 'G_L3_optStr21001_PromptMaps.mat']
        
        % Motion fields
        datasetInfoStruct.nDynamics = 35;
        datasetInfoStruct.refDynNum = 29;
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = [thesisN_basePath 'DB_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx.gipl'];
            datasetInfoStruct.motionFields.HF.fwd{it} = [thesisN_basePath 'DB_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my.gipl'];
            datasetInfoStruct.motionFields.AP.fwd{it} = [thesisN_basePath 'DB_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz.gipl'];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [];
            datasetInfoStruct.motionFields.HF.bck{it} = [];
            datasetInfoStruct.motionFields.AP.bck{it} = [];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = [thesisN_basePath 'DB_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx_inv_x100.gipl'];
            datasetInfoStruct.motionFields.HF.inv{it} = [thesisN_basePath 'DB_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my_inv_x100.gipl'];
            datasetInfoStruct.motionFields.AP.inv{it} = [thesisN_basePath 'DB_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz_inv_x100.gipl'];
        end
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = [thesisN_basePath 'MotionModels\DB\' num2str(datasetInfoStruct.realPET.nGates) 'G\'];
        % - forward
        datasetInfoStruct.motionModels.all.fwd = [];
        datasetInfoStruct.motionModels.RL.fwd = [];
        datasetInfoStruct.motionModels.HF.fwd = [];
        datasetInfoStruct.motionModels.AP.fwd = [];
        % - backward
        datasetInfoStruct.motionModels.all.bck = [];
        datasetInfoStruct.motionModels.RL.bck = [];
        datasetInfoStruct.motionModels.HF.bck = [];
        datasetInfoStruct.motionModels.AP.bck = [];
        % - inverse
        datasetInfoStruct.motionModels.all.inv = [];
        datasetInfoStruct.motionModels.RL.inv = [];
        datasetInfoStruct.motionModels.HF.inv = [];
        datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = [thesisN_basePath 'MotionModels\DB\' num2str(datasetInfoStruct.realPET.nGates) 'G\DB_6G_RespSigInfo_v1.mat'];
        
        %% THESIS2
        % Info: Second volunteer dataset, used in thesis. Low-to-mid--range motion displacements.
    case 'thesis2'
        datasetInfoStruct.dataStr = 'thesis2';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = [];
        datasetInfoStruct.maps.attenuation.human = [];
        datasetInfoStruct.maps.attenuation.hardware = [];
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = [];
        
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 6;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        % Alternatively, pre-generated PET simulations:
        datasetInfoStruct.sims.preExisting = ['\\nas14\motion_repository\PET_sims\Sinograms\DP\21001\DP_' num2str(datasetInfoStruct.realPET.nGates) 'G_L1_optStr21001_PromptMaps.mat'];
        
        
        % Motion fields
        datasetInfoStruct.nDynamics = 35;
        datasetInfoStruct.refDynNum = 35;
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = [thesisN_basePath 'DP_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx.gipl'];
            datasetInfoStruct.motionFields.HF.fwd{it} = [thesisN_basePath 'DP_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my.gipl'];
            datasetInfoStruct.motionFields.AP.fwd{it} = [thesisN_basePath 'DP_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz.gipl'];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [];
            datasetInfoStruct.motionFields.HF.bck{it} = [];
            datasetInfoStruct.motionFields.AP.bck{it} = [];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = [thesisN_basePath 'DP_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx_inv_x100.gipl'];
            datasetInfoStruct.motionFields.HF.inv{it} = [thesisN_basePath 'DP_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my_inv_x100.gipl'];
            datasetInfoStruct.motionFields.AP.inv{it} = [thesisN_basePath 'DP_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz_inv_x100.gipl'];
        end
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = [thesisN_basePath 'MotionModels\DP\' num2str(datasetInfoStruct.realPET.nGates) 'G\'];
        % - forward
        datasetInfoStruct.motionModels.all.fwd = [];
        datasetInfoStruct.motionModels.RL.fwd = [];
        datasetInfoStruct.motionModels.HF.fwd = [];
        datasetInfoStruct.motionModels.AP.fwd = [];
        % - backward
        datasetInfoStruct.motionModels.all.bck = [];
        datasetInfoStruct.motionModels.RL.bck = [];
        datasetInfoStruct.motionModels.HF.bck = [];
        datasetInfoStruct.motionModels.AP.bck = [];
        % - inverse
        datasetInfoStruct.motionModels.all.inv = [];
        datasetInfoStruct.motionModels.RL.inv = [];
        datasetInfoStruct.motionModels.HF.inv = [];
        datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = [thesisN_basePath 'MotionModels\DE\' num2str(datasetInfoStruct.realPET.nGates) 'G\DE_6G_RespSigInfo_v1.mat'];
        
        %% THESIS1
        % Info: Third volunteer dataset, used in thesis. Fairly small ranges of motion displacements.
    case 'thesis3'
        datasetInfoStruct.dataStr = 'thesis3';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = [];
        datasetInfoStruct.maps.attenuation.human = [];
        datasetInfoStruct.maps.attenuation.hardware = [];
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = [];
        
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 6;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        % Alternatively, pre-generated PET simulations:
        datasetInfoStruct.sims.preExisting = ['\\nas14\motion_repository\PET_sims\Sinograms\HT\21001\HT_' num2str(datasetInfoStruct.realPET.nGates) 'G_L1_optStr21001_PromptMaps.mat']
        
        % Motion fields
        datasetInfoStruct.nDynamics = 35;
        datasetInfoStruct.refDynNum = 10;
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = [thesisN_basePath 'HT_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx.gipl'];
            datasetInfoStruct.motionFields.HF.fwd{it} = [thesisN_basePath 'HT_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my.gipl'];
            datasetInfoStruct.motionFields.AP.fwd{it} = [thesisN_basePath 'HT_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz.gipl'];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [thesisN_basePath 'HT_forWindows\nb' num2str(it) 'tonb' num2str(datasetInfoStruct.refDynNum) '_mx.gipl'];
            datasetInfoStruct.motionFields.HF.bck{it} = [thesisN_basePath 'HT_forWindows\nb' num2str(it) 'tonb' num2str(datasetInfoStruct.refDynNum) '_my.gipl'];
            datasetInfoStruct.motionFields.AP.bck{it} = [thesisN_basePath 'HT_forWindows\nb' num2str(it) 'tonb' num2str(datasetInfoStruct.refDynNum) '_mz.gipl'];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = [thesisN_basePath 'HT_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx_inv_x100.gipl'];
            datasetInfoStruct.motionFields.HF.inv{it} = [thesisN_basePath 'HT_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my_inv_x100.gipl'];
            datasetInfoStruct.motionFields.AP.inv{it} = [thesisN_basePath 'HT_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz_inv_x100.gipl'];
        end
        
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = [thesisN_basePath 'MotionModels\HT\' num2str(datasetInfoStruct.realPET.nGates) 'G\'];
        % - forward
        datasetInfoStruct.motionModels.all.fwd = [];
        datasetInfoStruct.motionModels.RL.fwd = [];
        datasetInfoStruct.motionModels.HF.fwd = [];
        datasetInfoStruct.motionModels.AP.fwd = [];
        % - backward
        datasetInfoStruct.motionModels.all.bck = [];
        datasetInfoStruct.motionModels.RL.bck = [];
        datasetInfoStruct.motionModels.HF.bck = [];
        datasetInfoStruct.motionModels.AP.bck = [];
        % - inverse
        datasetInfoStruct.motionModels.all.inv = [];
        datasetInfoStruct.motionModels.RL.inv = [];
        datasetInfoStruct.motionModels.HF.inv = [];
        datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = [thesisN_basePath 'MotionModels\HT\' num2str(datasetInfoStruct.realPET.nGates) 'G\HT_6G_RespSigInfo_v1.mat'];
        
        %% THESIS4
        % Info: Fourth volunteer dataset, used in thesis. Fairly large motion displacements.
        %		Has up to 70 motion positions recorded.
    case 'thesis4'
        datasetInfoStruct.dataStr = 'thesis4';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = [];
        datasetInfoStruct.maps.attenuation.human = [];
        datasetInfoStruct.maps.attenuation.hardware = [];
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = [];
        
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 6;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        % Alternatively, pre-generated PET simulations:
        datasetInfoStruct.sims.preExisting = ['\\nas14\motion_repository\PET_sims\Sinograms\DE\21001\DE_' num2str(datasetInfoStruct.realPET.nGates) 'G_L3_optStr21001_PromptMaps.mat']
        
        % Motion fields
        datasetInfoStruct.nDynamics = 70;
        datasetInfoStruct.refDynNum = 36;
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = [thesisN_basePath 'DE_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx.gipl'];
            datasetInfoStruct.motionFields.HF.fwd{it} = [thesisN_basePath 'DE_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my.gipl'];
            datasetInfoStruct.motionFields.AP.fwd{it} = [thesisN_basePath 'DE_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz.gipl'];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [thesisN_basePath 'DE_forWindows\nb' num2str(it) 'tonb' num2str(datasetInfoStruct.refDynNum) '_mx.gipl'];
            datasetInfoStruct.motionFields.HF.bck{it} = [thesisN_basePath 'DE_forWindows\nb' num2str(it) 'tonb' num2str(datasetInfoStruct.refDynNum) '_my.gipl'];
            datasetInfoStruct.motionFields.AP.bck{it} = [thesisN_basePath 'DE_forWindows\nb' num2str(it) 'tonb' num2str(datasetInfoStruct.refDynNum) '_mz.gipl'];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = [thesisN_basePath 'DE_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mx_inv_x100.gipl'];
            datasetInfoStruct.motionFields.HF.inv{it} = [thesisN_basePath 'DE_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_my_inv_x100.gipl'];
            datasetInfoStruct.motionFields.AP.inv{it} = [thesisN_basePath 'DE_forWindows\g' num2str(datasetInfoStruct.refDynNum) 'tog' num2str(it) '_mz_inv_x100.gipl'];
        end
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = [thesisN_basePath 'MotionModels\DE\' num2str(datasetInfoStruct.realPET.nGates) 'G\'];
        % - forward
        datasetInfoStruct.motionModels.all.fwd = [];
        datasetInfoStruct.motionModels.RL.fwd = [];
        datasetInfoStruct.motionModels.HF.fwd = [];
        datasetInfoStruct.motionModels.AP.fwd = [];
        % - backward
        datasetInfoStruct.motionModels.all.bck = [];
        datasetInfoStruct.motionModels.RL.bck = [];
        datasetInfoStruct.motionModels.HF.bck = [];
        datasetInfoStruct.motionModels.AP.bck = [];
        % - inverse
        datasetInfoStruct.motionModels.all.inv = [];
        datasetInfoStruct.motionModels.RL.inv = [];
        datasetInfoStruct.motionModels.HF.inv = [];
        datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = [thesisN_basePath 'MotionModels\DE\' num2str(datasetInfoStruct.realPET.nGates) 'G\DE_6G_RespSigInfo_v1.mat'];
        
        %% CK1
        % Info: First patient dataset provided by CK, also processed for simulation. Small range of motion displacements.
    case 'CK1'
        datasetInfoStruct.dataStr = 'CK1';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = 'C:\Users\db12\Repository_WIN\RealDataAsSim\FDGmap.gipl';
        datasetInfoStruct.maps.attenuation.human = 'C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\PET_ACQ_51_20160201134327_umap_human_00_RPE.v.hdr';
        datasetInfoStruct.maps.attenuation.hardware = 'C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\PET_ACQ_51_20160201134327_umap_hardware_00.v.hdr';
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = 'C:\Users\db12\Repository_WIN\RealDataAsSim\SimData_wAtt.mat';
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 1;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        
        % Motion fields
        datasetInfoStruct.nDynamics = 8;
        datasetInfoStruct.refDynNum = 1;
        warning('orientations of CK1_sim MFs have not been checked')
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = ['C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Motion Fields\CKdat_gMR_g' num2str(datasetInfoStruct.refDynNum) '_to_g' num2str(it) '_mx.gipl'];
            datasetInfoStruct.motionFields.HF.fwd{it} = ['C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Motion Fields\CKdat_gMR_g' num2str(datasetInfoStruct.refDynNum) '_to_g' num2str(it) '_my.gipl'];
            datasetInfoStruct.motionFields.AP.fwd{it} = ['C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Motion Fields\CKdat_gMR_g' num2str(datasetInfoStruct.refDynNum) '_to_g' num2str(it) '_mz.gipl'];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [];
            datasetInfoStruct.motionFields.HF.bck{it} = [];
            datasetInfoStruct.motionFields.AP.bck{it} = [];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = ['C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Motion Fields\CKdat_gMR_g' num2str(datasetInfoStruct.refDynNum) '_to_g' num2str(it) '_inv_x1000_mx.gipl'];
            datasetInfoStruct.motionFields.HF.inv{it} = ['C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Motion Fields\CKdat_gMR_g' num2str(datasetInfoStruct.refDynNum) '_to_g' num2str(it) '_inv_x1000_mx.gipl'];
            datasetInfoStruct.motionFields.AP.inv{it} = ['C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Motion Fields\CKdat_gMR_g' num2str(datasetInfoStruct.refDynNum) '_to_g' num2str(it) '_inv_x1000_mx.gipl'];
        end
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = 'C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\MotionModel and RespSigs\RespSigs.mat';
        % - forward
        datasetInfoStruct.motionModels.all.fwd = [];
        datasetInfoStruct.motionModels.RL.fwd = [];
        datasetInfoStruct.motionModels.HF.fwd = [];
        datasetInfoStruct.motionModels.AP.fwd = [];
        % - backward
        datasetInfoStruct.motionModels.all.bck = [];
        datasetInfoStruct.motionModels.RL.bck = [];
        datasetInfoStruct.motionModels.HF.bck = [];
        datasetInfoStruct.motionModels.AP.bck = [];
        % - inverse
        datasetInfoStruct.motionModels.all.inv = [];
        datasetInfoStruct.motionModels.RL.inv = [];
        datasetInfoStruct.motionModels.HF.inv = [];
        datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = 'C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\NavSig.mat';
        
        %% CK2
        % Info: Second patient dataset provided by CK, processed for simulation. ? range of motion displacements.
    case 'CK2'
        datasetInfoStruct.dataStr = 'CK2';
        % Distribution maps from segmented UTE volumes:
        datasetInfoStruct.maps.emission.fn = [];
        datasetInfoStruct.maps.attenuation.fn = [];
        % Simulated data (if not generated, will be placed here)
        datasetInfoStruct.sims.fn = [];
        % Real PET data (if it exists)
        datasetInfoStruct.realPET.nGates = 1;
        for it = 1:datasetInfoStruct.realPET.nGates
            datasetInfoStruct.realPET.prompts{it} = [];
            datasetInfoStruct.realPET.ACFs{it} = [];
            datasetInfoStruct.realPET.NCFs{it} = [];
            datasetInfoStruct.realPET.scatters{it} = [];
            datasetInfoStruct.realPET.randoms{it} = [];
        end
        
        % Motion fields
        datasetInfoStruct.nDynamics = [];
        % - forward (ref -> postion)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.fwd{it} = [];
            datasetInfoStruct.motionFields.HF.fwd{it} = [];
            datasetInfoStruct.motionFields.AP.fwd{it} = [];
        end
        % - backward (position -> ref)
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.bck{it} = [];
            datasetInfoStruct.motionFields.HF.bck{it} = [];
            datasetInfoStruct.motionFields.AP.bck{it} = [];
        end
        % - inverse ( inv(ref -> position) )
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.motionFields.RL.inv{it} = [];
            datasetInfoStruct.motionFields.HF.inv{it} = [];
            datasetInfoStruct.motionFields.AP.inv{it} = [];
        end
        
        % Motion Modelling
        % - all (if exists)
        datasetInfoStruct.motionModels.all = [];
        % - forward
        datasetInfoStruct.motionModels.all.fwd = [];
        datasetInfoStruct.motionModels.RL.fwd = [];
        datasetInfoStruct.motionModels.HF.fwd = [];
        datasetInfoStruct.motionModels.AP.fwd = [];
        % - backward
        datasetInfoStruct.motionModels.all.bck = [];
        datasetInfoStruct.motionModels.RL.bck = [];
        datasetInfoStruct.motionModels.HF.bck = [];
        datasetInfoStruct.motionModels.AP.bck = [];
        % - inverse
        datasetInfoStruct.motionModels.all.inv = [];
        datasetInfoStruct.motionModels.RL.inv = [];
        datasetInfoStruct.motionModels.HF.inv = [];
        datasetInfoStruct.motionModels.AP.inv = [];
        
        % Dynamic MR image volumes
        for it = 1:datasetInfoStruct.nDynamics
            datasetInfoStruct.dynamicMR.fn{it} = [];
        end
        datasetInfoStruct.dynamicMR.ref = [];
        
        % Respiratory signals
        datasetInfoStruct.respSig.fn = [];
    otherwise
        error('Unrecognised dataset requested')
end

end