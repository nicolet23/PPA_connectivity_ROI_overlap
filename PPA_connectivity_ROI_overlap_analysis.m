% a few to dos
% 1) Exclude anterior and posterior ppa (i.e., the seed regions) from all calulations ✅
% 2) Get the average correlation values, Dice coeffients, etc from individual regions of interest ✅
% 3) Do this as a loop for all subjects and plot bar plots of results. ✅


addpath('~/Documents/MATLAB/afni_matlab/matlab/')

%% Helper functions 
function M = safe_read_1D(p)
    if exist('Read_1D','file')==2
        M = double(Read_1D(p));
    else
        % Fallback: skip lines starting with #
        fid = fopen(p, 'r');
        lines = {};
        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line), break; end
            if ~isempty(line) && line(1) ~= '#'
                lines{end+1} = line;
            end
        end
        fclose(fid);
        M = str2num(char(lines));
    end
end

function t = pct(v, p)
    v = v(isfinite(v(:)));
    if isempty(v), t = NaN; return; end
    v = sort(v);
    t = v(max(1, min(numel(v), ceil(p/100 * numel(v)))));
end

%% Configuration
subjects = {'scim001','scim002','scim008','scim012','scim014','scim076', ...
            'scim257','scim297','scim312','scim357','scim0989','scim995'};
hemis = {'lh','rh'};
target_rois = [39, 40, 41, 42];
activation_col = 10;
prc = 95;

proj_root = '/Users/nicoletang/Desktop/SteelLab';
roi_directory_root = fullfile(proj_root, 'projects/Nicole_PlaceMemory/data');

AllResults = struct();
ResultsTable = [];

%% Main loop
for subj_idx = 1:length(subjects)
    subj = subjects{subj_idx};
    
    for hemi_idx = 1:length(hemis)
        hemi = hemis{hemi_idx};
        
        fprintf('\n>>> %s - %s\n', subj, hemi);
        
        try
            %% Load data
            corr_1d_path = fullfile(proj_root, 'anterior_posterior/Adam/corrmaps/ppa.anteriorposterior', ...
                sprintf('%s.ppa.anteriorposterior.CONTRAST.smooth.GSR.%s.niml.dset.%s.1D.dset', subj, hemi, hemi));
            
            mimetic_dir = fullfile(roi_directory_root, subj, 'mimetic');
            act_files = dir(fullfile(mimetic_dir, sprintf('*.fullRun.MIMETIC.bucket_GAM.%s.1D.dset', hemi)));
            if isempty(act_files), continue; end
            
            cdat = safe_read_1D(corr_1d_path);
            adat = safe_read_1D(fullfile(mimetic_dir, act_files(1).name));
            
            %% Process data
            if size(cdat,2)==1, cdat = [(1:size(cdat,1))' cdat]; end
            if size(adat,2)==1, adat = [(1:size(adat,1))' adat]; end
            
            [idx_common, ic, ia] = intersect(cdat(:,1), adat(:,1), 'stable');
            c_use = cdat(ic, min(2,size(cdat,2)));
            a_use = adat(ia, activation_col);
            
            valid = isfinite(c_use) & isfinite(a_use);
            c_use = c_use(valid);
            a_use = a_use(valid);
            idx_common = idx_common(valid);
            
            %% Exclude seed regions
            seed_roi_path = fullfile(proj_root, 'anterior_posterior/Adam', subj, ...
                sprintf('%s_%s.ppa.anteriorposterior.1D.roi', subj, hemi));
            
            if exist(seed_roi_path, 'file')
                seed_dat = safe_read_1D(seed_roi_path);
                seed_nodes = seed_dat(ismember(seed_dat(:,2), [1, 2]), 1);
                exclude = ismember(idx_common, seed_nodes);
                idx_common = idx_common(~exclude);
                c_use = c_use(~exclude);
                a_use = a_use(~exclude);
                fprintf('   Excluded %d seeds (N=%d)\n', nnz(exclude), numel(idx_common));
            end
            
            %% Load ROI
            roi_path = fullfile(roi_directory_root, subj, 'mimetic', sprintf('ref.nicole-final-roi-%s.1D.roi', hemi));
            if ~exist(roi_path, 'file')
                roi_path = fullfile(roi_directory_root, subj, 'dynloc', sprintf('ref.nicole-final-roi-%s.1D.roi', hemi));
            end
            if ~exist(roi_path, 'file'), continue; end
            
            roi_dat = safe_read_1D(roi_path);
            r_idx = roi_dat(:,1);
            r_lab = roi_dat(:,2);
            
            %% Calculate metrics
            top_mask = c_use >= pct(c_use, prc);
            
            for roi_id = unique(r_lab(r_lab > 0))'
                [~, pos] = ismember(idx_common, r_idx);
                roi_mask = false(numel(idx_common), 1);
                hit = pos > 0;
                roi_mask(hit) = (r_lab(pos(hit)) == roi_id);
                
                if nnz(roi_mask) == 0, continue; end
                
                inter = nnz(top_mask & roi_mask);
                roi_size = nnz(roi_mask);
                mean_corr = mean(c_use(roi_mask), 'omitnan');
                median_corr = median(c_use(roi_mask), 'omitnan');
                dice_roi = 2*inter/(nnz(top_mask)+nnz(roi_mask));
                overlap_pct = 100 * inter / roi_size;
                
                AllResults.(subj).(hemi).(['ROI' num2str(roi_id)]) = struct(...
                    'size', roi_size, 'mean', mean_corr, 'median', median_corr, ...
                    'dice', dice_roi, 'overlap', overlap_pct);
                
                ResultsTable = [ResultsTable; {subj, hemi, roi_id, roi_size, ...
                    mean_corr, median_corr, dice_roi, overlap_pct}];
            end
            
        catch
            continue
        end
    end
end

%% Save results
ResultsTable = cell2table(ResultsTable, 'VariableNames', ...
    {'Subject','Hemisphere','ROI','Size','MeanCorr','MedianCorr','Dice','OverlapPct'});

save('AllSubjects_ROI_Analysis.mat', 'AllResults', 'ResultsTable');
writetable(ResultsTable, 'AllSubjects_ROI_Analysis.csv');
fprintf('\n✅ Analysis complete\n');

%% Individual subject plots
for roi_id = target_rois
    for hemi_str = {'lh', 'rh'}
        roi_data = ResultsTable(ResultsTable.ROI == roi_id & strcmp(ResultsTable.Hemisphere, hemi_str{1}), :);
        if isempty(roi_data), continue; end
        
        figure('Name', sprintf('ROI%d %s', roi_id, upper(hemi_str{1})), 'Position', [100 100 1400 900]);
        
        subplot(2,2,1); bar(roi_data.MeanCorr);
        set(gca, 'XTick', 1:height(roi_data), 'XTickLabel', roi_data.Subject, 'XTickLabelRotation', 45);
        ylabel('Mean Corr'); title(sprintf('ROI%d %s: Mean', roi_id, upper(hemi_str{1}))); grid on;
        
        subplot(2,2,2); bar(roi_data.MedianCorr);
        set(gca, 'XTick', 1:height(roi_data), 'XTickLabel', roi_data.Subject, 'XTickLabelRotation', 45);
        ylabel('Median Corr'); title(sprintf('ROI%d %s: Median', roi_id, upper(hemi_str{1}))); grid on;
        
        subplot(2,2,3); bar(roi_data.Dice);
        set(gca, 'XTick', 1:height(roi_data), 'XTickLabel', roi_data.Subject, 'XTickLabelRotation', 45);
        ylabel('Dice'); title(sprintf('ROI%d %s: Dice', roi_id, upper(hemi_str{1}))); grid on;
        
        subplot(2,2,4); bar(roi_data.OverlapPct);
        set(gca, 'XTick', 1:height(roi_data), 'XTickLabel', roi_data.Subject, 'XTickLabelRotation', 45);
        ylabel('Overlap %'); title(sprintf('ROI%d %s: Overlap', roi_id, upper(hemi_str{1}))); grid on;
        
        sgtitle(sprintf('ROI%d %s', roi_id, upper(hemi_str{1})), 'FontSize', 14, 'FontWeight', 'bold');
    end
end

%% Group average plots
for hemi_str = {'lh', 'rh'}
    figure('Name', sprintf('Group Average - %s', upper(hemi_str{1})), 'Position', [100 100 1000 600]);
    
    means = []; sems = []; ns = []; ids = [];
    for roi_id = target_rois
        data = ResultsTable(ResultsTable.ROI == roi_id & strcmp(ResultsTable.Hemisphere, hemi_str{1}), :);
        if isempty(data), continue; end
        means = [means; mean(data.OverlapPct)];
        sems = [sems; std(data.OverlapPct)/sqrt(height(data))];
        ns = [ns; height(data)];
        ids = [ids; roi_id];
    end
    
    if isempty(means), continue; end
    
    bar(means); hold on;
    errorbar(1:length(means), means, sems, 'k.', 'LineWidth', 1.5, 'MarkerSize', 15);
    set(gca, 'XTick', 1:length(ids), 'XTickLabel', arrayfun(@(x) sprintf('ROI%d', x), ids, 'UniformOutput', false));
    
    for i = 1:length(means)
        text(i, means(i)+sems(i)+2, sprintf('N=%d', ns(i)), 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    
    ylabel('Overlap % (Mean ± SEM)', 'FontSize', 12); xlabel('ROI', 'FontSize', 12);
    title(sprintf('%s: Group Average Overlap (95th %%ile)', upper(hemi_str{1})), 'FontWeight', 'bold');
    ylim([0 max(means+sems)*1.3]); grid on; box on;
end