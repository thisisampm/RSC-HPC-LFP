function plotwaveforms( waveforms, varargin )
%PLOTWAVEFORMS plots wave forms for each cell in the row vector clusters
% if you'd rather plot features (i.e. from sorting), specify which features
% with a [1,2] or [1,3] vect of numbers 1-8, with 1:4 indicating peaks and 
% 5-8 indicating valleys. Clusters may be input as a third variable to 
% label the legend 
%
% if second input is the string 'auto', the features with the highest
% variability among cluster means will be plotted
%
% if second input is the string 'pca', then all features will be input into
% a pca and the first 3 principal components will be ploted

%check for max min
if nargin == 2
    features = varargin{1};
    clusters = (1:length(waveforms))';
elseif nargin == 3
    features = varargin{1};
    clusters = varargin{2};
else
    clusters = (1:length(waveforms))';
end


%check for features strings
if nargin>1
    if strcmp(features, 'auto')
        features = auto_features(waveforms)
    elseif strcmp(features, 'pca')
        pca_features(waveforms, clusters)
        return
    end
end

figure; hold on; 
for icluster = 1:length(waveforms)
    legend_chars{icluster} = num2str(clusters(icluster,1));
    
    if exist('features', 'var') %plot features
        
        %features
        clearvars features_all

        %peak or valley
        for ifeatures = 1:length(features)
            if features(ifeatures) <=4 && features(ifeatures) >=1
                features_all(:,ifeatures) = max(waveforms{icluster}(:,features(ifeatures),:));                
            elseif features(ifeatures) >=5 && features(ifeatures) <=8
                features_all(:,ifeatures) = min(waveforms{icluster}(:,features(ifeatures)-4,:));
            end
        end
        
        %2D or 3D plot
        if length(features) == 2
            plot(features_all(:,1), features_all(:,2), '.', 'color', ([1 1 1]./length(waveforms).*(icluster-1)));
        elseif length(features) == 3
            plot3(features_all(:,1), features_all(:,2), features_all(:,3), '.', 'color', ([1 1 1]./length(waveforms).*(icluster-1)));
        else
            error('incorrect number of features input')
        end

    else %plot full waveforms
        
        %for event = 1:size(waveforms{icluster},3) 
            for electrode = 1:4
                subplot(2,2, electrode) 
                hold on
                %waveform = waveforms{icluster}(:,electrode,event);
                waveform = mean(waveforms{icluster}(:,electrode,:), 3);
                plot(waveform, 'color', ([1 1 1]./length(waveforms).*(icluster-1))); 
            end
        %end 
    end
    
end

%aesthetics
if exist('features', 'var')
    set(gca, 'TickLength', [0 0])
    legend(legend_chars)
    axis square

    feature_chars = {'Max 1', 'Max 2', 'Max 3', 'Max 4',...
        'Min 1', 'Min 2', 'Min 3', 'Min 4'};

    xlabel(feature_chars(features(1)))
    ylabel(feature_chars(features(2)))
    if length(features) == 3
        zlabel(feature_chars(features(3))); 
    end
end
end



function features = auto_features(waveforms)

    %preallocate means
    feature_means = nan(length(waveforms), 8);
    
    %for each feature
    for ifeature = 1:8

        %for each waveform
        for iclust = 1:length(waveforms)
            
            %calculate mean of that feature
            if ifeature <=4
                feature_means(iclust,ifeature) = mean(max(waveforms{iclust}(:,ifeature,:)));                
            else
                feature_means(iclust,ifeature) = mean(min(waveforms{iclust}(:,ifeature-4,:)));
            end
        end
    end
    
    %find variances
    feature_vars = var(feature_means);
    
    %find best features
    [~,sortfeat] = sort(feature_vars, 'descend');
    features = sortfeat(1:3);
    
end

function pca_features(waveforms, clusters)

    %preallocate soft
    feature_elements_all = [];
    waveform_idx = [];
    
    %for each waveform
    for iclust = 1:length(waveforms)
        legend_chars{iclust} = num2str(clusters(iclust,1));
        
        %preallocate soft
        feature_elements_local = nan(size(waveforms{iclust},3), 8);
    
        %for each feature
        for ifeature = 1:8

            %calculate that feature
            if ifeature <=4
                feature_elements_local(:, ifeature) = squeeze(max(waveforms{iclust}(:,ifeature,:))); 
            else
                feature_elements_local(:, ifeature) = squeeze(min(waveforms{iclust}(:,ifeature-4,:)));
            end
            
        end
        
        %load
        feature_elements_all = [feature_elements_all; feature_elements_local];
        waveform_idx = [waveform_idx; repmat(iclust, size(feature_elements_local,1), 1)];
    end
    
    feature_elements_all_size = size(feature_elements_all);
    downsamp_size = 5000; %larger vectors take too long to compute
    downsamp = round(linspace(1, size(feature_elements_all,1), downsamp_size));
    
    
    feature_elements_all = feature_elements_all(downsamp, :);
    feature_elements_all = zscore_mtx(feature_elements_all);
    waveform_idx = waveform_idx(downsamp, :);
    
    
    feature_elements_all = feature_elements_all(:,[1 2 3 5 6 7]);
    
    %pca
    [feature_elements_pca] = pca(feature_elements_all');

    %plot
    pca_axes = [1 2 3];
    figure; hold on
    for iclust = 1:length(waveforms)
        wf_idx = waveform_idx == iclust;
        %plot3(feature_elements_pca(wf_idx,pca_axes(1)), feature_elements_pca(wf_idx,pca_axes(2)),...
            %feature_elements_pca(wf_idx,pca_axes(3)), '.', 'color', ([1 1 1]./length(waveforms).*(iclust-1)))
        plot3(feature_elements_pca(wf_idx,pca_axes(1)), feature_elements_pca(wf_idx,pca_axes(2)),...
            feature_elements_pca(wf_idx,pca_axes(3)), '.')
    end
    
    %asthetics
    set(gca, 'TickLength', [0 0])
    legend(legend_chars)
    xlabel('PC 1')
    ylabel('PC 2')
    zlabel('PC 3');
    axis square
    
end

