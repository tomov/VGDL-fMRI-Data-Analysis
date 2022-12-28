function dir = get_mat_dir(location)

    if ~exist('location', 'var')
        location = 0; % fasse storage
    end

    [~, name] = system('hostname');
    if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
        % local
        dir = 'mat/';
    elseif  ~isempty( strfind(name,'fasse')) || ~isempty( strfind(name,'holy'))
        % fasse
        switch location
            case 1
                % ncf Mount
                dir = fullfile(getenv('MY_NCF_LAB'), 'Lab/scripts/matlab/VGDL_fMRI/mat/');
            case 0
                % regular storage
                %dir = fullfile(getenv('MY_LAB'), 'VGDL', 'mat/');
                dir = fullfile(getenv('MY_SCRATCH'), 'VGDL', 'mat_from_lab/');
            case 2
                %  fasse scratch
                dir = fullfile(getenv('MY_SCRATCH'), 'VGDL', 'mat/');
            otherwise
                assert(false);
        end
    else
        % ncf
        dir = 'mat/';
    end
    
