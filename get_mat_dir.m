function dir = get_mat_dir(fasse_ncf)

    if ~exist('fasse_ncf', 'var')
        fasse_ncf = false;
    end

    [~, name] = system('hostname');
    if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
        % local
        dir = 'mat/';
    elseif  ~isempty( strfind(name,'fasse')) || ~isempty( strfind(name,'holy'))
        % fasse
        if fasse_ncf
            % ncf Mount
            dir = fullfile(getenv('MY_NCF_LAB'), 'Lab/scripts/matlab/VGDL_fMRI/mat/');
        else
            % regular
            dir = fullfile(getenv('MY_LAB'), 'VGDL', 'mat/');
        end
    else
        % ncf
        dir = 'mat/';
    end
    
