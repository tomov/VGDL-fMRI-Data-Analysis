function dir = get_mat_dir()


    [~, name] = system('hostname');
    if ~isempty(strfind(name, 'omchil')) || ~isempty(strfind(name, 'dhcp-'))
        % local
        dir = 'mat/';
    elseif  ~isempty( strfind(name,'fasse')) || ~isempty( strfind(name,'holy'))
        % fasse
        dir = fullfile(getenv('MY_LAB'), 'VGDL', 'mat/');
    else
        % ncf
        dir = 'mat/';
    end
    
