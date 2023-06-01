function [load_table,data_db] = parse_db(db_dir)
% Tara Urner, Buckley Lab 2020
% load_table: the table of files as loaded
% patient_db: a table parsed on patient id and session id
% structure:
% patient_db{:,1} = patient id
% patient_db{:,2} = locations of patient data in load_table
% patient_db{:,3} = sessions data
% patient_db{:,3}{:,1} = session id
% patient_db{:,3}{:,2} = location in load table
% patient_db{:,3}{:,3} = filepaths
% patient_db{:,3}{:,3}.data_type = type of data contained in .csv
% patient_db{:,3}{:,3}.folder = folder location of data
% patient_db{:,3}{:,3}.name = name of file
% patient_db{:,3}{:,3}.path = full path for loading.
    cd(db_dir)
    pdb_dir = dir('*mat');
    pdb_t = struct2table(pdb_dir);
    fs = pdb_t.name;
    a = string(fs);
    %%
    % Find all patient/session names
    split_names = arrayfun( @(x) strsplit( x, '--' ), a, 'UniformOutput', false );
    names = vertcat(split_names{:});
    subject_idx = convertStringsToChars(names(:,1));
    pdb_t.subject_idx = subject_idx;
    pdb_t.session_idx = subject_idx;
    data_type = cellfun(@(x) strsplit(x,'.'),cellfun(@(x) x{2},split_names,'UniformOutput',false),'UniformOutput',false);
    pdb_t.data_type = cellfun(@(x) x{1},data_type,'UniformOutput',false);
    unique_sidx = unique(pdb_t.subject_idx);
    for i=1:length(unique_sidx)
        unique_sidx{i,2} = cell2mat(cellfun(@(x) any(strfind(x,unique_sidx{i})),fs,'UniformOutput',false));
        sessions = table2cell(unique(pdb_t(unique_sidx{i,2},'session_idx')));
        %session_id = pdb_t.session_idx;
        for j=1:size(sessions,1)
            % Unique sessions
            sessions{j,2} = cell2mat(cellfun(@(x) any(strfind(x,sessions{j})),pdb_t.session_idx,'UniformOutput',false));
            session_info = string(table2cell(pdb_t(sessions{j,2},'data_type')));
            session_info(:,2) = string(table2cell(pdb_t(sessions{j,2},'folder')));
            session_info(:,3) = string(table2cell(pdb_t(sessions{j,2},'name')));
            session_info(:,4) = strcat(session_info(:,2),filesep,session_info(:,3));
            session_info = array2table(session_info,'VariableNames',{'data_type','folder','name','path'});
            sessions{j,3} = session_info;
        end
        unique_sidx{i,3} = sessions;
        load_table = pdb_t;
        data_db = unique_sidx;
    end

end