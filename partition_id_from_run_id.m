function partition_id = partition_id_from_run_id(run_id)
    % convenience function for getting partitions for CV from run_ids

    partition_id = floor((run_id - 1) / 2) + 1);
