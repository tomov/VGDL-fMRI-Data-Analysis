function level = get_level_from_run_and_instance(run_id, instance_id)
    level = 1 + instance_id + floor((run_id - 1) / 2) * 3;
