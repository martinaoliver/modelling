select ao.model_param_id, ao."ssID", sp.simulation_param_uuid, ao."ssID", sp."T" from simulation_output
join model_param mp on mp.model_param_id = simulation_output.model_param_id
join analytical_output ao on mp.model_param_id = ao.model_param_id
join simulation_param sp on simulation_output.simulation_param_uuid = sp.simulation_param_uuid
where mp.variant='fitted7_gaussian4187715_nsr0.01' and sp."T"='110';
-- and ao.system_class='no steady state';


