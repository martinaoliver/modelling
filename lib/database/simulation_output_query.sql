select * from simulation_output
join model_param mp on mp.model_param_id = simulation_output.model_param_id
join analytical_output ao on mp.model_param_id = ao.model_param_id
join simulation_param sp on simulation_output.simulation_param_uuid = sp.simulation_param_uuid
where mp."parID"=5912403 ;
-- and ao.system_class='no steady state';


-- gives you counts of each system class for numerically solved parameters
SELECT ao.system_class, COUNT(*) AS count
FROM simulation_output so
JOIN model_param mp ON mp.model_param_id = so.model_param_id
JOIN analytical_output ao ON mp.model_param_id = ao.model_param_id
JOIN simulation_param sp ON so.simulation_param_uuid = sp.simulation_param_uuid
WHERE mp.variant = '9'
GROUP BY ao.system_class;



select * from simulation_output where false;


select * from simulation_output
join model_param mp on mp.model_param_id = simulation_output.model_param_id
where mp."parID"=5912403 ;


