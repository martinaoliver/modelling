select ao.model_param_id,ao.system_class,pco.pattern_class_nogrowth from pattern_class_output pco
join analytical_output ao on pco.model_param_id = ao.model_param_id
join model_param mp on mp.model_param_id = ao.model_param_id
where mp.variant='9'
and simulation_param_uuid = 'f557b922-67b0-4d93-aad1-4a6c362240c9';