select * from analytical_output ao
join model_param mp on mp.model_param_id = ao.model_param_id
where  "parID"=6;
-- and ao.system_class='no steady state';

SELECT system_class, COUNT(*) AS count
FROM analytical_output ao
join model_param mp on ao.model_param_id = mp.model_param_id
where mp.variant='8' or mp.variant='9'
    GROUP BY system_class
;

