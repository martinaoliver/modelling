select * from analytical_output ao
join model_param mp on mp.model_param_id = ao.model_param_id
where  "parID"=6;
-- and ao.system_class='no steady state';

SELECT system_class, COUNT(*) AS count
FROM analytical_output ao
join model_param mp on ao.model_param_id = mp.model_param_id
where (mp.variant='8' or mp.variant='9') and
 mp.n_samples=1000000
    GROUP BY system_class
;

select * from analytical_output ao
join model_param mp on mp.model_param_id = ao.model_param_id
where  mp."parID"=5912403;

select * from analytical_output ao
join model_param mp on mp.model_param_id = ao.model_param_id
where mp.circuit_n='turinghill' and mp.variant = '8' and mp.n_samples=1000000;

select * from (SELECT "ssID", "ss_list", "model_param_id"
               FROM analytical_output
               GROUP BY "ssID", "ss_list", "model_param_id"
               HAVING COUNT("ssID") > 1
) as sub
    --as sub where sub.ss_list = '{null}'

