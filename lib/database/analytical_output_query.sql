select * from analytical_output ao
join model_param mp on mp.model_param_id = ao.model_param_id
where  "parID"=6;
-- and ao.system_class='no steady state';

SELECT system_class, COUNT(*) AS count
FROM analytical_output ao
join model_param mp on ao.model_param_id = mp.model_param_id
where circuit_n='turinghill'
and mp.variant='12'  and
 mp.n_samples=2000000
    GROUP BY system_class
;

select * from analytical_output ao
join model_param mp on mp.model_param_id = ao.model_param_id
where  mp."parID"=388372;

select * from analytical_output ao
join model_param mp on mp.model_param_id = ao.model_param_id
where mp.circuit_n='turinghill' and mp.variant = '8' and mp.n_samples=1000000 and mp."parID"=388372;

select * from (SELECT "ssID", "ss_list", "model_param_id"
               FROM analytical_output
               GROUP BY "ssID", "ss_list", "model_param_id"
               HAVING COUNT("ssID") > 1
) as sub
    --as sub where sub.ss_list = '{null}'

select * from model_param mp
join analytical_output ao on mp.model_param_id = ao.model_param_id
-- where "parID" = 41 and mp.circuit_n='14';
where (mp.model_param_id = '55946_circuit:turinghill_variant:9_samples:2000000' or mp.model_param_id ='4892935_circuit:turinghill_variant:9_samples:1000000')

select count(*) from analytical_output
join model_param mp on mp.model_param_id = analytical_output.model_param_id
where circuit_n='turinghill'
and variant='12'
and n_samples=2000000;
-- and system_class='simple stable';
