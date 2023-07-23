SELECT COUNT(*)
FROM model_param t1
         JOIN analytical_output t2 ON t1.model_param_id = t2.model_param_id
WHERE t1.circuit_n = '14'
  and t1."variant" = '2nd';


select count(*)
from model_param p,
     analytical_output o
where p.model_param_id = o.model_param_id
  and p.circuit_n = '14'
  and p."variant" = '2nd';

select count(*)
from model_param p,
     simulation_param o
where p.model_param_id = o.model_param_id
  and p.circuit_n = 'turinghill'
  and p."variant" = '9';


-- two tables
SELECT count(*)
FROM simulation_output o
         JOIN model_param m on o.model_param_id = m.model_param_id
         JOIN simulation_param s ON o.simulation_param_id = s.simulation_param_id
-- WHERE  m.circuit_n = '14' and m.variant='2nd' and m.n_samples=1000000
-- AND s."T" = 2000 and s.growth='nogrowth' and "boundaryCoeff"=1;
-- AND s."shape" ='ca';
where s."shape" = 'ca';


select p."model_param_id"
from model_param p,
     analytical_output a
where a."ss_n" = 1
limit 10;
select p."model_param_id"
from model_param p,
     analytical_output a
where p.model_param_id = a.model_param_id
  AND p."circuit_n" = '14'
  AND p."variant" = '2nd'
  AND p."n_samples" = 1000000
  AND p."balance" = 'Balanced'


with params_ss1 as (select p."model_param_id"
                    from model_param p,
                         analytical_output a
                    where p.model_param_id = a.model_param_id
                      AND a.ss_n = 1
                      AND p."circuit_n" = '14'
                      AND p."variant" = '2nd'
                      AND p."n_samples" = 1000000
                      AND p."balance" = 'Balanced'),
     sim_params as (select "simulation_param_id" from simulation_param where shape = 'ca' and "T" = 25)
select o."U_final_1D"
from simulation_output o,
     params_ss1,
     sim_params
where o.model_param_id = params_ss1.model_param_id
  and o.simulation_param_id = sim_params.simulation_param_id
;

select count(*) from simulation_output o  where o."U_final_1D" is not null;



with params_ss1 as (select p."model_param_id" from model_param p, analytical_output a where p.model_param_id = a.model_param_id
--                       AND a.ss_n = 1
                      AND p."circuit_n" = '14'
                      AND p."variant" = '2nd'
                      AND p."n_samples" = 1000000
                      AND p."balance" = 'Balanced'),
     sim_params as (select "simulation_param_id" from simulation_param where shape='ca' and "T"=25)
--      notnull_sim_output as (select * from simulation_output where "U_final_1D" is not null)


-- select count(*)
select  o."model_param_id"
from simulation_output o,
     params_ss1,
     sim_params , model_param p
where o.model_param_id = params_ss1.model_param_id
  and o.simulation_param_id = sim_params.simulation_param_id
;

select count("U_final_1D") from simulation_output where model_param_id='0_circuit:14_variant:2nd_samples:1000000';
select model_param_id from simulation_output limit 1;
select model_param_id from simulation_output where "U_final_1D" is not null;
select count(*) from simulation_output where  model_param_id='32211_circuit:14_variant:2nd_samples:1000000';