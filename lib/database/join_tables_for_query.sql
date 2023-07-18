SELECT COUNT(*)
FROM model_param t1
JOIN analytical_output t2 ON t1.model_param_id = t2.model_param_id
WHERE t1.circuit_n = '14' and t1."variant"='2nd';


select count(*) from model_param p, analytical_output o
    where p.model_param_id = o.model_param_id and
          p.circuit_n = '14' and
          p."variant"='2nd';

select count(*) from model_param p, simulation_param o
    where p.model_param_id = o.model_param_id and
          p.circuit_n = 'turinghill' and
          p."variant"='9';


-- two tables
SELECT count(*)
FROM simulation_output o
JOIN model_param m on o.model_param_id = m.model_param_id
JOIN simulation_param s ON o.simulation_param_id = s.simulation_param_id
-- WHERE  m.circuit_n = '14' and m.variant='2nd' and m.n_samples=1000000
-- AND s."T" = 2000 and s.growth='nogrowth' and "boundaryCoeff"=1;
-- AND s."shape" ='ca';
where s."shape" ='ca';

