select *
from simulation_output
         join model_param mp on mp.model_param_id = simulation_output.model_param_id
         join analytical_output ao on mp.model_param_id = ao.model_param_id
         join simulation_param sp on simulation_output.simulation_param_uuid = sp.simulation_param_uuid
where mp."parID" = 5912403;
-- and ao.system_class='no steady state';


-- gives you counts of each system class for numerically solved parameters
SELECT ao.system_class, COUNT(*) AS count
FROM simulation_output so
         JOIN model_param mp ON mp.model_param_id = so.model_param_id
         JOIN analytical_output ao ON mp.model_param_id = ao.model_param_id
         JOIN simulation_param sp ON so.simulation_param_uuid = sp.simulation_param_uuid
WHERE mp.variant = '9'
GROUP BY ao.system_class;



select *
from simulation_output
where false;


select *
from simulation_output
         join model_param mp on mp.model_param_id = simulation_output.model_param_id
where mp."parID" = 5912403;


select *
from simulation_param
where simulation_param_uuid = 'f557b922-67b0-4d93-aad1-4a6c362240c9';

select ao.ss_n, mp.variant,simulation_param_uuid,  ao.model_param_id,ao.system_class,pco.pattern_class_nogrowth from pattern_class_output pco
join analytical_output ao on pco.model_param_id = ao.model_param_id
join model_param mp on mp.model_param_id = ao.model_param_id
where simulation_param_uuid = 'f557b922-67b0-4d93-aad1-4a6c362240c9'
          and ss_n=1
and( mp.variant='9' or mp.variant='8');



select ao.ss_n,
       mp.variant,
       simulation_param_uuid,
       ao.model_param_id,
       ao.system_class,
       pco.pattern_class_nogrowth,
       pco.pattern_class_openboundary,
       pco.pattern_class_edgegrowth2
from pattern_class_output pco
         join analytical_output ao on pco.model_param_id = ao.model_param_id
         join model_param mp on mp.model_param_id = ao.model_param_id
where simulation_param_uuid = 'f557b922-67b0-4d93-aad1-4a6c362240c9'
  and ss_n = 1
  and (mp.variant = '9' or mp.variant = '8');

SELECT pco.model_param_id,
       MAX(pco.pattern_class_nogrowth)     AS pattern_class_nogrowth,
       MAX(pco.pattern_class_openboundary) AS pattern_class_openboundary,
       MAX(pco.pattern_class_edgegrowth2)  AS pattern_class_edgegrowth2

FROM pattern_class_output pco
join model_param mp on mp.model_param_id = pco.model_param_id
join analytical_output ao on mp.model_param_id = ao.model_param_id
where ( mp.variant='9' or mp.variant='8')
and pco.simulation_param_uuid='f557b922-67b0-4d93-aad1-4a6c362240c9'
and ao.ss_n=1

GROUP BY pco.model_param_id ;

SELECT
            model_param_id,
            MAX(pattern_class_nogrowth) AS pattern_class_nogrowth,
            MAX(pattern_class_openboundary) AS pattern_class_openboundary,
            MAX(pattern_class_edgegrowth2) AS pattern_class_edgegrowth2

            FROM pattern_class_output


            GROUP BY model_param_id;


with cluster_pattern_class as (
SELECT
            model_param_id,
            MAX(pattern_class_nogrowth) AS pattern_class_nogrowth,
            MAX(pattern_class_openboundary) AS pattern_class_openboundary,
            MAX(pattern_class_edgegrowth2) AS pattern_class_edgegrowth2

            FROM pattern_class_output


            GROUP BY model_param_id)


select pco.model_param_id, pattern_class_nogrowth, pattern_class_openboundary, pattern_class_edgegrowth2 from cluster_pattern_class pco

join model_param mp on mp.model_param_id = pco.model_param_id
join analytical_output ao on mp.model_param_id = ao.model_param_id
where ( mp.variant='9' or mp.variant='8')
and ao.ss_n=1;



with lsa_vs_numerical as
    (select ao.ss_n, mp.variant,simulation_param_uuid,  ao.model_param_id,ao.system_class,pco.pattern_class_nogrowth from pattern_class_output pco
join analytical_output ao on pco.model_param_id = ao.model_param_id
join model_param mp on mp.model_param_id = ao.model_param_id
where simulation_param_uuid = 'f557b922-67b0-4d93-aad1-4a6c362240c9'
--           and ss_n=1
and( mp.variant='9' or mp.variant='8')
and ao.system_class='complex unstable')
select  * from lsa_vs_numerical

