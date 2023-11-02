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
-- and ao.system_class='complex unstable'
and mp.n_samples=1000000
and mp."parID"=388372)
select  * from lsa_vs_numerical;

with lsa_vs_numerical as
    (select ao.ss_n, mp.variant,simulation_param_uuid,  ao.model_param_id,ao.system_class,pco.pattern_class_nogrowth from pattern_class_output pco
join analytical_output ao on pco.model_param_id = ao.model_param_id
join model_param mp on mp.model_param_id = ao.model_param_id
where simulation_param_uuid = '8627037d-356d-4489-9215-1cab9c82638a'
          and ss_n=1
and( mp.variant='9' or mp.variant='8')
-- and ao.system_class='complex unstable'
and mp.n_samples=1000000)
select lsa_vs_numerical.system_class, count(*) as count from lsa_vs_numerical
group by lsa_vs_numerical.system_class;



select ao.system_class, count(*)
from simulation_output so
join model_param mp on mp.model_param_id = so.model_param_id
join analytical_output ao on so.model_param_id = ao.model_param_id
where simulation_param_uuid = '8627037d-356d-4489-9215-1cab9c82638a'
          and ss_n=1
and( mp.variant='9' or mp.variant='8')
-- and ao.system_class='complex unstable'
and mp.n_samples=1000000
group by ao.system_class ;


select ao.system_class,  mp."parID", pco.pattern_class_nogrowth
from pattern_class_output pco
join model_param mp on mp.model_param_id = pco.model_param_id
join analytical_output ao on ao.model_param_id = pco.model_param_id
-- where pco.simulation_param_uuid = '8627037d-356d-4489-9215-1cab9c82638a' --edgegrowth2
-- where pco.simulation_param_uuid = 'db183b7c-28dc-426f-9d3a-b4c7e6b3f260' --openboundary
where pco.simulation_param_uuid = 'f557b922-67b0-4d93-aad1-4a6c362240c9' --nogrowth

          and ss_n=1
and( mp.variant='9' or mp.variant='8')
-- and ao.system_class='complex unstable'
and mp.n_samples=1000000
and mp.circuit_n='turinghill'
  and pco.pattern_class_nogrowth='Homogeneous'
and ao.system_class='simple stable'
 ;


-- select mp."parID", so."ssID"  from simulation_output so
select count(*) from simulation_output so
inner join model_param mp on so.model_param_id = mp.model_param_id
-- inner join analytical_output ao on (ao.model_param_id,ao."ssID") = (so.model_param_id, so."ssID")

-- where ao.system_class in ('turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf', 'turing semi-hopf')

where mp.variant='12'
-- where simulation_param_uuid='b94c9e61-a717-4470-957b-a59ff727e948'
and simulation_param_uuid='132323a4-3f93-4287-aca9-d18e84848e37'
;




-- select mp."parID", so."ssID"  from simulation_output so
select count(*) from simulation_output so
-- select count(*) from pattern_class_output pco
-- select mp."parID", ao."ssID", ao.system_class, pco.wavelength, ao.estimated_wvl, pco.convergence_time, ao.maxeig from pattern_class_output pco
-- select * from pattern_class_output pco
inner join model_param mp on pco.model_param_id = mp.model_param_id
inner join analytical_output ao on (ao.model_param_id,ao."ssID") = (pco.model_param_id, pco."ssID")

where ao.system_class not in ('turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf', 'turing semi-hopf')

and mp.variant='11'
-- where simulation_param_uuid='b94c9e61-a717-4470-957b-a59ff727e948'
and simulation_param_uuid='132323a4-3f93-4287-aca9-d18e84848e37';



-- select mp."parID", so."ssID"  from simulation_output so
select count(*) from simulation_output so
inner join model_param mp on so.model_param_id = mp.model_param_id
inner join analytical_output ao on (ao.model_param_id,ao."ssID") = (so.model_param_id, so."ssID")
--
-- where ao.system_class in ('turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf', 'turing semi-hopf')
-- where ao.system_class in ('hopf')
and mp.variant='11'
and so.simulation_param_uuid='132323a4-3f93-4287-aca9-d18e84848e37'
and mp.n_samples=1000000;




-- select mp."parID", pco."ssID"  from pattern_class_output pco
select count(*) from pattern_class_output pco
inner join model_param mp on pco.model_param_id = mp.model_param_id
inner join analytical_output ao on (pco.model_param_id,pco."ssID") = (ao.model_param_id, ao."ssID")
-- where pattern_class_nogrowth = 'Temporal Oscillator'
-- and mp."variant"='12'
and ao.system_class='hopf';

SELECT DISTINCT pco.pattern_class_nogrowth
FROM pattern_class_output pco
INNER JOIN model_param mp ON pco.model_param_id = mp.model_param_id
INNER JOIN analytical_output ao ON (pco.model_param_id, pco."ssID") = (ao.model_param_id, ao."ssID")
WHERE ao.system_class = 'hopf';