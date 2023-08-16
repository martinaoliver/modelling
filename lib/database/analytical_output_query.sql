select count(*) from model_param mp
-- join model_param mp on mp.model_param_id = ao.model_param_id
where mp.variant='fitted7';-- and "parID"=12837401;
-- and ao.system_class='no steady state';


