select model_param.*, ao.system_class from model_param
         join analytical_output ao on model_param.model_param_id = ao.model_param_id
where "variant"='2nd'
and "parID" =598822;


--query in order
-- SELECT mp.*,ao.system_class, ao.ss_n
SELECT "parID", ao.system_class,"Dr", "Va", "Vb", "Vc", "Vd","Ve","Vf","Kub","Kvd","Kda","Kce","Kfe","Keb","Kee","muLVA","muASV","nub","nee","neb","nvd","nda","nce","nfe"
FROM model_param mp
-- JOIN unnest(ARRAY[356642,832446,598822,983486,963665,361095,162101,739015,555549,141318,278021,311062,146017,376154,252,630,434,93,549,595,100232,238,596]) WITH ORDINALITY AS idx(id, custom_order) ON mp."parID" = idx.id
-- JOIN unnest(ARRAY[252,356642,598822,983486,141318,361095]) WITH ORDINALITY AS idx(id, custom_order) ON mp."parID" = idx.id
JOIN unnest(ARRAY[555549, 141318, 739015, 311062,278021,146017,434,549, 93, 252, 630, 595,100232, 596, 238]) WITH ORDINALITY AS idx(id, custom_order) ON mp."parID" = idx.id
join analytical_output ao on mp.model_param_id = ao.model_param_id
where "variant"='2nd'
ORDER BY idx.custom_order;

SELECT mp.*
FROM model_param mp
JOIN unnest(ARRAY[6]) WITH ORDINALITY AS idx(id, custom_order) ON mp."parID" = idx.id
where "variant"='fitted7_gaussian4187715_nsr0.01'
ORDER BY idx.custom_order;

