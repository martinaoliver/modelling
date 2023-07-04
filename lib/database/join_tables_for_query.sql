SELECT COUNT(*)
FROM model_param t1
JOIN analytical_output t2 ON t1.model_param_id = t2.model_param_id
WHERE t1.circuit_n = 'turinghill';