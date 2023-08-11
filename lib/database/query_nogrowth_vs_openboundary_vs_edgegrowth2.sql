SELECT
    model_param_id,
    MAX(pattern_class_nogrowth) AS pattern_class_nogrowth,
    MAX(pattern_class_openboundary) AS pattern_class_openboundary,
    MAX(pattern_class_edgegrowth2) AS pattern_class_edgegrowth2
FROM pattern_class_output
GROUP BY model_param_id;

