-- this adds a 'generated' column (cannot be inserted) that exists for every row, which is a function of other
-- columns in the same row :)
alter table model_param
    add column "model_param_id" text
    GENERATED ALWAYS AS (("parID"::text || '_circuit:' || "circuit_n" || '_variant:' || "variant" || '_samples:' || "n_samples"::text)::text) STORED ;