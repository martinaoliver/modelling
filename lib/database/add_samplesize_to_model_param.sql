alter table model_param drop constraint model_param_pkey;

ALTER TABLE model_param
  ADD "n_samples" numeric NOT NULL default 5000;

ALTER TABLE model_param ADD PRIMARY KEY ("parID", "circuit_n", "variant", "n_samples");

ALTER table model_param
ALTER column "n_samples" DROP default;


