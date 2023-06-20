-- define foreign key in analytical output pointing to model_param

alter table analytical_output
    add "model_param_id" text NOT NULL default '';

update analytical_output
set "model_param_id" = p."model_param_id"

from model_param p
where p.n_samples = analytical_output.n_samples
  and p.variant = analytical_output.variant
  and p.circuit_n = analytical_output.circuit_n
  and p."parID" = analytical_output."parID";

alter table model_param
    add constraint model_param_id_unique unique ("model_param_id");

alter table analytical_output
    add constraint fk_model_param
        foreign key (model_param_id) references model_param (model_param_id);


alter table analytical_output
drop constraint analytical_output_pkey;

ALTER TABLE analytical_output ADD PRIMARY KEY (model_param_id, "ssID");

alter table analytical_output
drop column "parID",
drop column "circuit_n",
drop column "variant",
drop column "n_samples";